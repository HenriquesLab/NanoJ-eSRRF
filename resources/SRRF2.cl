//#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define MAX_FRAMES $MAX_FRAMES$
#define FWHM $FWHM$
#define RADIUS $RADIUS$
#define SIGMA $SIGMA$
#define MAGNIFICATION $MAGNIFICATION$
#define NTIMELAGS $NTIMELAGS$
#define DOBIN2 $DOBIN2$
#define DOBIN4 $DOBIN4$
#define w $WIDTH$
#define h $HEIGHT$
#define wh $WH$
#define wM $WM$
#define hM $HM$
#define whM $WHM$

// cubic function for interpolation
static float cubic(float x) {
    if (x < 0.0f) x = -x; // don't understand why? - do fabs instead?
    float z = 0.0f;
    if (x < 1.0f)
        z = x * x * (x * 1.5f - 2.5f) + 1.0f;
    else if (x < 2.0f)
        z = -.5f * x * x * x + 2.5f * x * x - 4.0f * x + 2.0f;
    return z;
}

// interpolation function: interpolate in continuous space with respect to the reference of the array
static float getInterpolatedValue(__global float* array, int const width, int const height, float const x, float const y, float const f) { // TODO: review the grid position in the interpolation (seems offset)
    const int fwh = f * width * height;
    const int u0 = (int) floor(x);
    const int v0 = (int) floor(y);

    float q = 0.0f;
    if (u0 > 0 && u0 < w - 2 && v0 > 0 && v0 < h - 2) { // do bicubic interpolation
        for (int j = 0; j <= 3; j++) {
            int v = v0 - 1 + j;
            float p = 0.0f;
            for (int i = 0; i <= 3; i++) {
                int u = u0 - 1 + i;
                p += array[fwh + v * width + u] * cubic(x - (float)u);
            }
            q = q + p * cubic(y - (float)v);
        }
    }
    else if (x >= 0 && x < w && y >= 0 && y < h) {
        int xbase = (int)x;
        int ybase = (int)y;
        int xbase1 = min(xbase+1, w-1);
        int ybase1 = min(ybase+1, h-1);

        float xFraction = x - (float)xbase;
        float yFraction = y - (float)ybase;
        xFraction = fmax(xFraction, 0);
        yFraction = fmax(yFraction, 0);

        float lowerLeft = array[fwh + ybase * width + xbase];
        float lowerRight = array[fwh + ybase * width + xbase1];
        float upperRight = array[fwh + ybase1 * width + xbase1];
        float upperLeft = array[fwh + ybase1 * width + xbase];
        float upperAverage = upperLeft + xFraction * (upperRight - upperLeft);
        float lowerAverage = lowerLeft + xFraction * (lowerRight - lowerLeft);
        q = lowerAverage + yFraction * (upperAverage - lowerAverage);
    }
    return q;
}

static float product(float* array, int start, int end) {
    float v = 1;
    for (int p=start; p<=end; p++) v *= array[p];
    return v;
}

// First kernel: evaluating gradient from image using Robert's cross ------------------------------------------
__kernel void calculateGradientRobX(
    __global float* pixels,
    __global float* GxArray,
    __global float* GyArray,
    int const nFrames
    ) {
    int x1 = get_global_id(0);
    int y1 = get_global_id(1);
    int xyOffset = y1 * w + x1;

    int x0 = max(x1-1, 0);
    int y0 = max(y1-1, 0);

    for (int f=0; f<nFrames; f++) {
        int fOffset = f * wh;

        int x0y0 = y0 * w + x0 + fOffset;
        int x0y1 = y1 * w + x0 + fOffset;
        int x1y0 = y0 * w + x1 + fOffset;
        int x1y1 = y1 * w + x1 + fOffset;

        // This calculates Robert's cross gradient and apply the rotation matrix 45 degrees to realign Gx and Gy to the image grid
        GxArray[fOffset+ xyOffset] = pixels[x1y0] - pixels[x0y1] + pixels[x1y1] - pixels[x0y0];
        GyArray[fOffset+ xyOffset] = - pixels[x1y0] + pixels[x0y1] + pixels[x1y1] - pixels[x0y0];
    }
}

__kernel void calculateSRRF(
    __global float* pixels,
    __global float* GxArray,
    __global float* GyArray,
    __global float* ShiftXArray,
    __global float* ShiftYArray,
    __global float* SRRFArray,
    int const nFrames
    ) {

    int xM = get_global_id(0);
    int yM = get_global_id(1);
    int xyMOffset = yM * wM + xM;

    float CGLH[MAX_FRAMES]; // note, MAX_FRAMES is passed before compile
    float vRAW_AVE = 0;
    float vSRRF_AVE = 0;
    float vSRRF_MAX = 0;

    // FIRST CALCULATE LOCAL GRADIENT CONVERGENCE (OLD RADIALITY!!)

    float sigma22 = 2 * SIGMA * SIGMA;
    float Gx, Gy;

    for (int f=0; f<nFrames; f++) {
        int fOffset = f * wh;
        float xc = (xM + 0.5f) / MAGNIFICATION + ShiftXArray[f] + 1; // continuous space position at the centre of magnified pixel
        float yc = (yM + 0.5f) / MAGNIFICATION + ShiftYArray[f] + 1; // this last +1 sum align the RGC to the interpolated magnified image - don't know why... but works

        float distanceWeightSum = 0;

        for (int j=-RADIUS; j<=RADIUS; j++) {
            for (int i=-RADIUS; i<=RADIUS; i++) {

                int vxPixelOrigin = (int) xc  + i;
                int vyPixelOrigin = (int) yc  + j;
                float vxPixelCentred = vxPixelOrigin + 0.5f;
                float vyPixelCentred = vyPixelOrigin + 0.5f;

                float dx = vxPixelCentred - xc;
                float dy = vyPixelCentred - yc;
                float distance = sqrt(dx * dx + dy * dy);

                if (distance != 0 && vxPixelCentred >=0 && vxPixelCentred < w && vyPixelCentred >=0 && vyPixelCentred < h) {
                    int p = fOffset + vyPixelOrigin * w + vxPixelOrigin;
                    Gx = GxArray[p];
                    Gy = GyArray[p];

                    float GMag = sqrt(Gx * Gx + Gy * Gy);

                    //float distanceWeight = exp(-0.5f*pow((float) (distanceWeight/SIGMA),2)); // gaussian weight
                    //float distanceWeight = exp(-0.5f*pow((float) (distanceWeight/(2*SIGMA)),4)); // gaussian flat top weight
                    float distanceWeight = distance*exp(-(distance*distance)/sigma22);  // TODO: dGauss: can use Taylor expansion there
                    distanceWeight = pow(distanceWeight, 2);

                    // Calculate perpendicular distance from (xc,yc) to gradient line through (vx,vy)
                    float Dk = fabs(Gy * (xc - vxPixelCentred) - Gx * (yc - vyPixelCentred)) / GMag;    // Dk = D*sin(theta)
                    if (isnan(Dk)) Dk = distance; // this makes Dk = 0 in the next line

                    Dk = 1 - Dk / distance; // Dk is now between 0 to 1, 1 if vector points precisely to (xc, yx), Dk = 1-sin(theta)
                    //Dk = fmax(Dk - 0.5f, 0)*2;

                    Dk *= distanceWeight;
                    distanceWeightSum += distanceWeight;

                    // Accumulate Variables
                    float GdotR = (Gx * i * MAGNIFICATION + Gy * j * MAGNIFICATION); // tells you if vector was pointing inward or outward
                    if (GdotR <= 0) CGLH[f] += Dk; // vector was pointing inwards
                    //else CGLH[f] -= Dk; // vector was pointing outwards
                }
            }
        }

        float v = getInterpolatedValue(pixels, w, h, ((float) xM)/MAGNIFICATION + ShiftXArray[f], ((float) yM)/MAGNIFICATION + ShiftYArray[f], f);
        CGLH[f] /= distanceWeightSum;
        //CGLH[f] = fmax(CGLH[f] - 0.2f, 0) * 1.2;

        // CALCULATE SRRF AVERAGE AND MAXIMUM
        int f1 = f+1;
        vRAW_AVE += (v - vRAW_AVE) / f1;
        vSRRF_MAX = fmax(CGLH[f], vSRRF_MAX);
        vSRRF_AVE += (CGLH[f] - vSRRF_AVE) / f1;
    }

    // CALCULATE SRRF BASIC TEMPORAL ANALYSIS
    SRRFArray[0 * whM + xyMOffset]  = vRAW_AVE;
    SRRFArray[1 * whM + xyMOffset]  = vSRRF_MAX * vRAW_AVE;
    SRRFArray[2 * whM + xyMOffset]  = vSRRF_AVE * vRAW_AVE;

    // CALCULATE SRRF TEMPORAL CORRELATIONS

    // mean subtract first
    for (int f=0; f<nFrames; f++) CGLH[f] = CGLH[f] - vSRRF_AVE;

    // initialize correlation array
    float SRRFCumTimeLags[NTIMELAGS]; // cumulative timelags
    for (int tl = 0; tl<NTIMELAGS; tl++) SRRFCumTimeLags[tl] = 0; // initialise SRRFCumTimeLags at 0

    ////////////////////////////////////////////////
    // CALCULATE SRRF TEMPORAL CORRELATIONS BIN 1 //
    ////////////////////////////////////////////////
    for (int f=0; f<nFrames; f++) {
        int f1 = f+1;
        SRRFCumTimeLags[0] += (CGLH[f] * CGLH[f] - SRRFCumTimeLags[0]) / f1;
        for (int n=1; n<NTIMELAGS; n++) {
            if (f < nFrames - n) SRRFCumTimeLags[n] += (fabs(product(CGLH, f, f+n)) - SRRFCumTimeLags[n]) / f1;
        }
    }

    int recOffset = 3;
    SRRFArray[recOffset * whM + xyMOffset]  = sqrt(SRRFCumTimeLags[0]) * vRAW_AVE;
    for (int n=1; n<NTIMELAGS; n++) {
        SRRFArray[(recOffset+n) * whM + xyMOffset] = rootn(SRRFCumTimeLags[n], 1+n) * vRAW_AVE;
    }

    ////////////////////////////////////////////////
    // CALCULATE SRRF TEMPORAL CORRELATIONS BIN 2 //
    ////////////////////////////////////////////////
    if (!DOBIN2) return; // if we don't want to calculate BIN 2, then return

    int nFrames2 = nFrames / 2;

    // bin data by 2 first
    for (int f=0; f<nFrames2; f++) CGLH[f] = CGLH[f*2] + CGLH[f*2+1];

    for (int tl = 0; tl<NTIMELAGS; tl++) SRRFCumTimeLags[tl] = 0; // initialise SRRFCumTimeLags at 0

    for (int f=0; f<nFrames2; f++) {
        int f1 = f+1;
        SRRFCumTimeLags[0] += (CGLH[f] * CGLH[f] - SRRFCumTimeLags[0]) / f1;
        for (int n=1; n<NTIMELAGS; n++) {
            if (f < nFrames - n) SRRFCumTimeLags[n] += (fabs(product(CGLH, f, f+n)) - SRRFCumTimeLags[n]) / f1;
        }
    }

    recOffset += NTIMELAGS;
    SRRFArray[(recOffset) * whM + xyMOffset]  = sqrt(SRRFCumTimeLags[0]) * vRAW_AVE;
    for (int n=1; n<NTIMELAGS; n++) {
        SRRFArray[(recOffset+n) * whM + xyMOffset] = rootn(SRRFCumTimeLags[n], 1+n) * vRAW_AVE;
    }

    ////////////////////////////////////////////////
    // CALCULATE SRRF TEMPORAL CORRELATIONS BIN 4 //
    ////////////////////////////////////////////////
    if (!DOBIN4) return; // if we don't want to calculate BIN 2, then return

    int nFrames4 = nFrames / 4;

    // bin data by 4 first
    for (int f=0; f<nFrames4; f++) CGLH[f] = CGLH[f*2] + CGLH[f*2+1];

    for (int tl = 0; tl<NTIMELAGS; tl++) SRRFCumTimeLags[tl] = 0; // initialise SRRFCumTimeLags at 0

    for (int f=0; f<nFrames4; f++) {
        int f1 = f+1;
        SRRFCumTimeLags[0] += (CGLH[f] * CGLH[f] - SRRFCumTimeLags[0]) / f1;
        for (int n=1; n<NTIMELAGS; n++) {
            if (f < nFrames - n) SRRFCumTimeLags[n] += (fabs(product(CGLH, f, f+n)) - SRRFCumTimeLags[n]) / f1;
        }
    }

    recOffset += NTIMELAGS;
    SRRFArray[recOffset * whM + xyMOffset]  = sqrt(SRRFCumTimeLags[0]) * vRAW_AVE;
    for (int n=1; n<NTIMELAGS; n++) {
        SRRFArray[(recOffset+n) * whM + xyMOffset] = rootn(SRRFCumTimeLags[n], 1+n) * vRAW_AVE;
    }
}

