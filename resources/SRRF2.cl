//#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define fwhm $FWHM$
#define sigma $SIGMA$
#define magnification $MAGNIFICATION$
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

    float CGLH[$MAX_FRAMES$]; // note, MAX_FRAMES is passed before compile
    float vRAW_AVE = 0;
    float vSRRF_AVE = 0;
    float vSRRF_MAX = 0;

    // FIRST CALCULATE LOCAL GRADIENT CONVERGENCE (OLD RADIALITY!!)

    float sigma22 = 2 * sigma * sigma;
    int radius = (int) (sigma * 2) + 1;    // radius can be set to something sensible like 3*Sigma
    float Gx, Gy;

    for (int f=0; f<nFrames; f++) {
        int fOffset = f * wh;
        float xc = (xM + 0.5f) / magnification + ShiftXArray[f] + 1; // continuous space position at the centre of magnified pixel
        float yc = (yM + 0.5f) / magnification + ShiftYArray[f] + 1; // this last +1 sum align the RGC to the interpolated magnified image - don't know why... but works

        float distanceWeightSum = 0;
        CGLH[f] = 0;

        for (int j=-radius; j<=radius; j++) {
            for (int i=-radius; i<=radius; i++) {

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

                    //float distanceWeight = exp(-0.5f*pow((float) (distanceWeight/sigma),2)); // gaussian weight
                    //float distanceWeight = exp(-0.5f*pow((float) (distanceWeight/(2*sigma)),4)); // gaussian flat top weight
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
                    float GdotR = (Gx * i * magnification + Gy * j * magnification); // tells you if vector was pointing inward or outward
                    if (GdotR <= 0) CGLH[f] += Dk; // vector was pointing inwards
                    //else CGLH[f] -= Dk; // vector was pointing outwards
                }
            }
        }

        float v = getInterpolatedValue(pixels, w, h, ((float) xM)/magnification + ShiftXArray[f], ((float) yM)/magnification + ShiftYArray[f], f);
        CGLH[f] /= distanceWeightSum;
        //CGLH[f] = fmax(CGLH[f] - 0.2f, 0) * 1.2;

        // CALCULATE SRRF AVERAGE AND MAXIMUM
        int f1 = f+1;
        vRAW_AVE += (v - vRAW_AVE) / f1;
        vSRRF_MAX = fmax(CGLH[f], vSRRF_MAX);
        vSRRF_AVE += (CGLH[f] - vSRRF_AVE) / f1;
    }

    // CALCULATE SRRF TEMPORAL CORRELATIONS

    // mean subtract first
    for (int f=0; f<nFrames; f++) CGLH[f] -= vSRRF_AVE;

    // calculate auto-correlation
    float SRRFCumTimeLags[6]; // cumulative timelags
    for (int tl = 0; tl<7; tl++) SRRFCumTimeLags[tl] = 0; // initialise SRRFCumTimeLags at 0
    int nFrames1 = nFrames-1;

    for (int f=0; f<nFrames; f++) {
        int f1 = f+1;
        float A = fabs(CGLH[f]);
        float B = fabs(CGLH[min(f+1, nFrames1)]);
        float C = fabs(CGLH[min(f+2, nFrames1)]);
        float D = fabs(CGLH[min(f+3, nFrames1)]);
        float E = fabs(CGLH[min(f+4, nFrames1)]);
        float F = fabs(CGLH[min(f+5, nFrames1)]);
        float G = fabs(CGLH[min(f+6, nFrames1)]);

        SRRFCumTimeLags[0] += (A * A - SRRFCumTimeLags[0]) / f1;

        if (f < nFrames - 1) SRRFCumTimeLags[1] += (A * B - SRRFCumTimeLags[1]) / f1;
        if (f < nFrames - 2) SRRFCumTimeLags[2] += (A * B * C - SRRFCumTimeLags[2]) / f1;
        if (f < nFrames - 3) SRRFCumTimeLags[3] += (A * B * C * D - SRRFCumTimeLags[3]) / f1;
        if (f < nFrames - 4) SRRFCumTimeLags[4] += (A * B * C * D * E - SRRFCumTimeLags[4]) / f1;
        if (f < nFrames - 5) SRRFCumTimeLags[5] += (A * B * C * D * E * F - SRRFCumTimeLags[5]) / f1;
        if (f < nFrames - 6) SRRFCumTimeLags[6] += (A * B * C * D * E * F * G - SRRFCumTimeLags[6]) / f1;
    }

    SRRFArray[0 * whM + xyMOffset] = vRAW_AVE;
    SRRFArray[1 * whM + xyMOffset] = vSRRF_MAX * vRAW_AVE;
    SRRFArray[2 * whM + xyMOffset] = vSRRF_AVE * vRAW_AVE;
    SRRFArray[3 * whM + xyMOffset] = sqrt(SRRFCumTimeLags[0]) * vRAW_AVE;
    SRRFArray[4 * whM + xyMOffset] = sqrt(SRRFCumTimeLags[1]) * vRAW_AVE;
    SRRFArray[5 * whM + xyMOffset] = pow(SRRFCumTimeLags[2], 1/3.f) * vRAW_AVE;
    SRRFArray[6 * whM + xyMOffset] = pow(SRRFCumTimeLags[3], 1/4.f) * vRAW_AVE;
    SRRFArray[7 * whM + xyMOffset] = pow(SRRFCumTimeLags[4], 1/5.f) * vRAW_AVE;
    SRRFArray[8 * whM + xyMOffset] = pow(SRRFCumTimeLags[5], 1/6.f) * vRAW_AVE;
}

