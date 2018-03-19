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
    if (x < 0.0f) x = -x;
    float z = 0.0f;
    if (x < 1.0f)
        z = x * x * (x * (1.5f) + (- 2.5f)) + 1.0f;
    else if (x < 2.0f)
        z = -.5f * x * x * x + 2.5f * x * x - 4.0f * x + 2.0f;
    return z;
}

// interpolation function: interpolate in continuous space with respect to the reference of the array
static float getInterpolatedValue(__global float* array, int const width, int const height, float const x, float const y, float const f) { // TODO: review the grid position in the interpolation (seems offset)
    const int fwh = f * width * height;
    const int u0 = (int) floor(x - 0.5f);
    const int v0 = (int) floor(y - 0.5f);
    float q = 0.0f;
    for (int j = 0; j <= 3; j++) {
        int v = min(max(v0 - 1 + j, 0), height-1);
        float p = 0.0f;
        for (int i = 0; i <= 3; i++) {
            int u = min(max(u0 - 1 + i, 0), width-1);
            p = p + array[fwh + v * width + u] * cubic(x - (u + 0.5f));
        }
        q = q + p * cubic(y - (v + 0.5f));
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
        float xc = (xM + 0.5f) / magnification + ShiftXArray[f] + .5f; // continuous space position at the centre of magnified pixel
        float yc = (yM + 0.5f) / magnification + ShiftYArray[f] + .5f; // this last .5f sum align the RGC to the interpolated magnified image - don't know why... but works

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

    SRRFArray[0 * whM + xyMOffset] = vRAW_AVE;
    SRRFArray[1 * whM + xyMOffset] = vSRRF_MAX * vRAW_AVE;
    SRRFArray[2 * whM + xyMOffset] = vSRRF_AVE * vRAW_AVE;

    // CALCULATE SRRF TEMPORAL CORRELATIONS

    // mean subtract first
    for (int f=0; f<nFrames; f++) CGLH[f] -= vSRRF_AVE;

    // calculate auto-correlation
    float vSRRF_1ST = 0;
    float vSRRF_2ND = 0;
    float vSRRF_3RD = 0;
    float vSRRF_4TH = 0;

    for (int f=0; f<nFrames; f++) {
        vSRRF_1ST += (CGLH[f] * CGLH[f] - vSRRF_1ST) / (f+1);
        if (f<nFrames-1)
            vSRRF_2ND += (fabs(CGLH[f] * CGLH[f+1]) - vSRRF_2ND) / (f+1);
        if (f<nFrames-2)
            vSRRF_3RD += (fabs(CGLH[f] * CGLH[f+1] * CGLH[f+2]) - vSRRF_3RD) / (f+1);
        if (f<nFrames-3) {
            float A = CGLH[f];
            float B = CGLH[f+1];
            float C = CGLH[f+2];
            float D = CGLH[f+3];
            float ABCD = A * B * C * D;
            float AB = A * B;
            float CD = C * D;
            float AC = A * C;
            float BD = B * D;
            float AD = A * D;
            float BC = B * C;
            vSRRF_4TH += (fabs(ABCD - AB * CD - AC * BD - AD * BC) - vSRRF_4TH) / (f+1);
        }

    }
    vSRRF_1ST = sqrt(vSRRF_1ST);
    vSRRF_2ND = sqrt(vSRRF_2ND);
    vSRRF_3RD = pow(vSRRF_3RD, 0.3333333f);
    vSRRF_4TH = pow(vSRRF_4TH, 0.25f);
    SRRFArray[3 * whM + xyMOffset] = vSRRF_1ST * vRAW_AVE;
    SRRFArray[4 * whM + xyMOffset] = vSRRF_2ND * vRAW_AVE;
    SRRFArray[5 * whM + xyMOffset] = vSRRF_3RD * vRAW_AVE;
    SRRFArray[6 * whM + xyMOffset] = vSRRF_4TH * vRAW_AVE;
}

