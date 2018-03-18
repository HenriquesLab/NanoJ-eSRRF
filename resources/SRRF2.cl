//#pragma OPENCL EXTENSION cl_khr_fp64: enable

// cubic function for interpolation
static float cubic(float x) {
    float a = 0.5f; // Catmull-Rom interpolation
    if (x < 0.0f) x = -x;
    float z = 0.0f;
    if (x < 1.0f)
        z = x * x * (x * (-a + 2.0f) + (a - 3.0f)) + 1.0f;
    else if (x < 2.0f)
        z = -a * x * x * x + 5.0f * a * x * x - 8.0f * a * x + 4.0f * a;
    return z;
}

// interpolation function: interpolate in continuous space with respect to the reference of the array
static float getInterpolatedValue(__global float* array, int const width, int const height, float const x, float const y, float const f) { // TODO: review the grid position in the interpolation (seems offset)
    int fwh = f * width * height;
    int u0 = (int) floor(x - 0.5f);
    int v0 = (int) floor(y - 0.5f);
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

// check boundaries of the image and returns the gradient value
static float getVBoundaryCheck(__global float* array, int const width, int const height, int const x, int const y) {
    int _x = min(max(x, 0), width-1);
    int _y = min(max(y, 0), height-1);
    return array[_y*width+_x];
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
    int w = get_global_size(0);
    int h = get_global_size(1);
    int wh = w * h;
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
    int const nFrames,
    int const magnification,
    float const fwhm
    ) {

    int xM = get_global_id(0);
    int yM = get_global_id(1);
    int wM = get_global_size(0);
    int hM = get_global_size(1);
    int w = wM / magnification;
    int h = hM / magnification;
    int wh = w * h;
    int whM = wM * hM;
    int xyMOffset = yM * wM + xM;

    float CGLH[MAX_FRAMES]; // note, MAX_FRAMES is passed before compile
    float vSRRF_AVE = 0;
    float vSRRF_MAX = 0;

    // FIRST CALCULATE LOCAL GRADIENT CONVERGENCE (OLD RADIALITY!!)

    float sigma = fwhm / 2.354f; // Sigma = 0.21 * lambda/NA in theory
    float sigma22 = 2 * sigma * sigma;
    int radius = (int) (sigma * 2) + 1;    // radius can be set to something sensible like 3*Sigma
    float Gx, Gy;

    for (int f=0; f<nFrames; f++) {
        int fOffset = f * wh;
        float xc = (xM + 0.5f) / magnification + ShiftXArray[f]; // continuous space position at the centre of magnified pixel
        float yc = (yM + 0.5f) / magnification + ShiftYArray[f];

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
                float distance = sqrt(dx * dx+ dy * dy);

                if (distance != 0 && vxPixelCentred >=0 && vxPixelCentred < w && vyPixelCentred >=0 && vyPixelCentred < h) {
                    int p = fOffset + vyPixelOrigin * w + vxPixelOrigin;
                    Gx = GxArray[p];
                    Gy = GyArray[p];

                    float GMag = sqrt(Gx * Gx + Gy * Gy);

                    float distanceWeight = distance*exp(-(distance*distance)/sigma22);  // TODO: dGauss: can use Taylor expansion there

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
                    //else CGLH -= Dk; // vector was pointing outwards
                }
            }
        }

        CGLH[f] /= distanceWeightSum;
        CGLH[f] *= getInterpolatedValue(pixels, w, h, ((float) xM)/magnification + ShiftXArray[f], ((float) yM)/magnification + ShiftYArray[f], f);

        // CALCULATE SRRF AVERAGE AND MAXIMUM
        vSRRF_MAX = fmax(CGLH[f], vSRRF_MAX);
        vSRRF_AVE += (CGLH[f] - vSRRF_AVE) / (f+1);
    }

    SRRFArray[0 * whM + xyMOffset] = vSRRF_MAX;
    SRRFArray[1 * whM + xyMOffset] = vSRRF_AVE;

    // CALCULATE SRRF TEMPORAL CORRELATIONS

    // mean subtract first
    for (int f=0; f<nFrames; f++) CGLH[f] -= vSRRF_AVE;

    // now calculate temporal correlations
    for (int tl = 0; tl < 8; tl++) { // assuming a max of 8 time-lags
        int counter = 0;
        double covariance = 0;

        for (int f=tl; f<nFrames; f++) {
            counter++;
            float v0 = CGLH[f-tl];
            float v1 = CGLH[f];
            covariance += (fabs(v0 * v1) - covariance) / counter;
        }
        //covariance = sqrt(covariance);
        SRRFArray[(2+tl) * whM + xyMOffset] = covariance;
    }
}

