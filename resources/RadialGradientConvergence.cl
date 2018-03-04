//#pragma OPENCL EXTENSION cl_khr_fp64: enable

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

static float getInterpolatedValue(__global float* array, int const width, int const height, float const x, float const y) {
    int u0 = (int) floor(x - 0.5f);
    int v0 = (int) floor(y - 0.5f);
    float q = 0.0f;
    for (int j = 0; j <= 3; j++) {
        int v = min(max(v0 - 1 + j, 0), height-1);
        float p = 0.0f;
        for (int i = 0; i <= 3; i++) {
            int u = min(max(u0 - 1 + i, 0), width-1);
            p = p + array[v*width+u] * cubic(x - (u + 0.5f));
        }
        q = q + p * cubic(y - (v + 0.5f));
    }
    return q;
}

static float getVBoundaryCheck(__global float* array, int const width, int const height, int const x, int const y) {
    int _x = min(max(x, 0), width-1);
    int _y = min(max(y, 0), height-1);
    return array[_y*width+_x];
}


// First kernel: evaluating gradient from image --------------------------------------------------------------
__kernel void calculateGradient(
    __global float* pixels,
    __global float* GxArray,
    __global float* GyArray
    ) {
    int x = get_global_id(0);
    int y = get_global_id(1);
    int w = get_global_size(0);
    int h = get_global_size(1);
    int offset = y * w + x;

    int x0 = max(x-1, 0);
    int x1 = min(x+1, w-1);
    int y0 = max(y-1, 0);
    int y1 = min(y+1, h-1);

    GxArray[offset] = - pixels[y * w + x0] + pixels[y * w + x1];
    GyArray[offset] = - pixels[y0 * w + x] + pixels[y1 * w + x];
}


// First kernel: evaluating gradient from image using Robert's cross ------------------------------------------
__kernel void calculateGradientRobX(
    __global float* pixels,
    __global float* GxArray,
    __global float* GyArray
    ) {
    int x0 = get_global_id(0);
    int y0 = get_global_id(1);
    int w = get_global_size(0);
    int h = get_global_size(1);
    int offset = y0 * w + x0;

    int x1 = min(x0+1, w-1);
    int y1 = min(y0+1, h-1);

    // This calculates Robert's cross gradient and apply the rotation matrix 45 degrees to realign Gx and Gy to the image grid
    GxArray[offset] = pixels[y0 * w + x1] - pixels[y1 * w + x0] - pixels[y1 * w + x1] + pixels[y0 * w + x0];
    GyArray[offset] = pixels[y0 * w + x1] - pixels[y1 * w + x0] + pixels[y1 * w + x1] - pixels[y0 * w + x0];
}


// Second kernel: evaluating RGC from gradient -----------------------------------------------------------------
__kernel void calculateRadialGradientConvergence(
    __global float* pixels,
    __global float* GxArray,
    __global float* GyArray,
    __global float* RGCArray,
    int const magnification,
    float const fwhm,
    float const shiftX,
    float const shiftY,
    float const vxy_offset
    ) {

    int xM = get_global_id(0);
    int yM = get_global_id(1);
    int wM = get_global_size(0);
    int hM = get_global_size(1);
    int w = wM / magnification;
    int h = hM / magnification;
    int offset = yM * wM + xM;

    float xc = (xM + 0.5) / magnification + shiftX;
    float yc = (yM + 0.5) / magnification + shiftY;

    float CGLH = 0; // CGLH stands for Radiality original name - Culley-Gustafsson-Laine-Henriques transform
    float distanceWeightSum = 0;

    float vx, vy, Gx, Gy;
    float sigma = fwhm / 2.354f; // Sigma = 0.21 * lambda/NA in theory
    float fradius = sigma * 2;
    int radius = (int) fradius + 1;    // radius can be set to something sensible like 3*Sigma

    for (int j=-radius; j<=radius; j++) {
        for (int i=-radius; i<=radius; i++) {
            vx = (int) xc + i + 0.5 + vxy_offset;
            vy = (int) yc + j + 0.5 + vxy_offset;

            float distance = sqrt(pow(vx - xc, 2)+pow(vy - yc, 2));    // Distance D

            if (distance != 0) {
                Gx = getVBoundaryCheck(GxArray, w, h, vx-vxy_offset, vy-vxy_offset);
                Gy = getVBoundaryCheck(GyArray, w, h, vx-vxy_offset, vy-vxy_offset);

                float GMag = sqrt(Gx * Gx + Gy * Gy);

                float distanceWeight = (distance/(sigma*sigma*sigma))*exp(-(distance*distance)/(2*sigma*sigma));  // can use Taylor expansion there
                distanceWeight = distanceWeight * distanceWeight;

                // Calculate perpendicular distance from (xc,yc) to gradient line through (vx,vy)
                float Dk = fabs(Gy * (xc - vx) - Gx * (yc - vy)) / GMag;    // Dk = D*sin(theta)
                if (isnan(Dk)) Dk = distance; // this makes Dk = 0 in the next line

                Dk = 1 - Dk / distance; // Dk is now between 0 to 1, 1 if vector points precisely to (xc, yx), Dk = 1-sin(theta)
                Dk = fmax(Dk - 0.5f, 0)*2;
                //Dk = Dk*Dk*Dk*Dk;   // i think it's better to apply non-linear functions at the CGH level
//                if (Dk >= 0.75) Dk = 1;
//                else Dk = 0;
                Dk *= distanceWeight;
                distanceWeightSum += distanceWeight;

                // Accumulate Variables
                float GdotR = (Gx * i * magnification + Gy * j * magnification); // tells you if vector was pointing inward or outward
                if (GdotR <= 0) CGLH += Dk; // vector was pointing inwards
                //else CGLH -= Dk; // vector was pointing outwards
            }
        }
    }
    CGLH /= distanceWeightSum;
//    if (CGLH >= 0) CGLH = pow(CGLH, radialitySensitivity);
//    else CGLH = 0;

    float v = getInterpolatedValue(pixels, w, h, ((float) xM)/magnification + shiftX, ((float) yM)/magnification + shiftY);
    RGCArray[offset] = v * CGLH;
}


