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
static float getInterpolatedValue(__global float* array, int const width, int const height, float const x, float const y) { // TODO: review the grid position in the interpolation (seems offset)
    int u0 = (int) floor(x);
    int v0 = (int) floor(y);
    float q = 0.0f;
    for (int j = 0; j <= 3; j++) {
        int v = min(max(v0 - 1 + j, 0), height-1);
        float p = 0.0f;
        for (int i = 0; i <= 3; i++) {
            int u = min(max(u0 - 1 + i, 0), width-1);
            p = p + array[v*width+u] * cubic(x - (float) (u));
        }
        q = q + p * cubic(y - (float) (v));
    }
    return q;
}

// check boundaries of the image and returns the gradient value
static float getVBoundaryCheck(__global float* array, int const width, int const height, int const x, int const y) {
    int _x = min(max(x, 0), width-1);
    int _y = min(max(y, 0), height-1);
    return array[_y*width+_x];
}


// First kernel: evaluating gradient from image: 3-point gradient --------------------------------------------------------------
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
    int x1 = get_global_id(0);
    int y1 = get_global_id(1);
    int w = get_global_size(0);
    int h = get_global_size(1);
    int offset = y1 * w + x1;

    int x0 = max(x1-1, 0);
    int y0 = max(y1-1, 0);

    // This calculates Robert's cross gradient and apply the rotation matrix 45 degrees to realign Gx and Gy to the image grid
    GxArray[offset] = pixels[y0 * w + x1] - pixels[y1 * w + x0] + pixels[y1 * w + x1] - pixels[y0 * w + x0];
    GyArray[offset] = - pixels[y0 * w + x1] + pixels[y1 * w + x0] + pixels[y1 * w + x1] - pixels[y0 * w + x0];
}

// First kernel: evaluating gradient from image using 2-point gradient --------------------------------------------------
__kernel void calculateGradient_2point(
    __global float* pixels,
    __global float* GxArray,
    __global float* GyArray
    ) {
    int x1 = get_global_id(0);
    int y1 = get_global_id(1);
    int w = get_global_size(0);
    int h = get_global_size(1);
    int offset = y1 * w + x1;

    int x0 = max(x1-1, 0);
    int y0 = max(y1-1, 0);

    // 2-point gradient
    GxArray[offset] = pixels[y1 * w + x1] - pixels[y1 * w + x0];
    GyArray[offset] = pixels[y1 * w + x1] - pixels[y0 * w + x1];
}

// First kernel follow up: interpolating gradient from 2-point gradient image -------------------------------------------
__kernel void calculateGradient2p_Interpolation(
    __global float* GxArray,
    __global float* GyArray,
    __global float* GxIntArray,
    __global float* GyIntArray
    ) {
    int x = get_global_id(0);
    int y = get_global_id(1);
    int wInt = get_global_size(0);
    int hInt = get_global_size(1);
    int offset = y * wInt + x;

    // Two-fold interpolation of the gradients
    GxIntArray[offset] = getInterpolatedValue(GxArray, (int) (wInt/2), (int) (hInt/2), (float) (x)/2.0f, (float) (y)/2.0f);
    GyIntArray[offset] = getInterpolatedValue(GyArray, (int) (wInt/2), (int) (hInt/2), (float) (x)/2.0f, (float) (y)/2.0f);
}



// Second kernel: evaluating RGC from gradient -----------------------------------------------------------------
__kernel void calculateRadialGradientConvergence(
    __global float* pixels,
    __global float* GxArray,
    __global float* GyArray,
    __global float* RGCArray,

//    __global float* wSumArray,
//    __global float* debugFunArray,

    int const magnification,
    int const GxGyMagnification,
    float const fwhm,
    int const sensitivity,
    float const shiftX,
    float const shiftY,
    float const vxy_offset,
    int const vxy_ArrayShift,
    float const vxy_PixelShift,
    int const intWeighting

    ) {

    int xM = get_global_id(0);
    int yM = get_global_id(1);
    int wM = get_global_size(0);
    int hM = get_global_size(1);
    int w = wM / magnification;
    int h = hM / magnification;
    int offset = yM * wM + xM;

    float xc = (xM + 0.5) / magnification + shiftX; // continuous space position at the centre of magnified pixel
    float yc = (yM + 0.5) / magnification + shiftY;

    float CGLH = 0; // CGLH stands for Radiality original name - Culley-Gustafsson-Laine-Henriques transform
    float distanceWeightSum = 0;

//    int countdebugFun = 0;

    float vx, vy, Gx, Gy;
    float sigma = fwhm / 2.354f; // Sigma = 0.21 * lambda/NA in theory
    float fradius = sigma * 2;
    float radius = ((float) ((int) (GxGyMagnification*fradius)))/GxGyMagnification + 1;    // this reduces the radius for speed, works when using dGauss^4 and 2p+I
//    int radius = (int) (fradius) + 1;    // this should be used otherwise


    for (int j=-GxGyMagnification*radius; j<=(GxGyMagnification*radius+1); j++) {
//        vy = ((int) (GxGyMagnification*(yc - vxy_PixelShift)) + j)/GxGyMagnification + vxy_PixelShift; // position in continuous space TODO: problems with negative values and (int)?
        vy = ((float) ((int) (GxGyMagnification*(yc - vxy_PixelShift))) + j)/GxGyMagnification + vxy_PixelShift; // position in continuous space TODO: problems with negative values and (int)?


        for (int i=-GxGyMagnification*radius; i<=(GxGyMagnification*radius+1); i++) {
//            vx = ((int) (GxGyMagnification*(xc - vxy_PixelShift)) + i)/GxGyMagnification + vxy_PixelShift; // position in continuous space TODO: problems with negative values and (int)?
            vx = ((float) ((int) (GxGyMagnification*(xc - vxy_PixelShift))) + i)/GxGyMagnification + vxy_PixelShift; // position in continuous space TODO: problems with negative values and (int)?

            float distance = sqrt((vx - xc)*(vx - xc) + (vy - yc)*(vy - yc));    // Distance D

//            if (j<countdebugFun) countdebugFun = j;

//            if (distance != 0 && distance <= (fwhm/2)) {
//            if (distance != 0 && distance <= (2*sigma+1) && GdotR <= 0) {
            if (distance != 0 && distance <= (2*sigma+1)) {
//            countdebugFun +=1;


                Gx = getVBoundaryCheck(GxArray, GxGyMagnification*w, GxGyMagnification*h, GxGyMagnification*(vx - vxy_offset) + vxy_ArrayShift, GxGyMagnification*(vy - vxy_offset));
                Gy = getVBoundaryCheck(GyArray, GxGyMagnification*w, GxGyMagnification*h, GxGyMagnification*(vx - vxy_offset), GxGyMagnification*(vy - vxy_offset) + vxy_ArrayShift);

                float distanceWeight = distance*exp(-(distance*distance)/(2*sigma*sigma));  // TODO: dGauss: can use Taylor expansion there
                distanceWeight = distanceWeight * distanceWeight * distanceWeight * distanceWeight ;  // TODO: dGauss: what power is best? Let's FRC !
                distanceWeightSum += distanceWeight;
                float GdotR = (Gx * (vx - xc) + Gy * (vy - yc)); // tells you if vector was pointing inward or outward


                if (GdotR < 0) {
                    // Calculate perpendicular distance from (xc,yc) to gradient line through (vx,vy)
                    float GMag = sqrt(Gx * Gx + Gy * Gy);
                    float Dk = fabs(Gy * (xc - vx) - Gx * (yc - vy)) / GMag;    // Dk = D*sin(theta) obtained from cross-product
                    if (isnan(Dk)) Dk = distance; // this makes Dk = 0 in the next line


                    // Linear function ----------------------
                    Dk = 1 - Dk / distance; // Dk is now between 0 to 1, 1 if vector points precisely to (xc, yx), Dk = 1-sin(theta)

                    //Linear truncated ----------------------
//                  Dk = 1 - Dk / distance; // Dk is now between 0 to 1, 1 if vector points precisely to (xc, yx), Dk = 1-sin(theta)
//                  Dk = fmax(Dk - 0.5f, 0)*2;

                    // Higher order of Dk (quadratic and power 4) -------------------
//                  Dk = 1 - Dk / distance; // Dk is now between 0 to 1, 1 if vector points precisely to (xc, yx), Dk = 1-sin(theta)
//                  Dk = Dk*Dk;   // i think it's better to apply non-linear functions at the CGH level

                    // Hard edge function -----------------------
//                  Dk = 1 - Dk / distance; // Dk is now between 0 to 1, 1 if vector points precisely to (xc, yx), Dk = 1-sin(theta)
//                  float edgePos = 0.75;
//                  if (Dk >= edgePos) Dk = 1;
//                  else Dk = 0;

                    // Gaussian function --------------------------
//                  Dk = Dk / distance; // Dk is now between 0 to 1, Dk = sin(theta) ~ theta
//                  float SigmaTheta = 0.3;
//                  Dk = exp(-0.5f*pow((float) (Dk/SigmaTheta), 2));

                    // Rational function -------------------------
//                  float SigmaTheta = 0.3;
//                  float Power = 4;
//                  Dk = Dk / distance; // Dk is now between 0 to 1, Dk = sin(theta) ~ theta
//                  Dk = 1/(1+pow((float) (Dk/SigmaTheta), Power));

                    // Exponential function --------------------------
//                  Dk = Dk / distance; // Dk is now between 0 to 1, Dk = sin(theta) ~ theta
//                  float SigmaTheta = 0.25;
//                  Dk = exp(-1.0f*(float) (Dk/SigmaTheta));

                    CGLH += Dk*distanceWeight;
                }


//                Dk *= distanceWeight;
//
//                // Accumulate Variables
//                float GdotR = (Gx * i + Gy * j); // tells you if vector was pointing inward or outward
//                if (GdotR <= 0) CGLH += Dk; // vector was pointing inwards
//                //else CGLH -= Dk; // vector was pointing outwards
            }
        }
    }

//    debugFunArray[offset] = countdebugFun;
//    wSumArray[offset] = distanceWeightSum;

    CGLH /= distanceWeightSum;
    if (CGLH >= 0) CGLH = pow(CGLH, sensitivity);
    else CGLH = 0;

    if (intWeighting == 1) {
        float v = getInterpolatedValue(pixels, w, h, ((float) xM)/magnification + shiftX - 0.5f, ((float) yM)/magnification + shiftY - 0.5f);
        RGCArray[offset] = v * CGLH;}
    else RGCArray[offset] = CGLH;
}


// Third kernel: calculating the interpolated intensity image -----------------------------------------------------------------
__kernel void kernelCalculateInterpolatedIntensity(
    __global float* pixels,
    __global float* intArray,
    int const magnification,
    float const shiftX,
    float const shiftY
    ){

    int xM = get_global_id(0);
    int yM = get_global_id(1);
    int wM = get_global_size(0);
    int hM = get_global_size(1);
    int w = wM / magnification;
    int h = hM / magnification;
    int offset = yM * wM + xM;

    float v = getInterpolatedValue(pixels, w, h, ((float) xM)/magnification + shiftX - 0.5f, ((float) yM)/magnification + shiftY - 0.5f);
    intArray[offset] = v;

}