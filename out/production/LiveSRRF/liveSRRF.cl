//#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define magnification $MAGNIFICATION$
//#define fwhm $FWHM$
#define sensitivity $SENSITIVITY$
#define GxGyMagnification $GXGYMAGNIFICATION$
//#define sigma $SIGMA$
#define tSS $TWOSIGSQUARE$
#define tSO $TWOSIGpONE$
#define radius $RADIUS$
#define width $WIDTH$
#define height $HEIGHT$
#define wh $WH$
#define wInt $WINT$
#define hInt $HINT$
#define whInt $WHINT$
#define wM $WM$
#define hM $HM$
#define whM $WHM$
#define nFrameForSRRF $NFRAMEFORSRRF$
#define intWeighting $INTWEIGHTING$

#define vxy_offset $VXY_OFFSET$
#define vxy_ArrayShift $VXY_ARRAYSHIFT$

// Cubic function for interpolation
static float cubic(float x) {
    const float a = 0.5f; // Catmull-Rom interpolation
    if (x < 0.0f) x = -x;
    float z = 0.0f;
    if (x < 1.0f)
        z = x * x * (x * (-a + 2.0f) + (a - 3.0f)) + 1.0f;
    else if (x < 2.0f)
        z = -a * x * x * x + 5.0f * a * x * x - 8.0f * a * x + 4.0f * a;
    return z;
}

// Interpolation function: interpolate in continuous space with respect to the reference of the array // TODO: check confusion between width and w height and h
//static float getInterpolatedValue(__global float* array, int const width, int const height, float const x, float const y, int const f) {
static float getInterpolatedValue(__global float* array, float const x, float const y, int const f) {
    const int u0 = (int) floor(x);
    const int v0 = (int) floor(y);
//    const int whf = width*height*f;
    const int whf = wh*f;

    float q = 0.0f;

    // Bicubic interpolation
    if (u0 > 0 && u0 < width - 2 && v0 > 0 && v0 < height - 2) {
        for (int j = 0; j <= 3; j++) {
            int v = min(max(v0 - 1 + j, 0), height-1);
            float p = 0.0f;
            for (int i = 0; i <= 3; i++) {
                int u = min(max(u0 - 1 + i, 0), width-1);
                p = p + array[v*width + u + whf] * cubic(x - (float) (u));
            }
            q = q + p * cubic(y - (float) (v));
        }
    }

//    // Bilinear interpolation
//    else if (x >= 0 && x < w-1 && y >= 0 && y < h-1) {
//
//        int xbase = (int)x;
//        int ybase = (int)y;
//        int xbase1 = min(xbase+1, w-1);
//        int ybase1 = min(ybase+1, h-1);
//
//        float xFraction = x - (float) xbase;
//        float yFraction = y - (float) ybase;
//        xFraction = fmax(xFraction, 0);
//        yFraction = fmax(yFraction, 0);
//
//        float lowerLeft = array[whf + ybase * width + xbase];
//        float lowerRight = array[whf + ybase * width + xbase1];
//        float upperRight = array[whf + ybase1 * width + xbase1];
//        float upperLeft = array[whf + ybase1 * width + xbase];
//        float upperAverage = upperLeft + xFraction * (upperRight - upperLeft);
//        float lowerAverage = lowerLeft + xFraction * (lowerRight - lowerLeft);
//        q = lowerAverage + yFraction * (upperAverage - lowerAverage);
//    }

    // Extrapolation
    else {

//        if (x < 0){
//            xbase = 0;
//            xbase1 = 1
//            }
//        else if (x >= w-1){
//            xbase = w-2;
//            xbase1 = w-1;
//            }
//        else {
//            xbase = (int)x;
//            xbase1 = xbase+1;
//            }
//
//        if (y < 0){
//            ybase = 0;
//            ybase1 = 1
//            }
//        else if (y >= h-1){
//            ybase = h-2;
//            ybase1 = h-1;
//            }
//        else {
//            ybase = (int)y;
//            ybase1 = ybase+1;
//            }

        int xbase = (int) fmin((float) width-2, fmax(x,0.0f));
        int xbase1 = xbase+1;

        int ybase = (int) fmin((float) height-2, fmax(y,0.0f));
        int ybase1 = ybase+1;

        float xFraction = x - (float) xbase;
        float yFraction = y - (float) ybase;
//        xFraction = fmax(xFraction, 0);
//        yFraction = fmax(yFraction, 0);

        float lowerLeft = array[whf + ybase * width + xbase];
        float lowerRight = array[whf + ybase * width + xbase1];
        float upperRight = array[whf + ybase1 * width + xbase1];
        float upperLeft = array[whf + ybase1 * width + xbase];
        float upperAverage = upperLeft + xFraction * (upperRight - upperLeft);
        float lowerAverage = lowerLeft + xFraction * (lowerRight - lowerLeft);
        q = lowerAverage + yFraction * (upperAverage - lowerAverage);

    }

    return q;
}

// Check boundaries of the image and returns the gradient value // TODO: extrapolate instead of boundary check?
static float getVBoundaryCheck(__global float* array, int const thisWidth, int const thisHeight, int const x, int const y, int const f) {
    const int _x = min(max(x, 0), thisWidth-1);
    const int _y = min(max(y, 0), thisHeight-1);
    return array[_x + _y*thisWidth + thisWidth*thisHeight*f];
}


//// First kernel: evaluating gradient from image: 3-point gradient --------------------------------------------------------------
//__kernel void calculateGradient(
//    __global float* pixels,
//    __global float* GxArray,
//    __global float* GyArray
//    ) {
//    int x = get_global_id(0);
//    int y = get_global_id(1);
//    int w = get_global_size(0);
//    int h = get_global_size(1);
//    int offset = y * w + x;
//
//    int x0 = max(x-1, 0);
//    int x1 = min(x+1, w-1);
//    int y0 = max(y-1, 0);
//    int y1 = min(y+1, h-1);
//
//    GxArray[offset] = - pixels[y * w + x0] + pixels[y * w + x1];
//    GyArray[offset] = - pixels[y0 * w + x] + pixels[y1 * w + x];
//}


//// First kernel: evaluating gradient from image using Robert's cross ------------------------------------------
//__kernel void calculateGradientRobX(
//    __global float* pixels,
//    __global float* GxArray,
//    __global float* GyArray
//    ) {
//    int x1 = get_global_id(0);
//    int y1 = get_global_id(1);
//    int w = get_global_size(0);
//    int h = get_global_size(1);
//    int offset = y1 * w + x1;
//
//    int x0 = max(x1-1, 0);
//    int y0 = max(y1-1, 0);
//
//    // This calculates Robert's cross gradient and apply the rotation matrix 45 degrees to realign Gx and Gy to the image grid
//    GxArray[offset] = pixels[y0 * w + x1] - pixels[y1 * w + x0] + pixels[y1 * w + x1] - pixels[y0 * w + x0];
//    GyArray[offset] = - pixels[y0 * w + x1] + pixels[y1 * w + x0] + pixels[y1 * w + x1] - pixels[y0 * w + x0];
//}


// First kernel: evaluating gradient from image using 2-point gradient --------------------------------------------------
__kernel void calculateGradient_2point(
    __global float* pixels,
    __global float* GxArray,
    __global float* GyArray,
    __global int* nCurrentFrame

    ) {
    const int x1 = get_global_id(0);
    const int y1 = get_global_id(1);
    const int f = get_global_id(2);

    const int offset = y1 * width + x1 + wh * f;
    const int x0 = max(x1-1, 0);
    const int y0 = max(y1-1, 0);

    // 2-point gradient
    GxArray[offset] = pixels[offset] - pixels[y1 * width + x0 + wh * f];
    GyArray[offset] = pixels[offset] - pixels[y0 * width + x1 + wh * f];

    // Reset the local current frame
    nCurrentFrame[1] = 0;
//    if (nCurrentFrame[0] == nFrameForSRRF) nCurrentFrame[0] = 0; // reset the frame number if it's reached the end

        // Current frame is a 2 element Int buffer:
                // nCurrentFrame[0] is the global current frame in the current SRRF frame (reset every SRRF frame)
                // nCurrentFrame[1] is the local current frame in the current GPU-loaded dataset (reset every turn of the method calculateSRRF (within the gradient calculation))

}

// First kernel follow up: interpolating gradient from 2-point gradient image -------------------------------------------
__kernel void calculateGradientInterpolation(
    __global float* GxArray,
    __global float* GyArray,
    __global float* GxIntArray,
    __global float* GyIntArray

    ){
    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int f = get_global_id(2);

    const int offset = y * wInt + x + f * whInt;

    // Two-fold interpolation of the gradients
    GxIntArray[offset] = getInterpolatedValue(GxArray, (float) (x)/2.0f, (float) (y)/2.0f, f);
    GyIntArray[offset] = getInterpolatedValue(GyArray, (float) (x)/2.0f, (float) (y)/2.0f, f);
}



// Second kernel: evaluating RGC from gradient -----------------------------------------------------------------
__kernel void calculateRadialGradientConvergence(
    __global float* pixels,
    __global float* GxArray,
    __global float* GyArray,
    __global float* PreviousFrameArray,
    __global float* OutArray,
    __global float* shiftXY,
    __global int* nCurrentFrame
    // Current frame is a 2 element Int buffer:
            // nCurrentFrame[0] is the global current frame in the current SRRF frame (reset every SRRF frame)
            // nCurrentFrame[1] is the local current frame in the current GPU-loaded dataset (reset every turn of the method calculateSRRF (within the gradient calculation))

    ) {

//    const int xM = get_global_id(0);
//    const int yM = get_global_id(1);
//
//    const int offset = yM * wM + xM;

    const int offset = get_global_id(0);
    const int yM = offset/wM;
    const int xM = offset - yM*wM;

    const float shiftX = shiftXY[nCurrentFrame[0]];
    const float shiftY = shiftXY[nCurrentFrame[0] + nFrameForSRRF];
    const float xc = (xM + 0.5) / magnification + shiftX; // continuous space position at the centre of magnified pixel
    const float yc = (yM + 0.5) / magnification + shiftY;
//    const float sigma22 = 2 * sigma * sigma; // TODO: add as hardcoded value? Something wrong happens when doing that
//    const float sigma21 = 2 * sigma + 1;

    float CGLH = 0; // CGLH stands for Radiality original name - Culley-Gustafsson-Laine-Henriques transform
    float distanceWeightSum = 0;

    float vx, vy, Gx, Gy;
//    float sigma = fwhm / 2.354f; // Sigma = 0.21 * lambda/NA in theory
//    float fradius = sigma * 2;
//    float radius = ((float) ((int) (GxGyMagnification*fradius)))/GxGyMagnification + 1;    // this reduces the radius for speed, works when using dGauss^4 and 2p+I
//    int radius = (int) (fradius) + 1;    // this should be used otherwise
    float dx, dy; // TODO: define dx and dy as the differences of coordinates
    float distance, distanceWeight, GdotR, GMag, Dk;

    for (int j=-(int) ((float) GxGyMagnification*(float) radius); j<=(int)((float) GxGyMagnification*(float) radius+1); j++) {
        vy = ((float) ((int) (GxGyMagnification*yc)) + j)/(float) GxGyMagnification; // position in continuous space CHECKED: optimised for 2-point gradient

        if (vy > 0 && vy < height){
            for (int i=-(int) ((float) GxGyMagnification*(float) radius); i<=(int)((float) GxGyMagnification*(float) radius+1); i++) {
                vx = ((float) ((int) (GxGyMagnification*xc)) + i)/(float) GxGyMagnification; // position in continuous space
                if (vx > 0 && vx < width){

                    dx = vx - xc;
                    dy = vy - yc;
                    distance = sqrt(dx*dx + dy*dy);    // Distance D

                    if (distance != 0 && distance <= (float) tSO) {

                        Gx = getVBoundaryCheck(GxArray, wInt, hInt, GxGyMagnification*(vx - vxy_offset) + vxy_ArrayShift, GxGyMagnification*(vy - vxy_offset), nCurrentFrame[1]);
                        Gy = getVBoundaryCheck(GyArray, wInt, hInt, GxGyMagnification*(vx - vxy_offset), GxGyMagnification*(vy - vxy_offset) + vxy_ArrayShift, nCurrentFrame[1]);

                        distanceWeight = distance*exp(-(distance*distance)/(float) tSS);  // TODO: dGauss: can use Taylor expansion there
                        distanceWeight = distanceWeight * distanceWeight * distanceWeight * distanceWeight ;
                        distanceWeightSum += distanceWeight;
                        GdotR = (Gx*dx + Gy*dy); // tells you if vector was pointing inward or outward

                        if (GdotR < 0) {
                            // Calculate perpendicular distance from (xc,yc) to gradient line through (vx,vy)
                            GMag = sqrt(Gx*Gx + Gy*Gy);
                            Dk = fabs(Gy*dx - Gx*dy) / GMag;    // Dk = D*sin(theta) obtained from cross-product
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
       }
    }

    CGLH /= distanceWeightSum;
    if (CGLH >= 0) CGLH = pow(CGLH, sensitivity);
    else CGLH = 0;

    float v = getInterpolatedValue(pixels, ((float) xM)/magnification + shiftX - 0.5f, ((float) yM)/magnification + shiftY - 0.5f, nCurrentFrame[1]);

    if (intWeighting == 1) {

            if (nCurrentFrame[0] == 0) { // Initialising the values of the array
                OutArray[offset] = v * CGLH / nFrameForSRRF;                    // AVG
                OutArray[offset + whM] = v * CGLH * v * CGLH / nFrameForSRRF;   // VAR
                OutArray[offset + 2 * whM] = 0;                                 // 2nd Order SOFI Tau = 1
                OutArray[offset + 3 * whM] = v / nFrameForSRRF;                 // interpolated
            }
            else {
                OutArray[offset] = OutArray[offset] + v * CGLH / nFrameForSRRF;
                OutArray[offset + whM] = OutArray[offset + whM] + v * CGLH * v * CGLH / nFrameForSRRF;
                OutArray[offset + 2 * whM] = OutArray[offset + 2 * whM] + v * CGLH * PreviousFrameArray[offset] / (nFrameForSRRF-1);
                OutArray[offset + 3 * whM] = OutArray[offset + 3 * whM] + v / nFrameForSRRF;
            }
        PreviousFrameArray[offset] = v * CGLH / (nFrameForSRRF-1); // update the value of the
    }
    else{
            if (nCurrentFrame[0] == 0) {
                OutArray[offset] = CGLH / nFrameForSRRF;
                OutArray[offset + whM] = CGLH * CGLH / nFrameForSRRF;
                OutArray[offset + 2 * whM] = 0;
                OutArray[offset + 3 * whM] = v / nFrameForSRRF;
            }
            else {
                OutArray[offset] = OutArray[offset] + CGLH / nFrameForSRRF;
                OutArray[offset + whM] = OutArray[offset + whM] + CGLH * CGLH / nFrameForSRRF;
                OutArray[offset + 2 * whM] = OutArray[offset + 2 * whM] + CGLH * PreviousFrameArray[offset] / (nFrameForSRRF-1);
                OutArray[offset + 3 * whM] = OutArray[offset + 3 * whM] + v / nFrameForSRRF;
            }
        PreviousFrameArray[offset] = CGLH / (nFrameForSRRF-1);
    }

//        }


//    else RGCArray[offset] = CGLH;


}


//// Third kernel: calculating the interpolated intensity image -----------------------------------------------------------------
//__kernel void kernelCalculateInterpolatedIntensity(
//    __global float* pixels,
//    __global float* intArray,
////    int const magnification,
//    float const shiftX,
//    float const shiftY
//    ){
//
//    int xM = get_global_id(0);
//    int yM = get_global_id(1);
////    int wM = get_global_size(0);
////    int hM = get_global_size(1);
////    int w = wM / magnification;
////    int h = hM / magnification;
//    int offset = yM * wM + xM;
//
//    float v = getInterpolatedValue(pixels, w, h, ((float) xM)/magnification + shiftX - 0.5f, ((float) yM)/magnification + shiftY - 0.5f);
//    intArray[offset] = v;
//}


// Kernel: calculate VAR image from the OutputArray
__kernel void kernelCalculateVar(
    __global float* OutArray
    ){
    const int offset = get_global_id(0);
    const float av2 = OutArray[offset] * OutArray[offset];
    OutArray[offset + whM] = OutArray[offset + whM] - av2;     // Var[X] = E[X^2] - (E[X])^2
    OutArray[offset + 2*whM] = OutArray[offset + 2*whM] - av2; // TAC2[X] = E[X(t)*X(t+1)] - (E[X])^2, TODO: this is only correct to a first approximation, could explain the negative values?

}


// Fourth kernel: increment the current frame number -----------------------------------------------------------------
__kernel void kernelIncrementFramePosition(
    __global int* nCurrentFrame
    ){
    const int i = get_global_id(0);
    nCurrentFrame[i] = nCurrentFrame[i] + 1;

        // Current frame is a 2 element Int buffer:
                // nCurrentFrame[0] is the global current frame in the current SRRF frame (reset every SRRF frame)
                // nCurrentFrame[1] is the local current frame in the current GPU-loaded dataset (reset every turn of the method calculateSRRF (within the gradient calculation))

}

// Fifth kernel: reset the frame number -----------------------------------------------------------------
__kernel void kernelResetFramePosition(
    __global int* nCurrentFrame
    ){
    nCurrentFrame[0] = 0;

        // Current frame is a 2 element Int buffer:
                // nCurrentFrame[0] is the global current frame in the current SRRF frame (reset every SRRF frame)
                // nCurrentFrame[1] is the local current frame in the current GPU-loaded dataset (reset every turn of the method calculateSRRF (within the gradient calculation))

}

// kernel: calculate the Macro-Pixel artefact map -----------------------------------------------------------------
__kernel void kernelCalculateMPmap(
    __global float* OutArray,
    __global float* MPmap
    ){

    const int offset_MPmap = get_global_id(0);
    const int frame = offset_MPmap/(magnification*magnification); // 0, 1 or 2 depending on whether we're dealing with AVG or VAR or SOFI2ndTau1
    const int y_MPmap = (offset_MPmap - frame*(magnification*magnification))/magnification;
    const int x_MPmap = offset_MPmap - y_MPmap*magnification - frame*(magnification*magnification);

    float thisMPmapValue = 0;
    int x;
    int y;
    int offset;

    for (int i=0; i<wh; i++) {
        y = i/width;
        x = i - y*width;
        offset = x_MPmap + x * magnification + wM * (y_MPmap + y * magnification) + frame * whM;
        thisMPmapValue += OutArray[offset];
    }
    MPmap[offset_MPmap] = thisMPmapValue/wh;

}

// kernel: correct for Macro-pixel artefacts -----------------------------------------------------------------
__kernel void kernelCorrectMPmap(
     __global float* OutArray,
     __global float* MPmap
){

    const int offset = get_global_id(0);
    const int frame = offset/whM; // 0, 1 or 2 depending on whether we're dealing with AVG or VAR or SOFI2ndTau1
    const int yM = (offset - frame*(whM))/wM;
    const int xM = offset - yM*wM - frame*whM;
    const int x = xM/magnification;
    const int y = yM/magnification;

    int offset_MPmap = magnification*magnification*frame + magnification*(yM - y*magnification) + xM - x*magnification;
    OutArray[offset] = OutArray[offset]/MPmap[offset_MPmap];

}
