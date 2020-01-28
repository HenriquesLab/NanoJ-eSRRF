//#pragma OPENCL EXTENSION cl_khr_fp64: enable
//#define magnification $MAGNIFICATION$
//#define fwhm $FWHM$
#define sensitivity $SENSITIVITY$
//#define GxGyMagnification $GXGYMAGNIFICATION$
//#define sigma $SIGMA$
#define tSS $TWOSIGSQUARE$
#define tSO $TWOSIGpONE$
#define radius $RADIUS$
#define width $WIDTH$
#define height $HEIGHT$
#define wh $WH$
//#define wInt $WINT$
//#define hInt $HINT$
//#define wM $WM$
//#define hM $HM$
//#define whM $WHM$
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
static float getInterpolatedValue(__global float* array, int const thisWidth, int const thisHeight, float const x, float const y, int const f) {
    const int u0 = (int) floor(x);
    const int v0 = (int) floor(y);
    const int thisWhf = thisWidth*thisHeight*f;

    float q = 0.0f;


    // Bicubic interpolation
    if (u0 > 0 && u0 < thisWidth - 2 && v0 > 0 && v0 < thisHeight - 2) {
        for (int j = 0; j <= 3; j++) {
            int v = min(max(v0 - 1 + j, 0), thisHeight-1);
            float p = 0.0f;
            for (int i = 0; i <= 3; i++) {
                int u = min(max(u0 - 1 + i, 0), thisWidth-1);
                p = p + array[v*thisWidth + u + thisWhf] * cubic(x - (float) (u));
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

        int xbase = (int) fmin((float) thisWidth-2, fmax(x,0.0f));
        int xbase1 = xbase+1;

        int ybase = (int) fmin((float) thisHeight-2, fmax(y,0.0f));
        int ybase1 = ybase+1;

        float xFraction = x - (float) xbase;
        float yFraction = y - (float) ybase;
//        xFraction = fmax(xFraction, 0);
//        yFraction = fmax(yFraction, 0);

        float lowerLeft = array[thisWhf + ybase * thisWidth + xbase];
        float lowerRight = array[thisWhf + ybase * thisWidth + xbase1];
        float upperRight = array[thisWhf + ybase1 * thisWidth + xbase1];
        float upperLeft = array[thisWhf + ybase1 * thisWidth + xbase];
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


// First kernel: evaluating gradient from image: 3-point gradient + --------------------------------------------------------------
__kernel void calculateGradient3pPlus(
    __global float* pixels,
    __global float* GxArray,
    __global float* GyArray,
    __global int* nCurrentFrame
    ) {
    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int f = get_global_id(2);

    const int frameOffset = wh*f;
    const int x0 = max(x-1, 0);
    const int x1 = min(x+1, width-1);
    const int y0 = max(y-1, 0);
    const int y1 = min(y+1, height-1);
    const int offset = y * width + x + frameOffset;

    // TODO: use the correct Forward finite difference coefficients when on the edges?
    GxArray[offset] = - pixels[y * width + x0 + frameOffset] + pixels[y * width + x1 + frameOffset];
    GyArray[offset] = - pixels[y0 * width + x + frameOffset] + pixels[y1 * width + x + frameOffset];

    // Reset the local current frame
    nCurrentFrame[1] = 0;
}

// First kernel: evaluating gradient from image: 3-point gradient x --------------------------------------------------------------
__kernel void calculateGradient3pX(
    __global float* pixels,
    __global float* GxArray,
    __global float* GyArray,
    __global int* nCurrentFrame
    ) {
    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int f = get_global_id(2);

    const int frameOffset = wh*f;
    const int x0 = max(x-1, 0);
    const int x1 = min(x+1, width-1);
    const int y0 = max(y-1, 0);
    const int y1 = min(y+1, height-1);
    const int offset = y * width + x + frameOffset;

    GxArray[offset] = - pixels[y0 * width + x0 + frameOffset] - pixels[y1 * width + x0 + frameOffset] + pixels[y0 * width + x1 + frameOffset] + pixels[y1 * width + x1 + frameOffset];
    GyArray[offset] = - pixels[y0 * width + x0 + frameOffset] - pixels[y0 * width + x1 + frameOffset] + pixels[y1 * width + x0 + frameOffset] + pixels[y1 * width + x1 + frameOffset];

    // Reset the local current frame
    nCurrentFrame[1] = 0;
}


// First kernel: evaluating gradient from image using Robert's cross ------------------------------------------
__kernel void calculateGradientRobX(
    __global float* pixels,
    __global float* GxArray,
    __global float* GyArray,
    __global int* nCurrentFrame
    ) {
    const int x1 = get_global_id(0);
    const int y1 = get_global_id(1);
    const int f = get_global_id(2);

    const int frameOffset = wh*f;
    const int x0 = max(x1-1, 0);
    const int y0 = max(y1-1, 0);
    const int offset = y1 * width + x1 + frameOffset;

    // This calculates Robert's cross gradient and apply the rotation matrix 45 degrees to realign Gx and Gy to the image grid
    GxArray[offset] = pixels[y0 * width + x1 + frameOffset] - pixels[y1 * width + x0 + frameOffset] + pixels[y1 * width + x1 + frameOffset] - pixels[y0 * width + x0 + frameOffset];
    GyArray[offset] = - pixels[y0 * width + x1 + frameOffset] + pixels[y1 * width + x0 + frameOffset] + pixels[y1 * width + x1 + frameOffset] - pixels[y0 * width + x0 + frameOffset];

    // Reset the local current frame
    nCurrentFrame[1] = 0;
}

// First kernel: evaluating gradient from image using Robert's cross ------------------------------------------
__kernel void calculateGradient5pPlus(
    __global float* pixels,
    __global float* GxArray,
    __global float* GyArray,
    __global int* nCurrentFrame
    ) {
    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int f = get_global_id(2);

    const int frameOffset = wh*f;
//    const int x0 = max(x1-1, 0);
//    const int y0 = max(y1-1, 0);
    const int offset = y * width + x + frameOffset;

    // This calculates 5 points gradient using the finite difference coefficients
    if (x > 1 && x < width-2) GxArray[offset] = (1.0f/12.0f)*(pixels[y*width + x-2 + frameOffset] - 8.0f*pixels[y*width + x-1 + frameOffset] + 8.0f*pixels[y*width + x+1 + frameOffset] - pixels[y*width + x+2 + frameOffset]);
    else if (x == 1) GxArray[offset] = (1.0f/6.0f)*(-2.0f*pixels[y*width + x-1 + frameOffset] - 3.0f*pixels[y*width + x + frameOffset] + 6.0f*pixels[y*width + x+1 + frameOffset] - pixels[y*width + x+2 + frameOffset]);
    else if (x == 0) GxArray[offset] = (1.0f/2.0f)*(-3.0f*pixels[y*width + x + frameOffset] + 4.0f*pixels[y*width + x+1 + frameOffset] - pixels[y*width + x+2 + frameOffset]);
    else if (x == width-2) GxArray[offset] = (1.0f/6.0f)*(pixels[y*width + x-2 + frameOffset] - 6.0f*pixels[y*width + x-1 + frameOffset] + 3.0f*pixels[y*width + x + frameOffset] + 2.0f*pixels[y*width + x+1 + frameOffset]);
    else if (x == width-1) GxArray[offset] = (1.0f/2.0f)*(pixels[y*width + x-2 + frameOffset] - 4.0f*pixels[y*width + x-1 + frameOffset] + 3.0f*pixels[y*width + x + frameOffset]);

    if (y > 1 && y < height-2) GyArray[offset] = (1.0f/12.0f)*(pixels[(y-2)*width + x + frameOffset] - 8.0f*pixels[(y-1)*width + x + frameOffset] + 8.0f*pixels[(y+1)*width + x + frameOffset] - pixels[(y+2)*width + x + frameOffset]);
    else if (y == 1) GyArray[offset] = (1.0f/6.0f)*(-2.0f*pixels[(y-1)*width + x + frameOffset] - 3.0f*pixels[y*width + x + frameOffset] + 6.0f*pixels[(y+1)*width + x + frameOffset] - pixels[(y+2)*width + x + frameOffset]);
    else if (y == 0) GyArray[offset] = (1.0f/2.0f)*(-3.0f*pixels[y*width + x + frameOffset] + 4.0f*pixels[(y+1)*width + x + frameOffset] - pixels[(y+2)*width + x + frameOffset]);
    else if (y == height-2) GyArray[offset] = (1.0f/6.0f)*(pixels[(y-2)*width + x + frameOffset] - 6.0f*pixels[(y-1)*width + x + frameOffset] + 3.0f*pixels[y*width + x + frameOffset] + 2.0f*pixels[(y+1)*width + x + frameOffset]);
    else if (y == height-1) GyArray[offset] = (1.0f/2.0f)*(pixels[(y-2)*width + x + frameOffset] - 4.0f*pixels[(y-1)*width + x + frameOffset] + 3.0f*pixels[y*width + x + frameOffset]);

    // Reset the local current frame
    nCurrentFrame[1] = 0;
}


//// First kernel: evaluating gradient from image using 2-point gradient --------------------------------------------------
//__kernel void calculateGradient_2point(
//    __global float* pixels,
//    __global float* GxArray,
//    __global float* GyArray,
//    __global int* nCurrentFrame
//
//    ) {
//    const int x1 = get_global_id(0);
//    const int y1 = get_global_id(1);
//    const int f = get_global_id(2);
//
//    const int offset = y1 * w + x1 + w * h * f;
//    const int x0 = max(x1-1, 0);
//    const int y0 = max(y1-1, 0);
//
//    // 2-point gradient
//    GxArray[offset] = pixels[y1 * w + x1 + w * h * f] - pixels[y1 * w + x0 + w * h * f];
//    GyArray[offset] = pixels[y1 * w + x1 + w * h * f] - pixels[y0 * w + x1 + w * h * f];
//
//    // Reset the local current frame
//    nCurrentFrame[1] = 0;
////    if (nCurrentFrame[0] == nFrameForSRRF) nCurrentFrame[0] = 0; // reset the frame number if it's reached the end
//
//        // Current frame is a 2 element Int buffer:
//                // nCurrentFrame[0] is the global current frame in the current SRRF frame (reset every SRRF frame)
//                // nCurrentFrame[1] is the local current frame in the current GPU-loaded dataset (reset every turn of the method calculateSRRF (within the gradient calculation))
//
//}

//// First kernel follow up: interpolating gradient from 2-point gradient image -------------------------------------------
//__kernel void calculateGradientInterpolation(
//    __global float* GxArray,
//    __global float* GyArray,
//    __global float* GxIntArray,
//   __global float* GyIntArray
//
//    ){
//    const int x = get_global_id(0);
//    const int y = get_global_id(1);
//    const int f = get_global_id(2);
//
//    const int offset = y * wInt + x + f * wInt * hInt;
//
//    // Two-fold interpolation of the gradients
//    GxIntArray[offset] = getInterpolatedValue(GxArray, (int) (wInt/2), (int) (hInt/2), (float) (x)/2.0f, (float) (y)/2.0f, f);
//    GyIntArray[offset] = getInterpolatedValue(GyArray, (int) (wInt/2), (int) (hInt/2), (float) (x)/2.0f, (float) (y)/2.0f, f);
//}



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
//    const int yM = offset/wM;
//    const int xM = offset - yM*wM;
    const int yM = offset/width;
    const int xM = offset - yM*width;

    const float shiftX = shiftXY[nCurrentFrame[0]];
    const float shiftY = shiftXY[nCurrentFrame[0] + nFrameForSRRF];

    const float xc = (xM + 0.5) + shiftX; // continuous space position at the centre of magnified pixel, shiftXY need to be converted in unit of magnified space
    const float yc = (yM + 0.5) + shiftY;
//    const float sigma22 = 2 * sigma * sigma;

    float CGLH = 0; // CGLH stands for Radiality original name - Culley-Gustafsson-Laine-Henriques transform
    float distanceWeightSum = 0;

    float vx, vy, Gx, Gy;
//    float sigma = fwhm / 2.354f; // Sigma = 0.21 * lambda/NA in theory
//    float fradius = sigma * 2;
//    float radius = ((float) ((int) (GxGyMagnification*fradius)))/GxGyMagnification + 1;    // this reduces the radius for speed, works when using dGauss^4 and 2p+I
//    int radius = (int) (fradius) + 1;    // this should be used otherwise

    float dx, dy; // define dx and dy as the differences of coordinates
    float distance, distanceWeight, GdotR, GMag, Dk;

//// Estimate the mean gradient
//    float GxMean = 0;
//    float GyMean = 0;
//    int nGradient = 0;
//    for (int j=-(int) ((float) GxGyMagnification*(float) radius); j<=(int) ((float) GxGyMagnification*(float) radius+1); j++) { // stepping in pixels space of the Gradient matrix
//        vy = ((float) ((int) (GxGyMagnification*(yc-vxy_offset))) + j)/(float) GxGyMagnification + vxy_offset; // position in continuous space
//
//        if (vy > 0 && vy < h){
//            for (int i=-GxGyMagnification*radius; i<=(GxGyMagnification*radius+1); i++) {
//                vx = ((float) ((int) (GxGyMagnification*(xc-vxy_offset))) + i)/(float) GxGyMagnification + vxy_offset; // position in continuous space
//                if (vx > 0 && vx < w){
//
//                    dx = vx - xc;
//                    dy = vy - yc;
//                    distance = sqrt(dx*dx + dy*dy);    // Distance D
////                    float distance = sqrt((vx - xc)*(vx - xc) + (vy - yc)*(vy - yc));    // Distance D
//
////                    if (distance != 0 && distance <= (2*sigma+1)) {
//                    if (distance != 0 && distance <= (float) tS0) {
//
//                        nGradient ++;
//                        GxMean += getVBoundaryCheck(GxArray, wInt, hInt, GxGyMagnification*(vx - vxy_offset) + vxy_ArrayShift, GxGyMagnification*(vy - vxy_offset), nCurrentFrame[1]);
//                        GyMean += getVBoundaryCheck(GyArray, wInt, hInt, GxGyMagnification*(vx - vxy_offset), GxGyMagnification*(vy - vxy_offset) + vxy_ArrayShift, nCurrentFrame[1]);
//
//                    }
//                }
//            }
//        }
//    }
//    GxMean /= nGradient;
//    GyMean /= nGradient;


    for (int j=-(int) radius; j<=(int) (radius+1); j++) { // stepping in magnified pixel space
        vy = ((float) ((int) ((yc-vxy_offset))) + j) + vxy_offset; // position in continuous space

        if (vy > 0 && vy < height){
            for (int i=-(int) radius; i<=(int) (radius+1); i++) {
                vx = ((float) ((int) ((xc-vxy_offset))) + i) + vxy_offset; // position in continuous space

                if (vx > 0 && vx < width){

                    dx = vx - xc;
                    dy = vy - yc;
                    distance = sqrt(dx*dx + dy*dy);    // Distance D

//                    if (distance != 0 && distance <= (2*sigma+1)) {
                     if (distance != 0 && distance <= (float) tSO) {

                        Gx = getVBoundaryCheck(GxArray, width, height, vx - vxy_offset + vxy_ArrayShift, vy - vxy_offset, nCurrentFrame[1]);
                        Gy = getVBoundaryCheck(GyArray, width, height, vx - vxy_offset, vy - vxy_offset + vxy_ArrayShift, nCurrentFrame[1]);

//                     Gx = getVBoundaryCheck(GxArray, wInt, hInt, GxGyMagnification*(vx - vxy_offset) + vxy_ArrayShift, GxGyMagnification*(vy - vxy_offset), nCurrentFrame[1]) - GxMean;
//                     Gy = getVBoundaryCheck(GyArray, wInt, hInt, GxGyMagnification*(vx - vxy_offset), GxGyMagnification*(vy - vxy_offset) + vxy_ArrayShift, nCurrentFrame[1]) - GyMean;

                        distanceWeight = distance*exp(-(distance*distance)/((float) tSS));
                        distanceWeight = distanceWeight * distanceWeight * distanceWeight * distanceWeight; // dGauss^4
//distanceWeight = distanceWeight * distanceWeight;
                        distanceWeightSum += distanceWeight;
//distanceWeightSum ++;
                        GdotR = (Gx*dx + Gy*dy); // tells you if vector was pointing inward or outward

                        if (GdotR < 0) {
                            // Calculate perpendicular distance from (xc,yc) to gradient line through (vx,vy)
                            GMag = sqrt(Gx*Gx + Gy*Gy);
                            Dk = fabs(Gy*dx - Gx*dy) / GMag;    // Dk = D*sin(theta) obtained from cross-product
                            if (isnan(Dk)) Dk = distance; // this makes Dk = 0 in the next line

                            // Linear function ----------------------
                            Dk = 1 - Dk / distance; // Dk is now between 0 to 1, 1 if vector points precisely to (xc, yx), Dk = 1-sin(theta)

//                            Dk = GMag*(1 - Dk / distance); // for testing


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

//                        if (GdotR < 0) {
//                            CGLH += Dk*distanceWeight; // ----------------WEIGHTING----------------
////                            CGLH += Dk;
//                        }
//                        else {
//                            CGLH -= Dk*distanceWeight;
//                        }

                    CGLH += Dk*distanceWeight;
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
    }


//    CGLH = distanceWeightSum; // for testing
    CGLH /= distanceWeightSum; // ----------------WEIGHTING NORMALIZATION----------------
    if (CGLH >= 0) CGLH = pow(CGLH, sensitivity);
    else CGLH = 0;

    float v;
    if (shiftX == 0 && shiftY == 0) v = pixels[wh*nCurrentFrame[1] + yM*width + xM]; // no interpolation
    else v = getInterpolatedValue(pixels, width, height, xM + shiftX, yM + shiftY, nCurrentFrame[1]); // interpolation

    if (intWeighting == 1) {

            if (nCurrentFrame[0] == 0) { // Initialising the values of the array
                OutArray[offset] = v * CGLH / nFrameForSRRF;                    // AVG
                OutArray[offset + wh] = v * CGLH * v * CGLH / nFrameForSRRF;   // VAR
                OutArray[offset + 2 * wh] = 0;                                 // 2nd Order SOFI Tau = 1
                OutArray[offset + 3 * wh] = v / nFrameForSRRF;                 // interpolated
            }
            else {
                OutArray[offset] = OutArray[offset] + v * CGLH / nFrameForSRRF;
                OutArray[offset + wh] = OutArray[offset + wh] + v * CGLH * v * CGLH / nFrameForSRRF;
                OutArray[offset + 2 * wh] = OutArray[offset + 2 * wh] + v * CGLH * PreviousFrameArray[offset] / (nFrameForSRRF-1);
                OutArray[offset + 3 * wh] = OutArray[offset + 3 * wh] + v / nFrameForSRRF;
            }
        PreviousFrameArray[offset] = v * CGLH; // update the value of the previous frame

    }
    else{
            if (nCurrentFrame[0] == 0) {
                OutArray[offset] = CGLH / nFrameForSRRF;
                OutArray[offset + wh] = CGLH * CGLH / nFrameForSRRF;
                OutArray[offset + 2 * wh] = 0;
                OutArray[offset + 3 * wh] = v / nFrameForSRRF;
            }
            else {
                OutArray[offset] = OutArray[offset] + CGLH / nFrameForSRRF;
                OutArray[offset + wh] = OutArray[offset + wh] + CGLH * CGLH / nFrameForSRRF;
                OutArray[offset + 2 * wh] = OutArray[offset + 2 * wh] + CGLH * PreviousFrameArray[offset] / (nFrameForSRRF-1);
                OutArray[offset + 3 * wh] = OutArray[offset + 3 * wh] + v / nFrameForSRRF;
            }
        PreviousFrameArray[offset] = CGLH;
    }



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


// Kernel: calculate VARIANCE image from the OutputArray
__kernel void kernelCalculateVar(
    __global float* OutArray
    ){
    const int offset = get_global_id(0);

    const float av2 = OutArray[offset]*OutArray[offset];
    OutArray[offset + wh] = OutArray[offset + wh] - av2;     // Var[X] = E[X^2] - (E[X])^2
    OutArray[offset + 2*wh] = OutArray[offset + 2*wh] - av2; // TAC2[X] = E[Xt*Xt+1] - (E[X])^2

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

//// kernel: calculate the Macro-Pixel artefact map -----------------------------------------------------------------
//__kernel void kernelCalculateMPmap(
//    __global float* OutArray,
//    __global float* MPmap
//    ){
//
//    const int offset_MPmap = get_global_id(0);
//    const int frame = offset_MPmap/(magnification*magnification); // 0, 1 or 2 depending on whether we're dealing with AVG or VAR or SOFI2ndTau1
//    const int y_MPmap = (offset_MPmap - frame*(magnification*magnification))/magnification;
//    const int x_MPmap = offset_MPmap - y_MPmap*magnification - frame*(magnification*magnification);
//
//    float thisMPmapValue = 0;
//    int x;
//   int y;
//    int offset;
//
//    for (int i=0; i<wh; i++) {
//        y = i/w;
//        x = i - y*w;
//        offset = x_MPmap + x * magnification + wM * (y_MPmap + y * magnification) + frame * whM;
//        thisMPmapValue += OutArray[offset];
//    }
//    MPmap[offset_MPmap] = thisMPmapValue/wh;
//
//}

//// kernel: correct for Macro-pixel artefacts
//__kernel void kernelCorrectMPmap(
//     __global float* OutArray,
//     __global float* MPmap
//){
//
//    const int offset = get_global_id(0);
//    const int frame = offset/whM; // 0, 1 or 2 depending on whether we're dealing with AVG or VAR or SOFI2ndTau1
//    const int yM = (offset - frame*(whM))/wM;
//    const int xM = offset - yM*wM - frame*whM;
//    const int x = xM/magnification;
//    const int y = yM/magnification;
//
//    int offset_MPmap = magnification*magnification*frame + magnification*(yM - y*magnification) + xM - x*magnification;
//    OutArray[offset] = OutArray[offset]/MPmap[offset_MPmap];
//
//}
