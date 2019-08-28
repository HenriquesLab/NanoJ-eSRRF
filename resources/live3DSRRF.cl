//#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define magnification $MAGNIFICATION$
#define fwhm $FWHM$
#define sensitivity $SENSITIVITY$
#define GxGyMagnification $GXGYMAGNIFICATION$
#define sigma $SIGMA$
#define radius $RADIUS$
#define w $WIDTH$ // these are width and height of the individual split images
#define h $HEIGHT$
#define wh $WH$
#define wInt $WINT$
#define hInt $HINT$
#define whInt $WHINT$
#define wM $WM$
#define hM $HM$
#define whM $WHM$
#define whdM $WHDM$
#define nFrameForSRRF $NFRAMEFORSRRF$
#define intWeighting $INTWEIGHTING$

//#define wS $WIDTHTHREED$
//#define hS $HEIGHTTHREED$
//#define whS $WHS$
//#define wSInt $WSINT$
//#define hSInt $HSINT$
//#define whSInt $WHSINT$
#define nPlanes $NPLANES$ // TODO: necessary when doing differnt number of planes! Big bug ahead with different pixels between full frame and sum of crops splits
#define dimRatio $DIMRATIO$
#define whd $WHD$ // w * h * nPlanes

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
static float getInterpolatedValue(__global float* array, float const x, float const y, int const f, int const z) {
    const int u0 = (int) floor(x);
    const int v0 = (int) floor(y);
    const int fdOffset = whd*f + wh*z;

    float q = 0.0f;

    // Bicubic interpolation
    if (u0 > 0 && u0 < w - 2 && v0 > 0 && v0 < h - 2) {
        for (int j = 0; j <= 3; j++) {
            int v = min(max(v0 - 1 + j, 0), h-1);
            float p = 0.0f;
            for (int i = 0; i <= 3; i++) {
                int u = min(max(u0 - 1 + i, 0), w-1);
                p = p + array[u + v*w + fdOffset] * cubic(x - (float) (u));
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

        int xbase = (int) fmin((float) w-2, fmax(x,0.0f));
        int xbase1 = xbase+1;

        int ybase = (int) fmin((float) h-2, fmax(y,0.0f));
        int ybase1 = ybase+1;

        float xFraction = x - (float) xbase;
        float yFraction = y - (float) ybase;
//        xFraction = fmax(xFraction, 0);
//        yFraction = fmax(yFraction, 0);

        float lowerLeft = array[fdOffset + ybase * w + xbase];
        float lowerRight = array[fdOffset + ybase * w + xbase1];
        float upperRight = array[fdOffset + ybase1 * w + xbase1];
        float upperLeft = array[fdOffset + ybase1 * w + xbase];
        float upperAverage = upperLeft + xFraction * (upperRight - upperLeft);
        float lowerAverage = lowerLeft + xFraction * (lowerRight - lowerLeft);
        q = lowerAverage + yFraction * (upperAverage - lowerAverage);

    }

    return q;
}

// Check boundaries of the image and returns the gradient value // TODO: extrapolate instead of boundary check?
static float getVBoundaryCheck(__global float* array, int const thisWidth, int const thisHeight, int const thisDepth, int const x, int const y, int const z, int const f) {
    const int _x = min(max(x, 0), thisWidth-1);
    const int _y = min(max(y, 0), thisHeight-1);
    const int _z = min(max(z, 0), thisDepth-1);
    return array[thisDepth*thisWidth*thisHeight*f + thisWidth*thisHeight*_z + _y*thisWidth + _x];
}

// Calculate the amplitude of the cross-product between two 3D vectors
static float getCrossProductMagnitude(float const a1, float const a2, float const a3, float const b1, float const b2, float const b3){
    const float G1 = a2*b3 - a3*b2;
    const float G2 = a3*b1 - a1*b3;
    const float G3 = a1*b2 - a2*b1;
    const float CP = sqrt(G1*G1 + G2*G2 + G3*G3);
    return CP;
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
    __global float* GzArray,
    __global int* nCurrentFrame

    ) {
    const int offset = get_global_id(0);

    const int f = offset/whd;
    const int z1 = (offset - f*whd)/wh;
    const int y1 = (offset - f*whd - z1*wh)/w;
    const int x1 =  offset - f*whd - z1*wh - y1*w;

    const int x0 = max(x1-1, 0); // this sets the gradient to zero on the edges, assumes same value adjacent
    const int y0 = max(y1-1, 0);
    const int z0 = max(z1-1, 0);

    // 2-point gradient
    GxArray[offset] = pixels[offset] - pixels[x0 + y1*w + z1*wh + whd*f];
    GyArray[offset] = pixels[offset] - pixels[x1 + y0*w + z1*wh + whd*f];
    GzArray[offset] = pixels[offset] - pixels[x1 + y1*w + z0*wh + whd*f];

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
    __global float* GzArray,
    __global float* GxIntArray,
    __global float* GyIntArray,
    __global float* GzIntArray

    ){

    const int offset = get_global_id(0);

    const int whdInt = wInt*hInt*nPlanes; // TODO: add as #define?
    const int f = offset/whdInt;
    const int z = (offset - f*whdInt)/whInt;
    const int y = (offset - f*whdInt - z*whInt)/wInt;
    const int x =  offset - f*whdInt - z*whInt - y*wInt;

//    const int x = get_global_id(0);
//    const int y = get_global_id(1);
//    const int f = get_global_id(2);
//
//    const int offset = y * wInt + x + f * wInt * hInt;

    // LATERAL interpolation of the gradients
    GxIntArray[offset] = getInterpolatedValue(GxArray, (float) (x)/2.0f, (float) (y)/2.0f, f, z);
    GyIntArray[offset] = getInterpolatedValue(GyArray, (float) (x)/2.0f, (float) (y)/2.0f, f, z);
    GzIntArray[offset] = getInterpolatedValue(GzArray, (float) (x)/2.0f, (float) (y)/2.0f, f, z);
}



// Second kernel: evaluating RGC from gradient -----------------------------------------------------------------
__kernel void calculateRadialGradientConvergence(
    __global float* pixels,
    __global float* GxArray, // these do not have the Int in the name but are the interpolated versions
    __global float* GyArray,
    __global float* GzArray,
    __global float* OutArray,
    __global float* driftXY,
    __global float* shiftXY3D,
    __global int* nCurrentFrame
    // Current frame is a 2 element Int buffer:
            // nCurrentFrame[0] is the global current frame in the current SRRF frame (reset every SRRF frame)
            // nCurrentFrame[1] is the local current frame in the current GPU-loaded dataset (reset every turn of the method calculateSRRF (within the gradient calculation))


    ) {

    const int offset = get_global_id(0);

    const int f = offset/(whdM);
    const int zM = (offset - f*whdM)/whM;
//    const int z = (offset - f*whdM)/whM; // TODO: work out the z position for the shift and implement
    const int yM = (offset - f*whdM - zM*whM)/wM;
    const int xM =  offset - f*whdM - zM*whM - yM*wM;

    const float driftX = driftXY[nCurrentFrame[0]];
    const float driftY = driftXY[nCurrentFrame[0] + nFrameForSRRF];
    const float sigma21 = 2 * sigma + 1;

    // Coordinates of the pixel of interest in continuous space
    const float xc = (xM + 0.5) / magnification + driftX; // continuous space position at the centre of magnified pixel
    const float yc = (yM + 0.5) / magnification + driftY;
    const float zc = (zM + 0.5) / magnification; // TODO: consider drift in Z?
    const float sigSquare = 2 * sigma * sigma; // TODO: add as hardcoded value?

    float CGLH = 0; // CGLH stands for Radiality original name - Culley-Gustafsson-Laine-Henriques transform
    float distanceWeightSum = 0;

    float vx, vy, vz, Gx, Gy, Gz;
//    float sigma = fwhm / 2.354f; // Sigma = 0.21 * lambda/NA in theory
//    float fradius = sigma * 2;
//    float radius = ((float) ((int) (GxGyMagnification*fradius)))/GxGyMagnification + 1;    // this reduces the radius for speed, works when using dGauss^4 and 2p+I
//    int radius = (int) (fradius) + 1;    // this should be used otherwise

    for (int j=-GxGyMagnification*radius; j<=(GxGyMagnification*radius+1); j++) {
        vy = ((float) ((int) (GxGyMagnification*yc)) + j)/GxGyMagnification; // position in continuous space TODO: problems with negative values and (int)?

        if (vy > 0 && vy < h){
            for (int i=-GxGyMagnification*radius; i<=(GxGyMagnification*radius+1); i++) {
                vx = ((float) ((int) (GxGyMagnification*xc)) + i)/GxGyMagnification; // position in continuous space TODO: problems with negative values and (int)?
                if (vx > 0 && vx < w){

                    for (int k=-radius; k<=(radius); k++) {
                        vz = ((float) ((int) (zc)) + k); // position in continuous space TODO: problems with negative values and (int)?
                        if (vz > 0 && vz < nPlanes){

                            float distance = sqrt((vx - xc)*(vx - xc) + (vy - yc)*(vy - yc) + (vz - zc)*(vz - zc));    // Distance D

                            if (distance != 0 && distance <= sigma21) {

                                Gx = getVBoundaryCheck(GxArray, wInt, hInt, nPlanes, GxGyMagnification*(vx - vxy_offset) + vxy_ArrayShift, GxGyMagnification*(vy - vxy_offset), vz, nCurrentFrame[1]);
                                Gy = getVBoundaryCheck(GyArray, wInt, hInt, nPlanes, GxGyMagnification*(vx - vxy_offset), GxGyMagnification*(vy - vxy_offset) + vxy_ArrayShift, vz, nCurrentFrame[1]);
                                Gz = getVBoundaryCheck(GzArray, wInt, hInt, nPlanes, GxGyMagnification*(vx - vxy_offset), GxGyMagnification*(vy - vxy_offset), vz + vxy_ArrayShift, nCurrentFrame[1]);

                                float distanceWeight = distance*exp(-(distance*distance)/sigSquare);  // TODO: dGauss: can use Taylor expansion there
                                distanceWeight = distanceWeight * distanceWeight * distanceWeight * distanceWeight ;
                                distanceWeightSum += distanceWeight;
                                float GdotR = (Gx * (vx - xc) + Gy * (vy - yc) + Gz * (vz - zc)); // tells you if vector was pointing inward or outward

                                if (GdotR < 0) {
                                    // Calculate perpendicular distance from (xc,yc) to gradient line through (vx,vy)
                                    float GMag = sqrt(Gx * Gx + Gy * Gy + Gz * Gz);
                                    float Dk = getCrossProductMagnitude(xc-vx, yc-vy, zc-vz, Gx, Gy, Gz)/GMag; // Dk = D*sin(theta) obtained from cross-product
//                                  float Dk = fabs(Gy * (xc - vx) - Gx * (yc - vy)) / GMag;    // Dk = D*sin(theta) obtained from cross-product
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
                            }
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

    float v = getInterpolatedValue(pixels, ((float) xM)/magnification + driftX - 0.5f, ((float) yM)/magnification + driftY - 0.5f, nCurrentFrame[1], zM);

    if (intWeighting == 1) {
            if (nCurrentFrame[0] == 0) {
                OutArray[offset] = v * CGLH / nFrameForSRRF;
                OutArray[offset + whM] = v * CGLH * v * CGLH / nFrameForSRRF;
                OutArray[offset + 2 * whM] = v / nFrameForSRRF;
            }
            else {
                OutArray[offset] = OutArray[offset] + v * CGLH / nFrameForSRRF;
                OutArray[offset + whM] = OutArray[offset + whM] + v * CGLH * v * CGLH / nFrameForSRRF;
                OutArray[offset + 2 * whM] = OutArray[offset + 2 * whM] + v / nFrameForSRRF;
            }
    }
    else{
            if (nCurrentFrame[0] == 0) {
                OutArray[offset] = CGLH / nFrameForSRRF;
                OutArray[offset + whM] = CGLH * CGLH / nFrameForSRRF;
                OutArray[offset + 2 * whM] = v / nFrameForSRRF;
            }
            else {
                OutArray[offset] = OutArray[offset] + CGLH / nFrameForSRRF;
                OutArray[offset + whM] = OutArray[offset + whM] + CGLH * CGLH / nFrameForSRRF;
                OutArray[offset + 2 * whM] = OutArray[offset + 2 * whM] + v / nFrameForSRRF;
            }
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


// Kernel: calculate STD image from the OutputArray
__kernel void kernelCalculateStd(
    __global float* OutArray
    ){

    const int offset = get_global_id(0);
    OutArray[offset + whM] = sqrt(OutArray[offset + whM] - OutArray[offset]*OutArray[offset]);

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
    const int frame = offset_MPmap/(magnification*magnification); // 0 or 1 depending on whether we're dealing with AVG or STD
    const int y_MPmap = (offset_MPmap - frame*(magnification*magnification))/magnification;
    const int x_MPmap = offset_MPmap - y_MPmap*magnification - frame*(magnification*magnification);

    float thisMPmapValue = 0;
    int x;
    int y;
    int offset;

    for (int i=0; i<wh; i++) {
        y = i/w;
        x = i - y*w;
        offset = x_MPmap + x * magnification + wM * (y_MPmap + y * magnification) + frame * whM;
        thisMPmapValue += OutArray[offset];
    }
    MPmap[offset_MPmap] = thisMPmapValue/wh;

}

// kernel: correct for Macro-pixel artefacts
__kernel void kernelCorrectMPmap(
     __global float* OutArray,
     __global float* MPmap
){

    const int offset = get_global_id(0);
    const int frame = offset/whM; // 0 or 1 depending on whether we're dealing with AVG or STD
    const int yM = (offset - frame*(whM))/wM;
    const int xM = offset - yM*wM - frame*whM;
    const int x = xM/magnification;
    const int y = yM/magnification;

    int offset_MPmap = magnification*magnification*frame + magnification*(yM - y*magnification) + xM - x*magnification;
    OutArray[offset] = OutArray[offset]/MPmap[offset_MPmap];

}
