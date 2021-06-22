//#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define magnification $MAGNIFICATION$
#define w $WIDTH$
#define h $HEIGHT$
#define wh $WH$
#define wM $WM$
#define hM $HM$
#define whM $WHM$


// kernel: calculate the Macro-Pixel artefact map -----------------------------------------------------------------
__kernel void kernelCalculateMPmap(
    __global float* OutArray,
    __global float* MPmap
    ){

    const int offset_MPmap = get_global_id(0);
    const int y_MPmap = offset_MPmap/magnification;
    const int x_MPmap = offset_MPmap - y_MPmap*magnification;

    float thisMPmapValue = 0;
    int x;
    int y;
    int offset;

    for (int i=0; i<wh; i++) {
        y = i/w;
        x = i - y*w;
        offset = x_MPmap + x * magnification + wM * (y_MPmap + y * magnification);
        thisMPmapValue += OutArray[offset];
    }
    MPmap[offset_MPmap] = thisMPmapValue/wh;

//    MPmap[offset_MPmap] = 0.5;

}

// kernel: correct for Macro-pixel artefacts
__kernel void kernelCorrectMPmap(
     __global float* OutArray,
     __global float* MPmap
){

    const int offset = get_global_id(0);
    const int yM = offset/wM;
    const int xM = offset - yM*wM;
    const int x = xM/magnification;
    const int y = yM/magnification;

    int offset_MPmap = magnification*(yM - y*magnification) + xM - x*magnification;
    OutArray[offset] = OutArray[offset]/MPmap[offset_MPmap];

}
