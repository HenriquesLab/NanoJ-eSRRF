//#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define nLocs $NLOCS$

// kernel: calculate the NND^2
__kernel void calculateNND2(
     __global float* xyBuffer,
     __global float* NNDbuffer
){

    const int n = get_global_id(0);

    float dx2;
    float dy2;
    float thisD;
    float minD = 1000000000; // TODO: a really big number

    for (int i=0; i<nLocs; i++) {
        if (i != n){
                dx2 = (xyBuffer[i]-xyBuffer[n])*(xyBuffer[i]-xyBuffer[n]);
                if (dx2 < minD){
                    dy2 = (xyBuffer[i+nLocs]-xyBuffer[n+nLocs])*(xyBuffer[i+nLocs]-xyBuffer[n+nLocs]);
                    thisD = dx2 +dy2;
                    if ((thisD) < minD) minD = thisD;
                }
        }
    }

    NNDbuffer[n] = minD;

}


// kernel: calculate the SQRT of the NND^2 for each localizations
__kernel void calculateSqrtNND2(
     __global float* NNDbuffer
){
    const int n = get_global_id(0);
    NNDbuffer[n] = sqrt(NNDbuffer[n]);

}