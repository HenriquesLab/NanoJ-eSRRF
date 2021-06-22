//#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define nLocs $NLOCS$
#define densityRadius2 $DENSITYRADIUS2$

// kernel: calculate the NND^2
__kernel void calculateNND2(
     __global float* xyBuffer,
     __global float* NNDbuffer,
     __global int* densityBuffer
){

    const int n = get_global_id(0);

    float dx2;
    float thisD;
    float minD;
    bool isFirst = true;
    int nLocInArea = 0;


    for (int i=0; i<nLocs; i++) {
        if (i != n){
            dx2 = (xyBuffer[i]-xyBuffer[n])*(xyBuffer[i]-xyBuffer[n]);
            if (isFirst){
                minD = dx2 + (xyBuffer[i+nLocs]-xyBuffer[n+nLocs])*(xyBuffer[i+nLocs]-xyBuffer[n+nLocs]);
                if (minD < densityRadius2) nLocInArea++;
                isFirst = false;
            }
            else {
                if (dx2 < minD | dx2 < densityRadius2){
                    thisD = dx2 + (xyBuffer[i+nLocs]-xyBuffer[n+nLocs])*(xyBuffer[i+nLocs]-xyBuffer[n+nLocs]);
                    if (thisD < minD) minD = thisD;
                    if (thisD < densityRadius2) nLocInArea++;
                }
            }

        }
    }

    NNDbuffer[n] = minD;
    densityBuffer[n] = nLocInArea;

}



// kernel: calculate the SQRT of the NND^2 for each localizations
__kernel void calculateSqrtNND2(
     __global float* NNDbuffer
){
    const int n = get_global_id(0);
    NNDbuffer[n] = sqrt(NNDbuffer[n]);

}