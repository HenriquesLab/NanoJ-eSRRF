__kernel void fhtsw(
    __global float* data,
    __global const float* sw,
    __global const float* cw,
    int w,
    int h,
    int d,
    int numSM, // number of SMs available on this device
    __local float* sData, //shared memory used to store the data
    __local float* sS,
    __local float* sC )
{
    int SM_Processor_NUM = get_global_id(0)/d;
    int lid = get_local_id(0);
    int l, i = 0;
    int temp = 0;
    int Ad1 = 0;
    int Ad2 = 0;
    int Ad3 = 0;
    int Ad4 = 0;
    int Ad0 = 0;
    int CSAd = 0;
    int numBfs, stage, gpNum, bfNum, gpSize, numGps;
    float rt1 = 0;    float rt2 = 0;    float rt3 = 0;    float rt4 = 0;

    int blockIndex = 0;
    int numFHT = 0;
    int maxN = 0;
    int Nlog2 = 0;

    int iIndex = 0;
    int jIndex = 0;
    int kIndex = 0;


    // set maxN
    maxN = w;

    if( lid < maxN )  //start in bounds check
    {
        // find Nlog2
        Nlog2 = findNlog2( maxN );

        if( lid < maxN/4 )
        {
            sC[lid] = cw[lid];
            sS[lid] = sw[lid];
        }


        for ( iIndex = SM_Processor_NUM; iIndex < d; iIndex+=numSM )
        {
            for ( jIndex = 0; jIndex < h; jIndex++ )
            {
                for (i = 0; i <= Nlog2; i++)
                {
                    if (( lid & ( 1 << i )) != 0 )
                    {
                        temp  |= ( 1 << ( Nlog2 - i - 1 ));
                    }
                }

                l = temp & 0x0000ffff;

                sData[ lid ] = data[ iIndex*h*w + jIndex*w + l];

                //Barrier
                barrier( CLK_LOCAL_MEM_FENCE );

                gpSize = 2;
                numGps = maxN / 4;
                if( lid < numGps )
                {
                    Ad1 = lid * 4;
                    Ad2 = Ad1 + 1;
                    Ad3 = Ad1 + gpSize;
                    Ad4 = Ad2 + gpSize;
                    rt1 = sData[Ad1] + sData[Ad2];   // a + b
                    rt2 = sData[Ad1] - sData[Ad2];   // a - b
                    rt3 = sData[Ad3] + sData[Ad4];   // c + d
                    rt4 = sData[Ad3] - sData[Ad4];   // c - d
                    sData[Ad1] = rt1 + rt3;      // a + b + (c + d)
                    sData[Ad2] = rt2 + rt4;      // a - b + (c - d)
                    sData[Ad3] = rt1 - rt3;      // a + b - (c + d)
                    sData[Ad4] = rt2 - rt4;      // a - b - (c - d)
                }

                //Barrier
                barrier( CLK_LOCAL_MEM_FENCE );

                if (Nlog2 > 2)
                {
                    // third + stages computed here
                    gpSize = 4;
                    numBfs = 2;
                    numGps = numGps / 2;

                    for (stage = 2; stage < Nlog2; stage++)
                    {
                        for (gpNum=0; gpNum<numGps; gpNum++)
                        {
                            Ad0 = gpNum * gpSize * 2;
                            if(lid==0)
                            {
                                Ad1 = Ad0;
                                Ad2 = Ad1 + gpSize;
                                Ad3 = Ad1 + gpSize / 2;
                                Ad4 = Ad3 + gpSize;

                                rt1 = sData[Ad1];
                                sData[Ad1] = sData[Ad1] + sData[Ad2];
                                sData[Ad2] = rt1 - sData[Ad2];
                                rt1 = sData[Ad3];
                                sData[Ad3] = sData[Ad3] + sData[Ad4];
                                sData[Ad4] = rt1 - sData[Ad4];
                            }

                            barrier( CLK_LOCAL_MEM_FENCE );

                            if( lid > 0 && lid < numBfs )
                            {
                                Ad1 = lid + Ad0;
                                Ad2 = Ad1 + gpSize;
                                Ad3 = gpSize - lid + Ad0;
                                Ad4 = Ad3 + gpSize;

                                CSAd = lid * numGps;
                                rt1 = sData[Ad2] * sC[CSAd] + sData[Ad4] * sS[CSAd];
                                rt2 = sData[Ad4] * sC[CSAd] - sData[Ad2] * sS[CSAd];

                                sData[Ad2] = sData[Ad1] - rt1;
                                sData[Ad1] = sData[Ad1] + rt1;
                                sData[Ad4] = sData[Ad3] + rt2;
                                sData[Ad3] = sData[Ad3] - rt2;
                            }

                        } // end gpNum loop

                        gpSize *= 2;
                        numBfs *= 2;
                        numGps = numGps / 2;

                    } // end for all stages
                } // end if Nlog2 > 2

                barrier( CLK_LOCAL_MEM_FENCE );

                data[ iIndex*h*w + jIndex*w + lid] = sData[lid];

            }  // end blockIndex
        }  // end rcIndex
    }  // stop in bounds check
};

__kernel void fhtsh(
    __global float* data,
    __global const float* sh,
    __global const float* ch,
    int w,
    int h,
    int d,
    int numSM, // number of SMs available on this device
    __local float* sData, //shared memory used to store the data
    __local float* sS,
    __local float* sC )
{
    int SM_Processor_NUM = get_global_id(0)/d;
    int lid = get_local_id(0);
    int l, i = 0;
    int temp = 0;
    int Ad1 = 0;
    int Ad2 = 0;
    int Ad3 = 0;
    int Ad4 = 0;
    int Ad0 = 0;
    int CSAd = 0;
    int numBfs, stage, gpNum, bfNum, gpSize, numGps;
    float rt1 = 0;    float rt2 = 0;    float rt3 = 0;    float rt4 = 0;

    int blockIndex = 0;
    int numFHT = 0;
    int maxN = 0;
    int Nlog2 = 0;

    int iIndex = 0;
    int jIndex = 0;
    int kIndex = 0;

    // set maxN
    maxN = h;

    if( lid < maxN )  //start in bounds check
    {
        // find Nlog2
        Nlog2 = findNlog2( maxN );

        if( lid < maxN/4 )
        {
            sC[lid] = ch[lid];
            sS[lid] = sh[lid];
        }


        for ( iIndex = SM_Processor_NUM; iIndex < d; iIndex+=numSM )
        {
            for ( jIndex = 0; jIndex < h; jIndex++ )
            {
                for (i = 0; i <= Nlog2; i++)
                {
                    if (( lid & ( 1 << i )) != 0 )
                    {
                        temp  |= ( 1 << ( Nlog2 - i - 1 ));
                    }
                }

                l = temp & 0x0000ffff;

                sData[ lid ] = data[ iIndex*h*w + l*w + jIndex ];

                //Barrier
                barrier( CLK_LOCAL_MEM_FENCE );

                gpSize = 2;
                numGps = maxN / 4;
                if( lid < numGps )
                {
                    Ad1 = lid * 4;
                    Ad2 = Ad1 + 1;
                    Ad3 = Ad1 + gpSize;
                    Ad4 = Ad2 + gpSize;
                    rt1 = sData[Ad1] + sData[Ad2];   // a + b
                    rt2 = sData[Ad1] - sData[Ad2];   // a - b
                    rt3 = sData[Ad3] + sData[Ad4];   // c + d
                    rt4 = sData[Ad3] - sData[Ad4];   // c - d
                    sData[Ad1] = rt1 + rt3;      // a + b + (c + d)
                    sData[Ad2] = rt2 + rt4;      // a - b + (c - d)
                    sData[Ad3] = rt1 - rt3;      // a + b - (c + d)
                    sData[Ad4] = rt2 - rt4;      // a - b - (c - d)
                }

                //Barrier
                barrier( CLK_LOCAL_MEM_FENCE );

                if (Nlog2 > 2)
                {
                    // third + stages computed here
                    gpSize = 4;
                    numBfs = 2;
                    numGps = numGps / 2;

                    for (stage = 2; stage < Nlog2; stage++)
                    {
                        for (gpNum=0; gpNum<numGps; gpNum++)
                        {
                            Ad0 = gpNum * gpSize * 2;
                            if(lid==0)
                            {
                                Ad1 = Ad0;
                                Ad2 = Ad1 + gpSize;
                                Ad3 = Ad1 + gpSize / 2;
                                Ad4 = Ad3 + gpSize;

                                rt1 = sData[Ad1];
                                sData[Ad1] = sData[Ad1] + sData[Ad2];
                                sData[Ad2] = rt1 - sData[Ad2];
                                rt1 = sData[Ad3];
                                sData[Ad3] = sData[Ad3] + sData[Ad4];
                                sData[Ad4] = rt1 - sData[Ad4];
                            }

                            barrier( CLK_LOCAL_MEM_FENCE );

                            if( lid > 0 && lid < numBfs )
                            {
                                Ad1 = lid + Ad0;
                                Ad2 = Ad1 + gpSize;
                                Ad3 = gpSize - lid + Ad0;
                                Ad4 = Ad3 + gpSize;

                                CSAd = lid * numGps;
                                rt1 = sData[Ad2] * sC[CSAd] + sData[Ad4] * sS[CSAd];
                                rt2 = sData[Ad4] * sC[CSAd] - sData[Ad2] * sS[CSAd];

                                sData[Ad2] = sData[Ad1] - rt1;
                                sData[Ad1] = sData[Ad1] + rt1;
                                sData[Ad4] = sData[Ad3] + rt2;
                                sData[Ad3] = sData[Ad3] - rt2;
                            }

                        } // end gpNum loop

                        gpSize *= 2;
                        numBfs *= 2;
                        numGps = numGps / 2;

                    } // end for all stages
                } // end if Nlog2 > 2

                barrier( CLK_LOCAL_MEM_FENCE );

                data[ iIndex*h*w + lid*w + jIndex ] = sData[lid];

            }  // end blockIndex
        }  // end rcIndex
    }  // stop in bounds check
};

__kernel void fhts(
    __global float* data,
    __global const float* s,
    __global const float* c,
    int w,
    int h,
    int d,
    int numSM, // number of SMs available on this device
    __local float* sData, //shared memory used to store the data
    __local float* sS,
    __local float* sC )
{
    int SM_Processor_NUM = get_global_id(0)/d;
    int lid = get_local_id(0);
    int l, i = 0;
    int temp = 0;
    int Ad1 = 0;
    int Ad2 = 0;
    int Ad3 = 0;
    int Ad4 = 0;
    int Ad0 = 0;
    int CSAd = 0;
    int numBfs, stage, gpNum, bfNum, gpSize, numGps;
    float rt1 = 0;    float rt2 = 0;    float rt3 = 0;    float rt4 = 0;

    int blockIndex = 0;
    int numFHT = 0;
    int maxN = 0;
    int Nlog2 = 0;

    int iIndex = 0;
    int jIndex = 0;
    int kIndex = 0;

    // set maxN
    maxN = d;

    if( lid < maxN )  //start in bounds check
    {
        // find Nlog2
        Nlog2 = findNlog2( maxN );

        if( lid < maxN/4 )
        {
            sC[lid] = c[lid];
            sS[lid] = s[lid];
        }

        for ( iIndex = SM_Processor_NUM; iIndex < d; iIndex+=numSM )
        {
            for ( jIndex = 0; jIndex < h; jIndex++ )
            {
                for (i = 0; i <= Nlog2; i++)
                {
                    if (( lid & ( 1 << i )) != 0 )
                    {
                        temp  |= ( 1 << ( Nlog2 - i - 1 ));
                    }
                }

                l = temp & 0x0000ffff;

                sData[ lid ] = data[ iIndex*w + l*w*h + jIndex ];

                //Barrier
                barrier( CLK_LOCAL_MEM_FENCE );

                gpSize = 2;
                numGps = maxN / 4;
                if( lid < numGps )
                {
                    Ad1 = lid * 4;
                    Ad2 = Ad1 + 1;
                    Ad3 = Ad1 + gpSize;
                    Ad4 = Ad2 + gpSize;
                    rt1 = sData[Ad1] + sData[Ad2];   // a + b
                    rt2 = sData[Ad1] - sData[Ad2];   // a - b
                    rt3 = sData[Ad3] + sData[Ad4];   // c + d
                    rt4 = sData[Ad3] - sData[Ad4];   // c - d
                    sData[Ad1] = rt1 + rt3;      // a + b + (c + d)
                    sData[Ad2] = rt2 + rt4;      // a - b + (c - d)
                    sData[Ad3] = rt1 - rt3;      // a + b - (c + d)
                    sData[Ad4] = rt2 - rt4;      // a - b - (c - d)
                }

                //Barrier
                barrier( CLK_LOCAL_MEM_FENCE );

                if (Nlog2 > 2)
                {
                    // third + stages computed here
                    gpSize = 4;
                    numBfs = 2;
                    numGps = numGps / 2;

                    for (stage = 2; stage < Nlog2; stage++)
                    {
                        for (gpNum=0; gpNum<numGps; gpNum++)
                        {
                            Ad0 = gpNum * gpSize * 2;
                            if(lid==0)
                            {
                                Ad1 = Ad0;
                                Ad2 = Ad1 + gpSize;
                                Ad3 = Ad1 + gpSize / 2;
                                Ad4 = Ad3 + gpSize;

                                rt1 = sData[Ad1];
                                sData[Ad1] = sData[Ad1] + sData[Ad2];
                                sData[Ad2] = rt1 - sData[Ad2];
                                rt1 = sData[Ad3];
                                sData[Ad3] = sData[Ad3] + sData[Ad4];
                                sData[Ad4] = rt1 - sData[Ad4];
                            }

                            barrier( CLK_LOCAL_MEM_FENCE );

                            if( lid > 0 && lid < numBfs )
                            {
                                Ad1 = lid + Ad0;
                                Ad2 = Ad1 + gpSize;
                                Ad3 = gpSize - lid + Ad0;
                                Ad4 = Ad3 + gpSize;

                                CSAd = lid * numGps;
                                rt1 = sData[Ad2] * sC[CSAd] + sData[Ad4] * sS[CSAd];
                                rt2 = sData[Ad4] * sC[CSAd] - sData[Ad2] * sS[CSAd];

                                sData[Ad2] = sData[Ad1] - rt1;
                                sData[Ad1] = sData[Ad1] + rt1;
                                sData[Ad4] = sData[Ad3] + rt2;
                                sData[Ad3] = sData[Ad3] - rt2;
                            }

                        } // end gpNum loop

                        gpSize *= 2;
                        numBfs *= 2;
                        numGps = numGps / 2;

                    } // end for all stages
                } // end if Nlog2 > 2

                barrier( CLK_LOCAL_MEM_FENCE );

                data[ iIndex*w + lid*w*h + jIndex ] = sData[lid];

            }  // end blockIndex
        }  // end rcIndex
    }  // stop in bounds check

};

__kernel void fhtf(
    __global float* data,
    int w,
    int h,
    int d,
    int numSM, // number of SMs available on this device
    int inverse )
{
    int SM_Processor_NUM = get_global_id(0)/d;
    int lid = get_local_id(0);
    int l, i = 0;
    int temp = 0;
    int Ad1 = 0;
    int Ad2 = 0;
    int Ad3 = 0;
    int Ad4 = 0;
    int Ad0 = 0;
    int CSAd = 0;
    int numBfs, stage, gpNum, bfNum, gpSize, numGps;
    float rt1 = 0;    float rt2 = 0;    float rt3 = 0;    float rt4 = 0;
    float rt5 = 0;    float rt6 = 0;    float rt7 = 0;    float rt8 = 0;

    int blockIndex = 0;
    int numFHT = 0;
    int maxN = 0;
    int Nlog2 = 0;

    int iIndex = 0;
    int jIndex = 0;
    int kIndex = 0;

    int whd2 = 2*w*h*d;

    // set maxN
    maxN = w;

    //Convert to actual Hartley transform
    for(iIndex = SM_Processor_NUM; iIndex <= d/2; iIndex+=numSM)
    {
        Ad3 = (d - iIndex) % d;
        for(jIndex = 0; jIndex <= h/2; jIndex++)
        {
            Ad2 = (h - jIndex) % h;
            if ( lid <= w/2 )
            {
                Ad1 = (w - lid) % w;
                rt1 = data[ iIndex*w*h + lid + w*Ad2];
                rt2 = data[ iIndex*w*h + Ad1 + w*jIndex];
                rt3 = data[ Ad3*w*h + lid + w*jIndex];
                rt4 = data[ Ad3*w*h + Ad1 + w*Ad2];
                rt5 = data[ Ad3*w*h + lid + w*Ad2];
                rt6 = data[ Ad3*w*h + Ad1 + w*jIndex];
                rt7 = data[ iIndex*w*h + lid + w*jIndex];
                rt8 = data[ iIndex*w*h + Ad1 + w*Ad2];
                if(inverse==0)
                {
                    data[ iIndex*w*h + lid + w*jIndex] = (rt1+rt2+rt3-rt4)/2;
                    data[ Ad3*w*h + lid + w*jIndex] = (rt5+rt6+rt7-rt8)/2;
                    data[ iIndex*w*h + lid + w*Ad2] = (rt7+rt8+rt5-rt6)/2;
                    data[ Ad3*w*h + lid + w*Ad2] = (rt3+rt4+rt1-rt2)/2;
                    data[ iIndex*w*h + Ad1 + w*jIndex] = (rt8+rt7+rt6-rt5)/2;
                    data[ Ad3*w*h + Ad1 + w*jIndex] = (rt4+rt3+rt2-rt1)/2;
                    data[ iIndex*w*h + Ad1 + w*Ad2] = (rt2+rt1+rt4-rt3)/2;
                    data[ Ad3*w*h + Ad1 + w*Ad2] = (rt6+rt5+rt8-rt7)/2;
                } else {
                    data[ iIndex*w*h + lid + w*jIndex] = (rt1+rt2+rt3-rt4)/whd2;
                    data[ Ad3*w*h + lid + w*jIndex] = (rt5+rt6+rt7-rt8)/whd2;
                    data[ iIndex*w*h + lid + w*Ad2] = (rt7+rt8+rt5-rt6)/whd2;
                    data[ Ad3*w*h + lid + w*Ad2] = (rt3+rt4+rt1-rt2)/whd2;
                    data[ iIndex*w*h + Ad1 + w*jIndex] = (rt8+rt7+rt6-rt5)/whd2;
                    data[ Ad3*w*h + Ad1 + w*jIndex] = (rt4+rt3+rt2-rt1)/whd2;
                    data[ iIndex*w*h + Ad1 + w*Ad2] = (rt2+rt1+rt4-rt3)/whd2;
                    data[ Ad3*w*h + Ad1 + w*Ad2] = (rt6+rt5+rt8-rt7)/whd2;
                }
            }
        }
    }
};

int findNlog2( int maxN )
{
    int Nlog2 = 15;
    while (!((maxN & (1<<Nlog2)) != 0))
    {
        Nlog2--;
    }

    return Nlog2;
}