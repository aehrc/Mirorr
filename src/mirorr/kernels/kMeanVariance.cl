#ifndef __OPENCL_VERSION__ //defined by the kernel compiler
#define __kernel
#define __global
#define CLK_LOCAL_MEM_FENCE
#define uint unsigned int
#endif

/**
 *
 * --- WARNING ---
 * The .h version of this file is generated from the .cl file. You
 * MUST run convertCl2Header.sh after any modification to the .cl file!
 * --- WARNING ----
 *
 *
 * Code by Jeremy Coatelen, 11 Apr 2011
 * Commented and edited by D Rivest-Henault, 31 Jan 2013
 * 
 * * Copyright (c) 2009-15 CSIRO. All rights reserved.
 * 
 * For complete copyright, license and disclaimer of warranty
 * information see the LICENSE.txt file for details.
 * 
 * This software is distributed WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the above copyright notice for more information.
 * 
 */


/**
 * kMeanVariance:
 * This OpenCL kernel takes two image buffers and compute the mean
 * and variance of the BlockWidth^N, N in{2,3}, with starting index
 * defined by "ig = get_global_id(0)".
 *
 * @param inputL Base image buffer
 * @param inputR Search image buffer
 * @param outputLMean Output mean buffer for the inputL image
 * @param outputRMean Output mean buffer for the inputR image
 * @param outputLVariance Output variance buffer for the inputL image
 * @param outputRVariance Output variance buffer for the inputR image
 * @param w Width of the search image buffer
 * @param h Height of the search image buffer
 * @param d Depth of the search image buffer
 * @param BlockWidth size of the {2,3}D block
 */

//__attribute__((vec_type_hint(float4)))
//#pragma OPENCL EXTENSION cl_amd_printf : enable
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
__kernel
void kMeanVariance(
   __global float * inputL,    __global float * inputR,
   __global float * outputLMean,     __global float * outputRMean,
   __global float * outputLVariance,   __global float * outputRVariance,
   const uint w, const uint h, const uint d,
   uint BlockWidth)
{
  int ig = get_global_id(0);
#if 0
  float   max1  = 0.0; float max2  = 0.0;
  float   mean1 = 0.0; float mean2 = 0.0;
  float   normalisation;
#else
  double   max1  = 0.0; double max2  = 0.0;
  double   mean1 = 0.0; double mean2 = 0.0;
  double   normalisation;
#endif
  int     i = 0; int j = 0; int k = 0;
  int     I; int J; int K;
  int     wh = w * h;
  float   diff;
  int     tmp; int n;
  int     BlockWidth2 = BlockWidth*BlockWidth;
  int     BlockWidth3 = BlockWidth2*BlockWidth;

  if(ig < wh*d)
  {
    // Initialization of the matrices (edges)
    normalisation = (float)((d>1)?(BlockWidth3):(BlockWidth2));\
    outputLMean[ig] = 1.0/normalisation;
    outputRMean[ig] = 1.0/normalisation;
    outputLVariance[ig] = 1.0/normalisation;
    outputRVariance[ig] = 1.0/normalisation;

    // Set the indexes
    i = ig % wh;
    I = i % w;      //  I in [ 0 ... WIDTH  ]
    J = i / w;      //  J in [ 0 ... HEIGHT ]
    K = ig / wh;    //  K in [ 0 ... DEPTH  ]

    // Check the edges
    if( I > w-BlockWidth)           return;
    if( J > h-BlockWidth)           return;
    if( d>1 && (K > d-BlockWidth))  return;
    if( d < 1 )                     return;
    mean1 = 0.0; mean2 = 0.0;
    max1  = 0.0; max2  = 0.0;

    // Compute the Mean matrix
    if( d>1 )
    {
      for(n = 0 ; n < BlockWidth3; ++n)
      {
        // Set the moving indexes in the block of the current neigh point
        tmp    = n % BlockWidth2;
        i      = I + tmp % BlockWidth;
        j      = J + tmp / BlockWidth;
        k      = K + n / BlockWidth2;
        mean1 += inputL[k * wh + j * w + i];
        mean2 += inputR[k * wh + j * w + i];
      }
    }
    else
    {
      for(n = 0 ; n < BlockWidth2; ++n)
      {
        // Set the moving indexes in the block of the current neigh point
        i      = I + n % BlockWidth;
        j      = J + n / BlockWidth;
        mean1 += inputL[j * w + i];
        mean2 += inputR[j * w + i];
      }
    }
    mean1 = mean1 / normalisation;
    mean2 = mean2 / normalisation;
    outputLMean[ig] = (float)mean1;
    outputRMean[ig] = (float)mean2;

    // Compute the Variance matrix
    if( d>1 )
    {
     for(k = K ; k <= K+BlockWidth; k++)
       for(j = J ; j <= J+BlockWidth; j++)
         for(i = I ; i <= I+BlockWidth; i++)
           {
             max1 += (inputL[k * wh + j * w + i] - mean1)*(inputL[k * wh + j * w + i] - mean1);
             max2 += (inputR[k * wh + j * w + i] - mean2)*(inputR[k * wh + j * w + i] - mean2);
           }
    }
    else
    {
      for(j = J ; j <= J+BlockWidth; j++)
        for(i = I ; i <= I+BlockWidth; i++)
        {
          max1 += (inputL[j * w + i] - mean1)*(inputL[j * w + i] - mean1);
          max2 += (inputR[j * w + i] - mean2)*(inputR[j * w + i] - mean2);
        }
    }
    max1 = max1 / normalisation;
    max2 = max2 / normalisation;
    outputLVariance[ig] = (float)(max1);
    outputRVariance[ig] = (float)(max2);
  }

  // Finish the threads
  barrier(CLK_LOCAL_MEM_FENCE);
}
