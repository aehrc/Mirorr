#ifndef __OPENCL_VERSION__ //defined by the kernel compiler
#define __kernel
#define __global
#define CLK_LOCAL_MEM_FENCE
#define CLK_GLOBAL_MEM_FENCE
#define uint unsigned int
#define size_t unsigned int
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
 * kMatch:
 * This OpenCL kernel ...
 *
 * @param imgL Base image buffer
 * @param imgR Search image buffer
 * @param meanL Input mean buffer for the inputL image
 * @param meanR Input mean buffer for the inputR image
 * @param varL Input variance buffer for the inputL image
 * @param varR Input variance buffer for the inputR image
 * @param kmatches Output buffer to store the best matches
 * @param scores Output buffer to store the score of the best matches
 * @param scores Output buffer to store the score of the 2nd best matches
 * @param w Width of the base image buffer
 * @param h Height of the base image buffer
 * @param d Depth of the base image buffer
 * @param BlockWidth size of the {2,3}D block
 * @param NhoodWidth size of the neighbourhood to search of correspondences
 * @param padding -- not used
 * @param BlockGap
 */
//__attribute__((vec_type_hint(float4)))
//#pragma OPENCL EXTENSION cl_amd_printf : enable
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_extended_atomics : enable
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
__kernel
void kMatch(
    __global float * imgL,          __global float * imgR,
    __global float * meanL,         __global float * meanR,
    __global float * varL,          __global float * varR,
    __global uint * kmatches, volatile __global uint * scores,__global uint * scores2nd,
    const uint w, const uint h, const uint d,
    uint BlockWidth,  uint NhoodWidth,  uint padding,
    uint BlockGap)
{
  const float EPS = 1e-10;

  // Gpu parameters
  int     ig          = get_global_id(0);
  size_t  gsize       = get_global_size(0);

  // Internal variables
  float   score       = 0.0f;
  int     i;
  int     j;
  int     k;
  int     r;
  int     s;
  int     t;
  int     r2;
  int     s2;
  int     t2;
  int     jp;
  int     n;
  int     tmp;
  float   v1;
  float   v2;
  float   imR1Val;
  float   imR2Val;
  float   bestscore   = 0.0f;
  uint    bestMatch   = 0;
  uint    old;
  uint    mask;
  uint  * iflip;

  // Algorithm parameters
  int     wh          = w * h;
  int     BlockWidth2 = BlockWidth*BlockWidth;
  int     NhoodWidth2 = NhoodWidth*NhoodWidth;
  int     BlockWidth3 = BlockWidth2*BlockWidth;
  int     NhoodWidth3 = NhoodWidth2*NhoodWidth;
  int     Qi          = w * h * d;
  int     idivs       = ig % Qi;                        // Current voxel position
  int     jdivs       = ig / Qi;                        // Current part neighbourhood division
  int     Qneigh      = NhoodWidth3;
  int     Divs        = (int)((float)gsize / (float)Qi);
          Divs        = (Divs < Qneigh) ? Divs : Qneigh;
  int     dQ          = Qneigh / Divs;
  int     limit       = Divs * Qi;

  // Debug variables
  int     dbgI        = 67;
  int     dbgJ        = 107;
  int     dbgK        = 67;
  int     dbgJdivs    = 0;
  int     arg         = 0;
  int     dbgR        = 9;
  int     dbgS        = 42;
  int     dbgT        = 14;
#if 0
  float   moyg;
  float   moyd;
  float   d1;
  float   d2;
  float   baseSearchProd  = 0.0f;
  float   baseSumSqr      = 0.0f;
  float   searchSumSqr    = 0.0f;
#else
  double   moyg;
  double   moyd;
  double   d1;
  double   d2;
  double   baseSearchProd  = 0.0f;
  double   baseSumSqr      = 0.0f;
  double   searchSumSqr    = 0.0f;
#endif

#if 0
  int foundInfinity = 0;
  int foundInfinity2 = 0;
  int foundInfinity3 = 0;
  int foundInfinity4 = 0;
#endif

  if(ig < limit)
  {
    // Set the global indexes
    r = idivs % wh;
    i = r % w;
    j = r / w;
    k = idivs / wh;

    // initialize the best match
    bestMatch = (uint) (k * wh + j * w + i); //DRH: bestMatch == idivs

    // Check the edges
    if (!(i > w - BlockWidth || j > h - BlockWidth
        || (d > 1 && (k > d - BlockWidth)))) //DRH: uint bug if {w,h,d} < BlockWidth
    {
      // Processing on 2D images
      if (d == 1)
      {
        // For each point in the neighbourhood that the kernel
        // has to threat
        baseSearchProd = 0.0f;
        baseSumSqr = 0.0f;
        searchSumSqr = 0.0f;
        for (jp = jdivs * dQ; jp < (jdivs + 1) * dQ; ++jp)
        {
          // Set the moving indexes in the neigh
          tmp = jp % NhoodWidth2;
          r = i + BlockGap * (tmp % NhoodWidth - NhoodWidth / 2);
          s = j + BlockGap * (tmp / NhoodWidth - NhoodWidth / 2);
          // Compute the moving values
          // get the mean
          score = 0.0f;
          moyd = meanR[s * w + r];
          moyg = meanL[j * w + i];
          // Work on the current block centered on r,s,t
          for (n = 0; n < BlockWidth2; ++n)
          {
            // Set the moving indexes in the block of the current neigh point
            tmp = n % BlockWidth2;
            r2 = i + tmp % BlockWidth;
            s2 = j + tmp / BlockWidth;
            d1 = imgL[s2 * w + r2] - moyg;

            r2 = r + tmp % BlockWidth;
            s2 = s + tmp / BlockWidth;
            d2 = imgR[s2 * w + r2] - moyd;
            baseSearchProd += d1 * d2;
            baseSumSqr += d1 * d1;
            searchSumSqr += d2 * d2;
          }
          // compute the normalized score (DRH: squared Normalized Cross-Correlation)
          score =
              (baseSumSqr <= EPS || searchSumSqr <= EPS) ?
                  -2.f :
                  baseSearchProd * baseSearchProd / baseSumSqr / searchSumSqr;
          // do we keep the match ?
          if (score > bestscore)
          {
            bestscore = score;
            bestMatch = (uint) (s * w + r);
          }
        }
      }
      // Processing on 3D images
      else
      {
        // Compute the fixed values
        // get the mean
        //moyg  = meanL[k * wh + j * w + i];
        //d1    = fabs(imgL[k * wh + j * w + i] - moyg);
        // get the standard deviation
        //v1    = sqrt(varL[k * wh + j * w + i]);
        // Check if the fixed point is good enough
        // we can check with a median variance here
        //if ( v1 != 0.0f )
        //{
        // For each point in the neighbourhood that the kernel
        // has to threat
        for (jp = jdivs * dQ; jp < (jdivs + 1) * dQ; ++jp)
        {
          // Set the moving indexes in the neigh
          tmp = jp % NhoodWidth2;
          r = i + ((int)BlockGap) * ((int)(tmp % NhoodWidth) - (int)(NhoodWidth / 2));
          s = j + ((int)BlockGap) * ((int)(tmp / NhoodWidth) - (int)(NhoodWidth / 2));
          t = k + ((int)BlockGap) * ((int)(jp / NhoodWidth2) - (int)(NhoodWidth / 2));

          if ( (r<0) || (r > w - BlockWidth) ||
               (s<0) || (s > h - BlockWidth) ||
               (t<0) || (t > d - BlockWidth) ) {
            continue; //DRH - block outside the image domain
          }

          // Compute the moving values
          // get the mean
          score = 0.0f;
          moyd = meanR[t * wh + s * w + r];
          moyg = meanL[k * wh + j * w + i];

          /* // - - - - - - - - - -
          // DRH - compute double precision means - to remove:
          // - - - - - - - - - -
          double d_moyd = 0., d_moyg = 0.;
          for (n = 0; n < BlockWidth3; ++n)
          {
            tmp = n % BlockWidth2;

            r2 = i + tmp % BlockWidth;
            s2 = j + tmp / BlockWidth;
            t2 = k + n / BlockWidth2;
            d_moyg += imgL[t2 * wh + s2 * w + r2];

            r2 = r + tmp % BlockWidth;
            s2 = s + tmp / BlockWidth;
            t2 = t + n / BlockWidth2;
            d_moyd += imgR[t2 * wh + s2 * w + r2];
          }
          moyd = d_moyd / (double)BlockWidth3;
          moyg = d_moyg / (double)BlockWidth3;
          // - - - - - - - - - - */

          // Work on the current block centered on r,s,t
          baseSearchProd = 0.0f;
          baseSumSqr = 0.0f;
          searchSumSqr = 0.0f;
          for (n = 0; n < BlockWidth3; ++n)
          {
            tmp = n % BlockWidth2;

            r2 = i + tmp % BlockWidth;
            s2 = j + tmp / BlockWidth;
            t2 = k + n / BlockWidth2;
            d1 = imgL[t2 * wh + s2 * w + r2] - moyg;

            r2 = r + tmp % BlockWidth;
            s2 = s + tmp / BlockWidth;
            t2 = t + n / BlockWidth2;
            d2 = imgR[t2 * wh + s2 * w + r2] - moyd;

            baseSearchProd += d1 * d2;
            baseSumSqr += d1 * d1;
            searchSumSqr += d2 * d2;
          }
          // compute the normalized score (DRH: squared Normalized Cross-Correlation)
          score =
              (baseSumSqr <= EPS || searchSumSqr <= EPS) ?
                  -2.f :
                  (float) (baseSearchProd * baseSearchProd / baseSumSqr / searchSumSqr);

          //score = baseSumSqr;
          //score = moyd; ////DRH DEBUG - GPU-CPU equal
          //score = searchSumSqr;
          //score = baseSearchProd;

          //DRH: the above code will sometimes still result in score=Inf due
          //     to limited float 32 precision, we take care of those here
          if (score > 3.3e33) {
            score = -2.f;
          }

          // do we keep the match ?
          if (score > bestscore)
          {
            bestscore = score;
            bestMatch = (uint) (t * wh + s * w + r);
          }
        }
        //}
      }
    }


    // A very good trick to compute the max value of the score
    // on all the kernels directly (with float values) !
    // DRH: Create a 32bits uint representation of any 32bits floating point number.
    //      This representation preserves ordering.
    iflip = (uint*) &bestscore;
    mask = -(int) ((*iflip) >> 31) | 0x80000000;
    mask = (*iflip) ^ mask;

#if 0
    // DRH: Check for infinity (due to limited float precision)
    if (mask == 4286578688U) {
      mask = 3212836864U; // == 1.0f
      //mask = 2147483648; // == 0.f
    }
#endif

    //old = atom_max(&scores[idivs], mask);
    //barrier(CLK_LOCAL_MEM_FENCE);

    old = atom_max(&scores[idivs], mask);
    barrier(CLK_GLOBAL_MEM_FENCE);
    uint drh = scores[idivs];

    // If this kernel has the best score, the score is saved and
    //if (old < mask)
    if ( (drh == mask) && (old != mask) )
    {
      // we save the match (there should only have one update following the barrier)
      kmatches[k * wh + j * w + i] = (uint) bestMatch;
      //kmatches[k * wh + j * w + i] = (uint) (1001 * wh + 2 * w + 1);
      //kmatches[k * wh + j * w + i] = (uint) (1001 * wh + s * w + r);
      //kmatches[k * wh + j * w + i] = (uint) (t * wh + s * w + r);
      //scores2nd[k * wh + j * w + i] = (uint) old;
      atom_max(&scores2nd[idivs], old);
    }
    else //if we don't have the best, maybe we have the 2nd best
    {
      atom_max(&scores2nd[idivs], mask);
    }

    /*DRH - useless?
    if ((int) scores2nd[k * wh + j * w + i] == 0)
    {
      scores2nd[k * wh + j * w + i] = (uint) 1069547519; //Flag value == -3.f when converted as above
    }*/
  }
}
