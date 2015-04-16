#ifndef __OPENCL_VERSION__ //defined by the kernel compiler
#define __kernel
#define __global
#define CLK_LOCAL_MEM_FENCE
#define CLK_GLOBAL_MEM_FENCE
#define uint unsigned int
#define size_t unsigned int
#endif

// define NO_COMPUTE if you want to measure only the memory access

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

/**
 *
 * --- WARNING ---
 * The .h version of this file is generated from the .cl file. You
 * MUST run convertCl2Header.sh after any modification to the .cl file!
 * --- WARNING ----
 *
 *
 * Code by Maciej Golebiewski, 13 Jun 2014
 * 
 * Copyright (c) 2009-15 CSIRO. All rights reserved.
 * 
 * For complete copyright, license and disclaimer of warranty
 * information see the LICENSE.txt file for details.
 * 
 * This software is distributed WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the above copyright notice for more information.
 * 
 */


/* Parallel sum across a work group, using a share vector in local memory
 * This kernel is only called from newKmatch, not directly from host.
 *
 * @param n number of work items (threads) in work group
 * @param t serial/1D local index of this work item in its work group
 * @param rblock vector of n elements in local memory, shared across work group
 *
 * The result is in rblock[0]; if work items other than the first one (with t==0)
 * want to access the result they need to call memory fence first after parsum returns.
 */
__kernel
inline void parsum(const int n, const int t, double iv, __local double *rblock) {

	double s = iv;
	int i;

	for (i = n/2; i > 0; i = i/2) {
		if (t >= i) rblock[t] = s;
		barrier (CLK_LOCAL_MEM_FENCE);
		if (t < i) s += rblock[t+i];
	}
	if (t == 0) rblock[t] = s;

} // parsum


// voxel offset; must allow for negative values
// to represent voxel outside (to the left of) image
typedef signed int voff_t;


__kernel
void optKmatch (
	__read_only image3d_t imgL,
	__read_only image3d_t imgR,
    __global uint * kmatches, volatile __global float* scores,__global float* scores2nd,
    __local double *rblock,
    __local double *d1vec,
    const uint w, const uint h, const uint d,
    const uint BlockWidth,  const uint NhoodWidth,  const uint padding,
    const uint BlockGap, const uint NhoodGap,
    const uint padx, const uint pady, const uint padz)
{
	const float EPS = 1e-10;
	const sampler_t sampler = CLK_FILTER_NEAREST | CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE;

	int4   pxl;	    // temp coords variable for calling read_imagef
	float4 ImA;         // value of my pixel in the starting image
	float4 ImB;         // value of my pixel in the matching block

	size_t wgdim[3];    // dimensions of the work group
	size_t myblock0[3]; // 3D pixel coords of first pixel in starting block

	size_t mypixel[3];  // which pixel in a block this work item is responsible for?

	// to allow for negative coords
	voff_t nhood0[3];   // 3D pixel coords of first pixel in first matching block
	voff_t nhood1[3];   // 3D pixel coords of first pixel in last matching block

	// to allow for negative coords
	voff_t nblock0[3];  // 3D pixel coords of first pixel in matching block

        double   moyg;  // mean for starting image block
        double   moyd;  // mean for matching image block
	double   d1, d1sqr;
	double   d2, d2sqr;
	double   d1d2prod;
	double   baseSearchProd;
	double   baseSumSqr;
	double   searchSumSqr;

	float	score;

	float   bestscore = -3.0f;
	float   secondBest = -3.0f;

	size_t slid, snum, mybase;

	voff_t nhood2; // half pixel length of neighbourhood

	uint    bestMatch = 0;

	uint	block_index, ln_off;

	int	i;

	// length of neighbourhood edge in pixels
	nhood2 = (voff_t)((NhoodWidth-1)*BlockGap + BlockWidth)/2;

	for (i = 0; i < 3; i++) {

		// get coords of my block in starting image
		myblock0[i] = get_group_id(i);

		// get my pixel offsets in a block
		mypixel[i] = get_local_id(i);

		// dimensions of my work group
		wgdim[i] = get_local_size(i);
	}

	// linear index of my block for storing result in the output vectors
	block_index = myblock0[2]*get_num_groups(1)*get_num_groups(0) + myblock0[1]*get_num_groups(0) + myblock0[0];

	for (i = 0; i < 3; i++) {
		// compute coords of first pixel in my block
		// convert block coords to coords of first pixel in block
		myblock0[i] = myblock0[i] * NhoodGap;
	}

	myblock0[0] += padx;
	myblock0[1] += pady;
	myblock0[2] += padz;

	for (i = 0; i < 3; i++) {
		// compute coords of first pixel in matching neighborhood
		// (may fall outside image)
		nhood0[i] = (voff_t)(myblock0[i] + (BlockWidth/2))-nhood2;
	}

	// serial local index in working group
	mybase = mypixel[1]*wgdim[0] + mypixel[2]* wgdim[0]*wgdim[1];
	slid = mypixel[0] + mybase;
	// total number of pixels in a block (or work items in a group)
	snum = wgdim[0]*wgdim[1]*wgdim[2];

	if (slid == 0) {
		scores[block_index] = -1.0;
		scores2nd[block_index] = -2.0;
		kmatches[block_index] = (uint)0;
	}

	// read and remember the value of my pixel in the starting image
	pxl.x = myblock0[0] + mypixel[0];
	pxl.y = myblock0[1] + mypixel[1];
	pxl.z = myblock0[2] + mypixel[2];
	ImA = read_imagef(imgL, sampler, pxl);

	// parallel sum across work group to compute
	// meang of starting image block
	parsum(snum, slid, ImA.s0, rblock);
	barrier (CLK_LOCAL_MEM_FENCE);
	// mean for my starting block image
	moyg = rblock[0] / snum;

        d1 = ImA.s0 - moyg;

	// save d1 in shared memory for re-use by other work items
	// when computing d1d2prod
	d1vec[slid] = d1;

	// compute and remember base sum square
	d1sqr = d1 * d1;

	// parallel sum across work group to compute
	// baseSumSqr
	parsum(snum, slid, d1sqr, rblock);
	barrier (CLK_LOCAL_MEM_FENCE);
	baseSumSqr = rblock[0];

	// compute neighbourhood boundaries for the 3 dimensions
	for (i = 0; i < 3; i++) {
		// end boundary
		nhood1[i] = nhood0[i] + BlockGap * (NhoodWidth - 1);
		// adjust start boundary to inside the image
		for (; nhood0[i] < 0; nhood0[i] += BlockGap) ;
	}

	// ensure end boundaries are inside the image
	nhood1[0] = (nhood1[0] < w-BlockWidth) ? nhood1[0] : w-BlockWidth;
	nhood1[1] = (nhood1[1] < h-BlockWidth) ? nhood1[1] : h-BlockWidth;
	nhood1[2] = (nhood1[2] < d-BlockWidth) ? nhood1[2] : d-BlockWidth;

	// z coord of the first layer/slice
	pxl.z = nhood0[2] + mypixel[2];

	// loop over all layers of the matching neighbourhood
	for (nblock0[2] = nhood0[2]; nblock0[2] <= nhood1[2]; nblock0[2] += BlockGap) {

		// y coord of the first row in the current layer
		pxl.y = nhood0[1] + mypixel[1];

		// loop over all rows of the matching neighbourhood
		for (nblock0[1] = nhood0[1]; nblock0[1] <= nhood1[1]; nblock0[1] += BlockGap) {

			// x coord of the first column in the current row
			pxl.x = nhood0[0] + mypixel[0];

			// pixels in the first block read by all work items
			ImB = read_imagef(imgR, sampler, pxl);

			// loop over all columns of the matching neighbourhood
			// (same as moving matching block along the row)
			for (nblock0[0] = nhood0[0]; nblock0[0] <= nhood1[0]; nblock0[0] += BlockGap) {

				// is my pixel outside the current matching block?
				if (pxl.x < nblock0[0]) {
					// yes, calculate new pixel coord
					pxl.x += BlockWidth;
					// and read new value
					ImB = read_imagef(imgR, sampler, pxl);
				}

				// linear offset of my current pixel relative to start of the current
				// matching block
				// (used to find d1 corresponding to d2 for my current pixel)
				ln_off = pxl.x - nblock0[0] + mybase;
				d1d2prod = d1vec[ln_off];

				// parallel sum across work group to compute
				// meang of matching image block
				parsum(snum, slid, ImB.s0, rblock);
				barrier (CLK_LOCAL_MEM_FENCE);
				moyd = rblock[0] / snum;

				d2 = ImB.s0 - moyd;
				d2sqr = d2 * d2;
				d1d2prod = d1d2prod * d2;

				// parallel sums across work group to compute
				// searchSumSqr
				parsum(snum, slid, d2sqr, rblock);
				barrier (CLK_LOCAL_MEM_FENCE);
				searchSumSqr = rblock[0];

				// parallel sums across work group to compute
				// baseSearchProd
				parsum(snum, slid, d1d2prod, rblock);
				barrier (CLK_LOCAL_MEM_FENCE);
				baseSearchProd = rblock[0];

				// to avoid branches all work items in work group
				// compute score but only the first work item in
				// each group will store the result after the loops
				score = (baseSumSqr <= EPS || searchSumSqr <= EPS) ?
					-2.f :
				(float) (baseSearchProd * baseSearchProd / baseSumSqr / searchSumSqr);
				score = score < 3.3e33 ? score : -2.f;

				// keep per workgroup best score and best match
				if (score > bestscore) {
					secondBest = bestscore;
					bestscore = score;
					// host code expects index of first voxel in block;
					// not block index
					// index/offset of the first voxel in the block
					bestMatch = (uint)(nblock0[2]*w*h + nblock0[1]*w + nblock0[0]);
				}
				else {
					secondBest = (secondBest > score) ? secondBest : score;
				}

			} // x dimension

			// advance to next row
			pxl.y += BlockGap;

		} // y dimension

		// advance to next layer
		pxl.z += BlockGap;

	} // z dimension

	// store the best score in result vectors
	if (slid == 0) {
		scores[block_index] = bestscore;
		scores2nd[block_index] = secondBest;
		kmatches[block_index] = bestMatch;
	}


} // optKmatch
