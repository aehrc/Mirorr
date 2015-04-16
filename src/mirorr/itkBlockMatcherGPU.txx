/*=========================================================================
Program: mirorr
Module: itkBlockMatcherGPU2.h
Author: Maciej Golebiewski
Created: 13 Jun 2014
Based on:
  Module: itkBlockMatcherGPU.h by 
  Author: Jeremy Coatelen
  Created: 11 Apr 2011

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
// Revised: 30 Jan 2013 by D Rivest-Henault. Moved the kernel source code to external header files.
// Revised: 25 Jan 2013 by D Rivest-Henault. Patched Initialize() to make sure that the
//                                           kernel "kMeanVariance" is always the first kernel.
// Revised: Jun - Sep 2014 by M Golebiewski. Kernel creation and lauching (setting kernels arguments)
//					     modified for the new OpenCL kernel that combines mean
//					     variance and block matching computations in a single
//					     kernel.

#ifndef __itkBlockMatcherGPU_txx
#define __itkBlockMatcherGPU_txx

#include "itkBlockMatcherGPU.h"

#include <sstream>
#include <algorithm>


// Utility function and MACROS to help manage exceptions
const std::string makeExceptionString(std::string msg, std::string sourceFile, int line) {
  std::stringstream ss;
  ss << msg << "  (" << sourceFile << "-" << line << ")";
  return ss.str();
}

// ----------------------------------------------------------------------------------------
// Legacy macros that should be considered as deprecated, we need something better in 2015
#ifdef __linux__
#define BASENAME_TAG  basename(__FILE__)
#else
  #define BASENAME_TAG ""
#endif
#define RESET_COLOR ""
#define DEBUG_TAG "[DEBUG]"
#define INFO_TAG RESET_COLOR << "[INFO]"
#define ERROR_TAG "[ERROR]"
#define WARNING_TAG "[WARNING]"

#define BM_GPU_EXCEPTION(x) \
                BlockMatcherGPUException(makeExceptionString(x, BASENAME_TAG, __LINE__))
#define BM_GPU_EXCEPTION_ALGO(x) \
                BlockMatcherGPUExceptionAlgo(makeExceptionString(x, BASENAME_TAG, __LINE__))
#define BM_GPU_EXCEPTION_PLATFORM(x) \
                BlockMatcherGPUExceptionPlatform(makeExceptionString(x, BASENAME_TAG, __LINE__))

#define INFO_INLINE(x) \
                std::cout << x
#define INFO_MARK(x) \
                std::cout << INFO_TAG << "\t\t" << x << std::endl
#define INFO_GATHER(x) \
                std::cout << INFO_TAG << "\t\t"; \
                x \
                std::cout << std::endl
#define INFO_IFPRINT(ETI,x) if( ETI & x ) INFO_MARK(#ETI)
#define INFO_ELIFPRINT(ETI,x) else if( ETI & x ) INFO_MARK(#ETI)
#define INFO_ELSEPRINT(ETI) else INFO_MARK(#ETI)

#define ERROR_MARK(x) \
                std::cerr << ERROR_TAG << "\t\t" << x << "\t(" << BASENAME_TAG << "-" <<__LINE__ << ")" << RESET_COLOR << std::endl
#define ERROR_IFPRINT(ETI,x) if( ETI & x ) ERROR_MARK(#ETI)
#define ERROR_IFEQUALPRINT(ETI,x) if( ETI == x ) ERROR_MARK(#ETI)
// ----------------------------------------------------------------------------------------


// need uint32_t for compatibility with OpenCL uint
#if (defined(__GLIBCXX__) || (__cpluspluss >= 199901L))
#include <stdint.h>
#elif defined(_MSC_VER)
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
#endif

/**
 * Loading the source code for our two kernels at compile time.
 */
#ifdef GPU_BM_USE_MASK
  #include "kernels/optKmatchMask.h"
#else
  #include "kernels/optKmatch.h"
#endif

// symbolic constant to make code referring to this kernel more readable
const int kMatchID = 0;

namespace itk
{

  template < typename ImageType >
  BlockMatcherGPU<ImageType>::BlockMatcherGPU() :
  Superclass()
  {
    this->region[0] = 0;
    this->region[1] = 0;
    this->region[2] = 0;

    this->blockGrid[0] = 0;
    this->blockGrid[1] = 0;
    this->blockGrid[2] = 0;

    this->input_L_ = 0;
    this->input_R_ = 0;
#ifdef GPU_BM_USE_MASK
    this->input_ML_ = 0;
    this->input_MR_ = 0;
#endif
    this->m_Verbosity = 0;
    this->releasedAll_ = true;
    this->releasedBuffers_ = true;
    this->recycleImageBufferData_ = false;
    this->Initialize();
  }

  template < typename ImageType >
  BlockMatcherGPU<ImageType>::~BlockMatcherGPU()
  {
    this->release();
  }

  template < typename ImageType >
  void BlockMatcherGPU<ImageType>::release()
  {
    releaseBuffers();
    if(  ! this->releasedAll_ )
      {
        clReleaseKernel(kernels[kMatchID]);
        clReleaseProgram(program);
        clReleaseCommandQueue(commands);
        clReleaseContext(context);
        releasedAll_ = true;
      }
  }

  template < typename ImageType >
  void BlockMatcherGPU<ImageType>::releaseBuffers()
  {
    if(  ! this->releasedBuffers_ )
      {
        // image3D objects released using standard API?
        clReleaseMemObject(input_L_);
        clReleaseMemObject(input_R_);
        this->input_L_ = 0;
        this->input_R_ = 0;

#ifdef GPU_BM_USE_MASK
        clReleaseMemObject(input_ML_);
        clReleaseMemObject(input_MR_);
        this->input_ML_ = 0;
        this->input_MR_ = 0;
#endif
        clReleaseMemObject(output_Matches);
        clReleaseMemObject(output_Scores);
        clReleaseMemObject(output_Scores2nd);
        this->output_Matches = 0;
        this->output_Scores = 0;
        this->output_Scores2nd = 0;

        this->region[0] = 0;
        this->region[1] = 0;
        this->region[2] = 0;

        this->blockGrid[0] = 0;
        this->blockGrid[1] = 0;
        this->blockGrid[2] = 0;
        releasedBuffers_ = true;
      }
  }

  template < typename ImageType >
  void BlockMatcherGPU<ImageType>::SwapImageOrder()
  {
    AbstractBlockMatcher<ImageType>::SwapImageOrder();

    cl_mem tmp = this->input_L_;
    this->input_L_ = this->input_R_;
    this->input_R_ = tmp;

#ifdef GPU_BM_USE_MASK
    tmp = this->input_ML_;
    this->input_ML_ = this->input_MR_;
    this->input_MR_ = tmp;
#endif

    this->recycleImageBufferData_ = true;
  }

  template < typename ImageType >
  void BlockMatcherGPU<ImageType>::GetOffsets(
      PointListType &blockPositionsOut,       // outputs
      PointListType &blockMatchPositionsOut   //
  )
  {
    if( this->m_Verbosity >= 3 )
      INFO_MARK("..:: Match ::..");
    const int ImDim = ImageType::ImageDimension;

    if (! this->recycleImageBufferData_) {
      this->updateBuffers();
    } else {
      if( this->m_Verbosity >= 3 )
      INFO_MARK("Not updating buffers.");
    }
    this->recycleImageBufferData_ = false;

    ////////////////////////////////////////////////////////////////////////////// Kernels Arguments

    // image dimensions
    RegionType  baseRegion      = this->m_BaseImage->GetLargestPossibleRegion();
    SizeType sizs = baseRegion.GetSize();
    const unsigned int w = sizs[0];
    const unsigned int h = sizs[1];
    const unsigned int d = sizs[2];

    // number of blocks in starting image, along each dimension
    // will be used to calculate total number of global work items
    size_t nblocks[3];

    nblocks[0] = (w - this->m_BlockWidth) / this->m_NhoodGap + 1;
    nblocks[1] = (h - this->m_BlockWidth) / this->m_NhoodGap + 1;
    nblocks[2] = (d - this->m_BlockWidth) / this->m_NhoodGap + 1;

    unsigned int gapAtTheEnd[3];
    unsigned int startOffset[3];

    // compute start offset so at to centre the matched region in the image;
    // to avoid bias away from centre
    for (int i = 0; i < ImDim; i++) {
	    int regSize = this->m_NhoodGap * (nblocks[i]-1) + this->m_BlockWidth;
	    gapAtTheEnd[i] = sizs[i] - regSize;
	    startOffset[i] = gapAtTheEnd[i] / 2;
    }

    // number of voxels or work items in a block
    const size_t blocksize = this->m_BlockWidth * this->m_BlockWidth * this->m_BlockWidth;

    // number of voxels in an image
    //count = sizs[0] * sizs[1] * sizs[2];

    if( this->m_Verbosity >= 3 ) {
      INFO_MARK("Parameters "
          "[" << w << "x" << h << "x" << d << "]"
          << "  - Nhood:[gap="<< this->m_NhoodGap << "px, width=" << this->m_NhoodWidth
          << "blks] - Block:[gap="<< this->m_BlockGap << "px, width=" << this->m_BlockWidth
          << "px] - Padding "<< this->m_Padding );
    }

    int arg = 0;
    this->err = 0;
    this->err  = clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(cl_mem), &(this->input_L_));
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(cl_mem), &(this->input_R_));
#ifdef GPU_BM_USE_MASK
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(cl_mem), &(this->input_ML_));
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(cl_mem), &(this->input_MR_));
#endif
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(cl_mem), &(this->output_Matches));
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(cl_mem), &(this->output_Scores));
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(cl_mem), &(this->output_Scores2nd));
    // shared mem, must match the size of a work group which equal size of matching block
    // for parallel sum reduction
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, (blocksize)*sizeof(double), NULL); // shared mem
    // for pivoting d1 buffer
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, (blocksize)*sizeof(double), NULL); // shared mem
#ifdef GPU_BM_USE_MASK
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, (blocksize)*sizeof(unsigned int), NULL); // shared mem
#endif
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(unsigned int), &w);
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(unsigned int), &h);
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(unsigned int), &d);
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(unsigned int), &(this->m_BlockWidth));
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(unsigned int), &(this->m_NhoodWidth));
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(unsigned int), &(this->m_Padding));
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(unsigned int), &(this->m_BlockGap));
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(unsigned int), &(this->m_NhoodGap));
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(unsigned int), &(startOffset[0]));
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(unsigned int), &(startOffset[1]));
    this->err |= clSetKernelArg(this->kernels[kMatchID], arg++, sizeof(unsigned int), &(startOffset[2]));

    if (this->err != CL_SUCCESS) {
        ERROR_MARK( "Failed to set kernel arguments! " << this->err );

        ERROR_IFPRINT(CL_INVALID_KERNEL,this->err);
        ERROR_IFPRINT(CL_INVALID_ARG_INDEX,this->err);
        ERROR_IFPRINT(CL_INVALID_ARG_VALUE,this->err);
        ERROR_IFPRINT(CL_INVALID_MEM_OBJECT,this->err);
        ERROR_IFPRINT(CL_INVALID_SAMPLER,this->err);
        ERROR_IFPRINT(CL_INVALID_ARG_SIZE,this->err);

        throw BM_GPU_EXCEPTION_PLATFORM("Failed to set kernel arguments! ");
    }

    // total number of blocks in starting image
    size_t block_count = nblocks[0] * nblocks[1] * nblocks[2];
    
    if( this->m_Verbosity >= 2 ) {
      INFO_MARK("no of blocks in starting image: " << nblocks[0] << "x" << nblocks[1] << "x" << nblocks[2] << "=" << block_count);
      INFO_MARK("no of pixels in image: " << w << "x" << h << "x" << d << "= " << w*h*d);
    }
    

#ifdef DEBUG_GPU_CPU_CONSISTENCY
    this->scoresOut.clear();
    this->scoresOut.reserve(block_count);
    this->scoresOut2nd.clear();
    this->scoresOut2nd.reserve(block_count);
    std::cout << "USING itkBlockMatcherGPU2.txx" << std::endl;
#endif

    ////////////////////////////////////////////////////////////////////////////// Matches Computation

    // global numbers of work items along each dimension in the starting image
    size_t gworks[3];

    for (int i = 0; i < 3; i++) gworks[i] = nblocks[i] * this->m_BlockWidth;

    // number of work items along each dimension of work group
    size_t mblocks[3];
    mblocks[0] = mblocks[1] = mblocks[2] = this->m_BlockWidth;

    cl_event profile;

    this->err = clEnqueueNDRangeKernel(this->commands, this->kernels[kMatchID],
        3, NULL, gworks, mblocks,
        0, NULL, &profile); //Run the "kMatch" kernel
    if (this->err) {
        ERROR_MARK( "Failed to execute kernel!" );

        ERROR_IFEQUALPRINT(CL_INVALID_PROGRAM_EXECUTABLE,err);
        ERROR_IFEQUALPRINT(CL_INVALID_COMMAND_QUEUE,err);
        ERROR_IFEQUALPRINT(CL_INVALID_KERNEL,err);
        ERROR_IFEQUALPRINT(CL_INVALID_CONTEXT,err);
        ERROR_IFEQUALPRINT(CL_INVALID_KERNEL_ARGS,err);
        ERROR_IFEQUALPRINT(CL_INVALID_WORK_DIMENSION,err);
        ERROR_IFEQUALPRINT(CL_INVALID_WORK_GROUP_SIZE,err);
        ERROR_IFEQUALPRINT(CL_INVALID_WORK_ITEM_SIZE,err);
        ERROR_IFEQUALPRINT(CL_INVALID_GLOBAL_OFFSET,err);
        ERROR_IFEQUALPRINT(CL_OUT_OF_RESOURCES,err);
        ERROR_IFEQUALPRINT(CL_MEM_OBJECT_ALLOCATION_FAILURE,err);
        ERROR_IFEQUALPRINT(CL_INVALID_EVENT_WAIT_LIST,err);
        ERROR_IFEQUALPRINT(CL_OUT_OF_HOST_MEMORY,err);

        throw BM_GPU_EXCEPTION_PLATFORM("Failed to execute kernel!");
    }

    // allocate host buffer to copy results from device
    unsigned int *tmpMatches = new unsigned int[block_count];
    float *tmpScores = new float[block_count];
    float *tmpScores2nd = new float[block_count];

    // reset output vectors
    blockPositionsOut.clear();
    blockMatchPositionsOut.clear();

    // Wait for all commands to complete
    clFinish(this->commands);
    if (this->err != CL_SUCCESS) {
        ERROR_MARK( "Error waiting for kMatch kernel! " <<  this->err );
        ERROR_IFEQUALPRINT(CL_INVALID_COMMAND_QUEUE,this->err);
        ERROR_IFEQUALPRINT(CL_OUT_OF_HOST_MEMORY,this->err);
        ERROR_IFEQUALPRINT(CL_OUT_OF_RESOURCES,this->err);
        throw BM_GPU_EXCEPTION_PLATFORM("Error waiting for kMatch kernel! ");
    }

    // measure kernel's execution time
    if( this->m_Verbosity > 2 ) {
	    cl_ulong tstart, tend;
	    float elapsed;

	    clGetEventProfilingInfo (profile, CL_PROFILING_COMMAND_START,
         sizeof(cl_ulong), &tstart, NULL);
	    clGetEventProfilingInfo (profile, CL_PROFILING_COMMAND_END,
         sizeof(cl_ulong), &tend, NULL);
	    elapsed = (tend - tstart) * 1e-9f; // convert from nanosec to seconds
	    INFO_MARK ("optKmatch kernel walltime: " << elapsed << " seconds");
    }


    this->err  = clEnqueueReadBuffer( this->commands, this->output_Matches,
        CL_TRUE, 0, sizeof(unsigned int) * block_count,
        tmpMatches, 0, NULL, NULL );
    if (this->err != CL_SUCCESS) {
        ERROR_MARK( "Failed to read output array (output_Matches)! " <<  this->err );
        ERROR_IFEQUALPRINT(CL_OUT_OF_RESOURCES,this->err);
        ERROR_IFEQUALPRINT(CL_INVALID_COMMAND_QUEUE,this->err);
        ERROR_IFEQUALPRINT(CL_INVALID_CONTEXT,this->err);
        ERROR_IFEQUALPRINT(CL_INVALID_MEM_OBJECT,this->err);
        ERROR_IFEQUALPRINT(CL_INVALID_VALUE,this->err);
        ERROR_IFEQUALPRINT(CL_INVALID_EVENT_WAIT_LIST,this->err);
        ERROR_IFEQUALPRINT(CL_MEM_OBJECT_ALLOCATION_FAILURE,this->err);
        ERROR_IFEQUALPRINT(CL_OUT_OF_HOST_MEMORY,this->err);
        throw BM_GPU_EXCEPTION_PLATFORM("Failed to read output array (output_Matches)! ");
    }

    this->err  = clEnqueueReadBuffer( this->commands, this->output_Scores,
        CL_TRUE, 0, sizeof(float) * block_count,
        tmpScores, 0, NULL, NULL );
    if (err != CL_SUCCESS) {
        ERROR_MARK( "Failed to read output array (output_Scores)! " <<  err );
        throw BM_GPU_EXCEPTION_PLATFORM("Failed to read output array (output_Scores)! ");
    }
    this->err  = clEnqueueReadBuffer( this->commands, this->output_Scores2nd,
        CL_TRUE, 0, sizeof(float) * block_count,
        tmpScores2nd, 0, NULL, NULL );
    if (err != CL_SUCCESS) {
        ERROR_MARK( "Failed to read output array (output_Scores2nd)! " <<  err );
        throw BM_GPU_EXCEPTION_PLATFORM("Failed to read output array (output_Scores2nd)! ");
    }

    IndexType matchPoint;
    IndexType sourcePoint;
    unsigned int u;

    //Convert GPU points which are in voxel units to ITK points which are in mm units
    {
      int nLowScores = 0, nMultimatch = 0;
      double score_avg = 0.f, score_min = 1e200, score_max = -1;

      int matchblock = 0;
      int foundblocks = 0;
      for (unsigned int z = 0; z < nblocks[2]; z++)
      {
        for (unsigned int y = 0; y < nblocks[1]; y++)
        {
          for (unsigned int x = 0; x < nblocks[0]; x++)
          {
            double currentBestScore = tmpScores[matchblock];
            double current2ndBestScore = tmpScores2nd[matchblock];

            if (currentBestScore > 1e-7)
            {
              if (fabs(currentBestScore - current2ndBestScore) > 1e-4)
              {

                sourcePoint[0] = startOffset[0] + x * this->m_NhoodGap;
                sourcePoint[1] = startOffset[1] + y * this->m_NhoodGap;
                sourcePoint[2] = startOffset[2] + z * this->m_NhoodGap;

                u = tmpMatches[matchblock] % (w * h);
                matchPoint[0] = u % w;
                matchPoint[1] = u / w;
                matchPoint[2] = tmpMatches[matchblock] / (w * h);

#ifdef DEBUG_GPU_CPU_CONSISTENCY
                this->scoresOut.push_back(currentBestScore);
                this->scoresOut2nd.push_back(current2ndBestScore);
#endif
                blockPositionsOut.push_back(this->blockPositionToPoint(sourcePoint));
                blockMatchPositionsOut.push_back(this->blockPositionToPoint(matchPoint));

                score_avg += currentBestScore;
                score_min = std::min(currentBestScore, score_min);
                score_max = std::max(currentBestScore, score_max);

                foundblocks++;
              }
              else
              {
                nMultimatch++;
              }
            }
            else
              nLowScores++;

            matchblock++;
          }
        }
      }

      if( this->m_Verbosity >= 2 ) {
        INFO_MARK("blocks inserted: " << foundblocks);
        INFO_MARK("block_count="<<block_count<<"; matchblock="<<matchblock);
      }
      if (this->m_Verbosity >= 2)
      {
        // TODO: low variance i think is wrong; compare with host implementation
        std::cout << "Blocks: Initial: " << block_count << ", Cull "
            << 100 * (1.0 - this->m_PortionMatchesKept) << "% low variance (<="
            << 0 << "): " << block_count
            << ", Cull " << nLowScores << " dissimilar: "
            << block_count - nLowScores
            << ", Cull " << nMultimatch << " multi-match: "
            << blockPositionsOut.size();

        if (!blockPositionsOut.empty())
        {
          std::cout << "  Scores: " << score_min << " to " << score_max << " ("
              << score_avg / blockPositionsOut.size() << ")" << std::endl;
        }
        else
          std::cout << "  Scores: " << "None" << std::endl;
      }
      if (blockPositionsOut.size() == 0)
      {
        throw BM_GPU_EXCEPTION_ALGO( \
            "WARNING: No block matches occurred. Images probably have" \
            " zero overlap. Try supplying a better initial transform.");
      }
    }

    delete[] tmpMatches;
    delete[] tmpScores;
    delete[] tmpScores2nd;

    return;

  }

// TODO: remove (unused)
// TODO2: scove_avg/min/max compute add to the loop
// in GetOffsets above; same with verbose messages
template < typename ImageType >
void BlockMatcherGPU<ImageType>
::PostProcessOffsets(
    const IndexListType &blockPositions, const IndexListType &blockMatchPositions,
    const std::vector<double> &blockMatchScores, const std::vector<double> &blockMatchScores2ndBest,
    PointListType &blockPositionsOut, PointListType &blockMatchPositionsOut ) const
{
  //Copy all elements up to cutoff_similarity out - because we're selecting particular elements
  //We have to do this manually
  blockMatchPositionsOut.clear();
  blockMatchPositionsOut.reserve( blockMatchScores.size() );
  blockPositionsOut.clear();
  blockPositionsOut.reserve( blockMatchScores.size() );

  int nLowScores = 0, nMultimatch = 0;
  double score_avg=0.f, score_min = 1e200, score_max = -1;

  //std::cout << "Blocks: " << std::endl;
  for( unsigned int currentIndex = 0; currentIndex<blockMatchScores.size(); ++currentIndex )
  {
    double currentBestScore = blockMatchScores[currentIndex];
    double current2ndBestScore = blockMatchScores2ndBest[currentIndex];

    //If we don't have multiple good matches, keep the block
    if( currentBestScore > 1e-7 )
    {

      if( fabs( currentBestScore - current2ndBestScore ) > 1e-4 )
      {
        //Note these are absolute positions
        blockPositionsOut.push_back( this->blockPositionToPoint( blockPositions[currentIndex] ) );

        //Note these are absolute positions
        blockMatchPositionsOut.push_back( this->blockPositionToPoint( blockMatchPositions[currentIndex] ) );
        score_avg += currentBestScore;
        score_min = std::min(currentBestScore,score_min);
        score_max = std::max(currentBestScore,score_max);
      }
      else {
        nMultimatch++;
      }
    }
    else
      nLowScores++;
  }

  if( this->m_Verbosity >= 2 )
  {
    std::cout
    << "Blocks: Initial: " << blockPositions.size()
    << ", Cull " << 100*(1.0-this->m_PortionMatchesKept) << "% low variance (<="
    << 0 << "): " << blockPositions.size()
    << ", Cull " << nLowScores << " dissimilar: " << blockPositions.size()-nLowScores
    << ", Cull " << nMultimatch << " multi-match: " <<  blockPositionsOut.size();

    if( !blockPositionsOut.empty() )
      std::cout << "  Scores: " << score_min << " to " << score_max << " (" << score_avg/blockPositionsOut.size() << ")" << std::endl;
      //INFO_INLINE("  Scores: " << score_min << " to " << score_max << " (" << score_avg/blockPositionsOut.size() << ")" );
    else
      std::cout << "  Scores: " << "None" << std::endl;
  }

  if( blockPositions.size() == 0 ) //DRH - Cannot check for blockPositionsOut.size() here because this is likely to crash in the multithreaded implementation
  {
    throw BM_GPU_EXCEPTION_ALGO( \
            "WARNING: No block matches occurred. Images probably have" \
            " zero overlap. Try supplying a better initial transform.");
  }
}

template < typename ImageType >
void BlockMatcherGPU<ImageType>::Initialize()
{
  if( this->m_Verbosity >= 2 )
    INFO_MARK("..:: Init GPU matcher ::..");

  if( ! this->releasedAll_ )
    this->release();

  this->num_kernels = 2;
  this->SetDeviceType(GPU);

  // how many platforms are there?
  cl_uint nplat;
  this->err = clGetPlatformIDs(0, NULL, &nplat);
  if (err != CL_SUCCESS) {
      ERROR_MARK( "Failed to get number of platforms!" );
      throw BM_GPU_EXCEPTION_PLATFORM("Failed to get number of platforms!");
  }
  if( this->m_Verbosity >= 2 ) {INFO_MARK("Found " << nplat << " OpenCL platforms");}

  // get all platforms ids
  cl_platform_id *pls = new cl_platform_id[nplat];
  this->err = clGetPlatformIDs(nplat, pls, NULL);
  if (err != CL_SUCCESS) {
      ERROR_MARK( "Failed to get list of platforms!");
      throw BM_GPU_EXCEPTION_PLATFORM("Failed to get list of platforms!");
  }

  // loop over all platforms until found required type
  for (unsigned int i = 0; i < nplat; i++) {
    this->err = clGetDeviceIDs(pls[i], this->m_DeviceType, 1, &(this->device_id), NULL);
    if (err == CL_SUCCESS) {
      if( this->m_Verbosity >= 2 ) {INFO_MARK("Found required device type in platform " << nplat);}
      this->cpPlatform = pls[i];
      break;
    } else if (err != CL_DEVICE_NOT_FOUND) {
      ERROR_MARK( "Failed to create a device group!" );

      ERROR_IFEQUALPRINT(CL_INVALID_PLATFORM,this->err);
      ERROR_IFEQUALPRINT(CL_INVALID_DEVICE_TYPE,this->err);
      ERROR_IFEQUALPRINT(CL_INVALID_VALUE,this->err);
      //ERROR_IFEQUALPRINT(CL_DEVICE_NOT_FOUND,this->err);
      ERROR_IFEQUALPRINT(CL_OUT_OF_RESOURCES,this->err);
      ERROR_IFEQUALPRINT(CL_OUT_OF_HOST_MEMORY,this->err);
      throw BM_GPU_EXCEPTION_PLATFORM("Failed to create a device group!");
    }
  }
  if (err != CL_SUCCESS) {
      ERROR_MARK( "Failed to create a group with the required device!" );
      throw BM_GPU_EXCEPTION_PLATFORM("Failed to create a group with the required device!");
  }
  delete[] pls;

  this->checkDevice();

  this->context = clCreateContext(0, 1, &(this->device_id), NULL, NULL, &(this->err));
  if (!context) {
      ERROR_MARK( "Failed to create a compute context!" );
      throw BM_GPU_EXCEPTION_PLATFORM("Failed to create a compute context!");
  }

  // required to measure kernel execution times
  cl_command_queue_properties profiling = CL_QUEUE_PROFILING_ENABLE;

  this->commands = clCreateCommandQueue(this->context, this->device_id, profiling, &(this->err));
  if (!commands) {
      ERROR_MARK( "Failed to create a command commands!" );
      throw BM_GPU_EXCEPTION_PLATFORM("Failed to create a command commands!");
  }

  // Create the compute program from the source buffer
  char ** liste = this->loadKernels();
  
 /* The new version uses only one kernel called from host code
  * because means are calculated by it as they are needed
  */
  this->program = clCreateProgramWithSource(this->context, 1,
      (const char **) liste,
      NULL, &(this->err));
  if (!program) {
      ERROR_MARK( "Failed to create compute program!" );
      throw BM_GPU_EXCEPTION_PLATFORM("Failed to create compute program!");
  }

  // Build the program executable
  this->err = clBuildProgram(this->program, 0, NULL, "-Werror", NULL, NULL);
  if (this->err != CL_SUCCESS) {
      size_t len;
      char buffer[2048];
      ERROR_MARK( "Failed to build program executable!" );

      ERROR_IFEQUALPRINT(CL_INVALID_PROGRAM,this->err);
      ERROR_IFEQUALPRINT(CL_INVALID_VALUE,this->err);
      ERROR_IFEQUALPRINT(CL_INVALID_DEVICE,this->err);
      ERROR_IFEQUALPRINT(CL_INVALID_BINARY,this->err);
      ERROR_IFEQUALPRINT(CL_INVALID_BUILD_OPTIONS,this->err);
      ERROR_IFEQUALPRINT(CL_INVALID_OPERATION,this->err);
      ERROR_IFEQUALPRINT(CL_COMPILER_NOT_AVAILABLE,this->err);
      ERROR_IFEQUALPRINT(CL_BUILD_PROGRAM_FAILURE,this->err);
      ERROR_IFEQUALPRINT(CL_OUT_OF_HOST_MEMORY,this->err);

      clGetProgramBuildInfo(this->program, this->device_id, CL_PROGRAM_BUILD_LOG,
          sizeof(buffer), buffer, &len);
      ERROR_MARK( "Error message:\n" << buffer );
      throw BM_GPU_EXCEPTION_PLATFORM("Failed to build program executable!");
  }

  // Create the compute kernel in the program
  /* The new version uses only one kernel called from host code
   * because means are calculated by it as they are needed
   */
  this->kernels[kMatchID] = clCreateKernel(this->program, "optKmatch", &err);
  if (err != CL_SUCCESS) {
      ERROR_IFEQUALPRINT(CL_INVALID_PROGRAM,this->err);
      ERROR_IFEQUALPRINT(CL_INVALID_PROGRAM_EXECUTABLE,this->err);
      ERROR_IFEQUALPRINT(CL_INVALID_KERNEL_NAME,this->err);
      ERROR_IFEQUALPRINT(CL_INVALID_KERNEL_DEFINITION,this->err);
      ERROR_IFEQUALPRINT(CL_INVALID_VALUE,this->err);
      ERROR_IFEQUALPRINT(CL_OUT_OF_RESOURCES,this->err);
      ERROR_IFEQUALPRINT(CL_OUT_OF_HOST_MEMORY,this->err);

      throw BM_GPU_EXCEPTION_PLATFORM("Failed to create kernel!");
  }

  this->releasedAll_ = false;
}

template < typename ImageType >
void BlockMatcherGPU<ImageType>::updateBuffers()
{

  if( !( this->m_BaseImage && this->m_SearchImage && this->m_BaseMask && this->m_SearchMask) )
    return;

  SizeType sizs = this->m_SearchImage->GetLargestPossibleRegion().GetSize();

  // maximum number of blocks in starting image, along each dimension
  size_t nblocks[3];
  nblocks[0] = (sizs[0] - this->m_BlockWidth) / this->m_BlockGap + 1;
  nblocks[1] = (sizs[1] - this->m_BlockWidth) / this->m_BlockGap + 1;
  nblocks[2] = (sizs[2] - this->m_BlockWidth) / this->m_BlockGap + 1;

  // check if we need to allocate memory (no if the old buffer are of the same size)
  bool allocateMemory = true;
#ifndef GPU_BM_USE_MASK
  if ( this->input_L_ && this->input_R_
      && this->output_Matches && this->output_Scores && this->output_Scores2nd
      && (this->region[0]    == sizs[0])    && (this->region[1]    == sizs[1])    && (this->region[2] == sizs[2])
      && (this->blockGrid[0] == nblocks[0]) && (this->blockGrid[1] == nblocks[1]) && (this->blockGrid[2] == nblocks[2]))
  {
    allocateMemory = false;
  }
#else
  if ( this->input_L_ && this->input_R_ && this->input_ML_ && this->input_MR_
      && this->output_Matches && this->output_Scores && this->output_Scores2nd
      && (this->region[0]    == sizs[0])    && (this->region[1]    == sizs[1])    && (this->region[2] == sizs[2])
      && (this->blockGrid[0] == nblocks[0]) && (this->blockGrid[1] == nblocks[1]) && (this->blockGrid[2] == nblocks[2]))
  {
    allocateMemory = false;
  }
#endif

  //allocateMemory = true;
  if (allocateMemory)
  {
    if( ! this->releasedBuffers_ ) // Must be first call in this section
      this->releaseBuffers();

    this->SetPadding(1);

    this->region[0] = sizs[0];
    this->region[1] = sizs[1];
    this->region[2] = sizs[2];

    this->blockGrid[0] = nblocks[0];
    this->blockGrid[1] = nblocks[1];
    this->blockGrid[2] = nblocks[2];

    createImage3D(&(this->input_L_), CL_MEM_READ_ONLY, region[0],region[1],region[2], CL_FLOAT, NULL);
    createImage3D(&(this->input_R_), CL_MEM_READ_ONLY, region[0],region[1],region[2], CL_FLOAT, NULL);
#ifdef GPU_BM_USE_MASK
    createImage3D(&(this->input_ML_), CL_MEM_READ_ONLY, region[0],region[1],region[2], CL_UNSIGNED_INT8, NULL);
    createImage3D(&(this->input_MR_), CL_MEM_READ_ONLY, region[0],region[1],region[2], CL_UNSIGNED_INT8, NULL);
#endif

    size_t block_count = nblocks[0] * nblocks[1] * nblocks[2];
    createBuffer(&(this->output_Scores), CL_MEM_WRITE_ONLY, sizeof(float) * block_count, NULL);
    createBuffer(&(this->output_Scores2nd), CL_MEM_WRITE_ONLY, sizeof(float) * block_count, NULL);
    createBuffer(&(this->output_Matches), CL_MEM_WRITE_ONLY, sizeof(unsigned int) * block_count, NULL);

    this->releasedBuffers_ = false;

    if (!input_L_ || !input_R_ || !output_Scores || !output_Scores2nd || !output_Matches)
    {
        ERROR_MARK( "Failed to allocate device memory!" );
        throw BM_GPU_EXCEPTION_PLATFORM("Failed to allocate device memory!");
    }
#ifdef GPU_BM_USE_MASK
    if (!this->input_ML_ || !this->input_MR_)
    {
      ERROR_MARK( "Failed to allocate device memory! (mask)" );
      throw BM_GPU_EXCEPTION_PLATFORM("Failed to allocate device memory! (mask)");
    }
#endif
  }

  // copy image data to OpenCL buffers
  size_t origin[3];
  origin[0] = origin[1] = origin[2] = 0;
  // row and slice pitch could actually be 0 for automatic calculation by openCL runtime?
  this->err  = clEnqueueWriteImage(this->commands, this->input_L_, CL_TRUE, origin, region,
           sizeof(float) * region[0],
           sizeof(float) * region[0] * region[1],
                 this->m_BaseImage->GetBufferPointer(),
           0, NULL, NULL);
  this->err |= clEnqueueWriteImage(this->commands, this->input_R_, CL_TRUE, origin, region,
           sizeof(float) * region[0],
           sizeof(float) * region[0] * region[1],
           this->m_SearchImage->GetBufferPointer(),
           0, NULL, NULL);
#ifdef GPU_BM_USE_MASK
  this->err |= clEnqueueWriteImage(this->commands, this->input_ML_, CL_TRUE, origin, region,
           sizeof(unsigned char) * region[0],
           sizeof(unsigned char) * region[0] * region[1],
           this->m_BaseMask->GetBufferPointer(),
           0, NULL, NULL);
  this->err |= clEnqueueWriteImage(this->commands, this->input_MR_, CL_TRUE, origin, region,
           sizeof(unsigned char) * region[0],
           sizeof(unsigned char) * region[0] * region[1],
           this->m_SearchMask->GetBufferPointer(),
           0, NULL, NULL);
#endif
}

template < typename ImageType >
void BlockMatcherGPU<ImageType>::createBuffer(cl_mem * buffer, cl_mem_flags flags,size_t size, void *host_ptr)
{
  cl_int errcode_ret;
  //if( * buffer )
  //  clReleaseMemObject(*buffer);
  if( this->context )
    {
      *buffer = clCreateBuffer(this->context, flags, size, host_ptr, &errcode_ret);
      if( errcode_ret != CL_SUCCESS )
        {
          ERROR_IFPRINT(CL_INVALID_CONTEXT,errcode_ret);
          ERROR_IFPRINT(CL_INVALID_VALUE,errcode_ret);
          ERROR_IFPRINT(CL_INVALID_BUFFER_SIZE,errcode_ret);
          ERROR_IFPRINT(CL_INVALID_HOST_PTR,errcode_ret);
          ERROR_IFPRINT(CL_MEM_OBJECT_ALLOCATION_FAILURE,errcode_ret);
          ERROR_IFPRINT(CL_OUT_OF_HOST_MEMORY,errcode_ret);

          throw BM_GPU_EXCEPTION_PLATFORM("Impossible to create a CL buffer.");
        }
    }
  else
    throw BM_GPU_EXCEPTION_PLATFORM("No context for createBuffer.");
}

/*
 * cl_channel_type = CL_FLOAT, CL_UNSIGNED_INT8...
 */
template < typename ImageType >
void BlockMatcherGPU<ImageType>::createImage3D(cl_mem * buffer, cl_mem_flags flags,
             size_t w, size_t h, size_t d, cl_channel_type data_type,
           void *host_ptr)
{
  cl_int errcode_ret;
  cl_image_format i_format ;
  i_format.image_channel_data_type = data_type; //== CL_FLOAT, CL_UNSIGNED_INT8...
  i_format.image_channel_order = CL_INTENSITY;

  if( this->context ) {
    *buffer = clCreateImage3D(this->context, flags, &i_format, w, h, d, 0, 0, host_ptr,
    &errcode_ret);
    if( errcode_ret != CL_SUCCESS ) {
      ERROR_IFPRINT(CL_INVALID_CONTEXT,errcode_ret);
      ERROR_IFPRINT(CL_INVALID_VALUE,errcode_ret);
      ERROR_IFPRINT(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR,errcode_ret);
      ERROR_IFPRINT(CL_INVALID_IMAGE_SIZE,errcode_ret);
      ERROR_IFPRINT(CL_IMAGE_FORMAT_NOT_SUPPORTED,errcode_ret);
      ERROR_IFPRINT(CL_INVALID_OPERATION,errcode_ret);
      ERROR_IFPRINT(CL_INVALID_HOST_PTR,errcode_ret);
      ERROR_IFPRINT(CL_MEM_OBJECT_ALLOCATION_FAILURE,errcode_ret);
      ERROR_IFPRINT(CL_OUT_OF_RESOURCES,errcode_ret);
      ERROR_IFPRINT(CL_OUT_OF_HOST_MEMORY,errcode_ret);

      throw BM_GPU_EXCEPTION_PLATFORM("Impossible to create a CL 3D image.");
    }
  }
  else
    throw BM_GPU_EXCEPTION_PLATFORM("No context for createImage3D.");
}

template < typename ImageType >
void BlockMatcherGPU<ImageType>::checkDevice()
{
  if( this->m_Verbosity >= 2 )
    {

      INFO_MARK("..:: Check the device ::..");
      size_t 						p_size;
      size_t 						ret_size;			// returned size for type infos
      char 							param[1024];			// returned arg for type infos
      cl_device_type 		dev_type;			// type infos of the device
      cl_uint 					entries;
      size_t 						arr_tsize[3];
      cl_ulong 					long_entries;
      cl_bool 					bool_entries;
      cl_device_local_mem_type 			mem_type;
      cl_device_fp_config 					fp_conf;
      cl_device_exec_capabilities 	exec_cap;

      clGetDeviceInfo(device_id,CL_DEVICE_TYPE,sizeof(dev_type),&dev_type,&ret_size);
      INFO_MARK("Device Type :");
      INFO_IFPRINT(CL_DEVICE_TYPE_GPU,dev_type);
      INFO_IFPRINT(CL_DEVICE_TYPE_CPU,dev_type);
      INFO_IFPRINT(CL_DEVICE_TYPE_ACCELERATOR,dev_type);
      INFO_IFPRINT(CL_DEVICE_TYPE_DEFAULT,dev_type);


      clGetDeviceInfo(device_id,CL_DEVICE_NAME,sizeof(param),param,&ret_size);
      INFO_MARK("Name : \t\t\t\t" << param);

      clGetDeviceInfo(device_id,CL_DEVICE_VENDOR,sizeof(param),param,&ret_size);
      INFO_MARK("Vendor : \t\t\t" << param);

      clGetDeviceInfo(device_id,CL_DEVICE_VENDOR_ID,sizeof(cl_uint),&entries,&ret_size);
      INFO_MARK("Vendor ID : \t\t\t" << entries);

      clGetDeviceInfo(device_id,CL_DEVICE_VERSION,sizeof(param),param,&ret_size);
      INFO_MARK("Version : \t\t\t" << param);

      clGetDeviceInfo(device_id,CL_DEVICE_PROFILE,sizeof(param),param,&ret_size);
      INFO_MARK("Profile : \t\t\t" << param);

      clGetDeviceInfo(device_id,CL_DRIVER_VERSION,sizeof(param),param,&ret_size);
      INFO_MARK("Driver : \t\t\t" << param);

      clGetDeviceInfo(device_id,CL_DEVICE_EXTENSIONS,sizeof(param),param,&ret_size);
      INFO_MARK("Extensions : \t\t\t" << param);

      clGetDeviceInfo(device_id,CL_DEVICE_MAX_WORK_ITEM_SIZES,3*sizeof(size_t),arr_tsize,&ret_size);
      INFO_MARK("Max Work-Item Sizes : \t\t(" << arr_tsize[0] << "," << arr_tsize[1] << "," << arr_tsize[2] << ")");

      clGetDeviceInfo(device_id,CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(size_t),&p_size,&ret_size);
      INFO_MARK("Max Work Group Size:\t\t"<<p_size);

      clGetDeviceInfo(device_id,CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(cl_uint),&entries,&ret_size);
      INFO_MARK("Max Compute Units:\t\t"<<entries);

      clGetDeviceInfo(device_id,CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(cl_uint),&entries,&ret_size);
      INFO_MARK("Max Frequency (Mhz):\t\t"<<entries);

      clGetDeviceInfo(device_id,CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE,sizeof(cl_uint),&entries,&ret_size);
      INFO_MARK("Cache Line (bytes):\t\t"<<entries);

      clGetDeviceInfo(device_id,CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,&ret_size);
      INFO_GATHER(
          INFO_INLINE("Global Memory (MB):\t\t");
      printf("%llu",(long long unsigned int)(long_entries/1024/1024));
      );
      INFO_MARK("Global Memory (bytes):\t\t"<<long_entries);

      clGetDeviceInfo(device_id,CL_DEVICE_LOCAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,&ret_size);
      INFO_GATHER(
          INFO_INLINE("Local Memory (MB):\t\t");
      printf("%llu",(long long unsigned int)(long_entries/1024/1024));
      );
      INFO_MARK("Local Memory (bytes):\t\t"<<long_entries);

      clGetDeviceInfo(device_id,CL_DEVICE_LOCAL_MEM_TYPE,sizeof(cl_device_local_mem_type),&mem_type,&ret_size);
      INFO_MARK("Local Memory Type : ");
      INFO_IFPRINT(CL_LOCAL,mem_type);
      INFO_ELIFPRINT(CL_GLOBAL,mem_type);
      INFO_ELSEPRINT(UNKNOWN);

      clGetDeviceInfo(device_id,CL_DEVICE_MAX_MEM_ALLOC_SIZE,sizeof(cl_ulong),&long_entries,&ret_size);
      INFO_GATHER(
          INFO_INLINE("Max Mem Alloc (MB):\t\t");
      printf("%llu",(long long unsigned int)(long_entries/1024/1024));
      );
      INFO_GATHER(
          INFO_INLINE("Max Mem Alloc (B):\t\t");
      printf("%llu",(long long unsigned int)(long_entries));
      );

      clGetDeviceInfo(device_id,CL_DEVICE_MAX_PARAMETER_SIZE,sizeof(size_t),&p_size,&ret_size);
      INFO_MARK("Max Param Size (MB):\t\t"<<p_size);

      clGetDeviceInfo(device_id,CL_DEVICE_MEM_BASE_ADDR_ALIGN,sizeof(cl_uint),&entries,&ret_size);
      INFO_MARK("Base Mem Align (bits):\t\t"<<entries);

      clGetDeviceInfo(device_id,CL_DEVICE_ADDRESS_BITS,sizeof(cl_uint),&entries,&ret_size);
      INFO_MARK("Address Space (bits):\t\t"<<entries);

      clGetDeviceInfo(device_id,CL_DEVICE_IMAGE_SUPPORT,sizeof(cl_bool),&bool_entries,&ret_size);
      INFO_MARK("Image Support:\t\t\t"<<bool_entries);

      clGetDeviceInfo(device_id,CL_DEVICE_TYPE,sizeof(fp_conf),&fp_conf,&ret_size);
      INFO_MARK("Float Functionality:");
      INFO_IFPRINT(CL_FP_DENORM,fp_conf);
      INFO_IFPRINT(CL_FP_ROUND_TO_NEAREST,fp_conf);
      INFO_IFPRINT(CL_FP_ROUND_TO_ZERO,fp_conf);
      INFO_IFPRINT(CL_FP_ROUND_TO_INF,fp_conf);
      INFO_IFPRINT(CL_FP_FMA,fp_conf);
      INFO_IFPRINT(CL_FP_INF_NAN,fp_conf);

      clGetDeviceInfo(device_id,CL_DEVICE_ERROR_CORRECTION_SUPPORT,sizeof(cl_bool),&bool_entries,&ret_size);
      INFO_MARK("ECC Support:\t\t\t"<<bool_entries);

      clGetDeviceInfo(device_id,CL_DEVICE_EXECUTION_CAPABILITIES,sizeof(cl_device_exec_capabilities),&exec_cap,&ret_size);
      INFO_MARK("Exec Functionality:");
      INFO_IFPRINT(CL_EXEC_KERNEL,exec_cap);
      INFO_IFPRINT(CL_EXEC_NATIVE_KERNEL,exec_cap);

      clGetDeviceInfo(device_id,CL_DEVICE_ENDIAN_LITTLE,sizeof(cl_bool),&bool_entries,&ret_size);
      INFO_MARK("Little Endian Device:\t\t"<<bool_entries);

      clGetDeviceInfo(device_id,CL_DEVICE_PROFILING_TIMER_RESOLUTION,sizeof(size_t),&p_size,&ret_size);
      INFO_MARK("Profiling Res (ns):\t\t"<<p_size);

      clGetDeviceInfo(device_id,CL_DEVICE_AVAILABLE,sizeof(cl_bool),&bool_entries,&ret_size);
      INFO_MARK("Device Available:\t\t"<<bool_entries);

      clGetDeviceInfo(device_id,CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,sizeof(cl_ulong),&long_entries,&ret_size);
      INFO_MARK("CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE : " << long_entries);
    }
}


template < typename ImageType >
char ** BlockMatcherGPU<ImageType>::loadKernels()
{
  if( this->m_Verbosity >= 2 )
    INFO_MARK("..:: Load OpenCL Kernel ::..");
  char **liste = new char*[1];

  // ---------------------------------------------------------------

  liste[0] = new char[newKmatchString.size()+1];
  strncpy(liste[0],newKmatchString.c_str(),newKmatchString.size());
  liste[0][newKmatchString.size()] = '\0';

  return liste;
}

template < typename ImageType >
void BlockMatcherGPU<ImageType>::SetBaseImage( const ImageType * movingImage )
{
  Superclass::SetBaseImage(movingImage);
}

template < typename ImageType >
void BlockMatcherGPU<ImageType>::SetSearchImage( const ImageType * fixedImage )
{
  Superclass::SetSearchImage(fixedImage);
}

template < typename ImageType >
void BlockMatcherGPU<ImageType>::SetDeviceType(const BlockMatcherGPU<ImageType>::devType t)
{
  if( this->m_Verbosity >= 2 )
    INFO_MARK("..:: SetDeviceType ::..");
  if( t == GPU )
    this->m_DeviceType = CL_DEVICE_TYPE_GPU;
  else
    this->m_DeviceType = CL_DEVICE_TYPE_CPU;
}

template < typename ImageType >
void BlockMatcherGPU<ImageType>::displayInfo()
{
  Superclass::displayInfo();

  std::cout << "GPU normalized_correlation" << std::endl;
}

template < typename ImageType >
void BlockMatcherGPU<ImageType>::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

template < typename ImageType >
void BlockMatcherGPU<ImageType>::getImageFormatsInfo( cl_context ctxt, cl_mem_flags m, cl_mem_object_type te)
{
  if( this->m_Verbosity >= 2 )
    {
      INFO_MARK("..:: getImageFormatsInfo ::..");
      cl_uint num_entries;  cl_image_format *image_formats;
      cl_int status=clGetSupportedImageFormats (ctxt,m,te,0,NULL,&num_entries);
      if(status==CL_SUCCESS && num_entries>0)
        {
          image_formats=(cl_image_format*)malloc(num_entries*sizeof(cl_image_format));
          status=clGetSupportedImageFormats (ctxt,m,te,num_entries,image_formats,NULL);
          if(status == CL_SUCCESS)
            {
              int o,t;
              int i,j;
              cl_int orders[] = {CL_R,  CL_A,CL_INTENSITY, CL_LUMINANCE,CL_RG,  CL_RA,CL_RGB,CL_RGBA,CL_ARGB, CL_BGRA};
              std::string or2[] = {"CL_R", "CL_A","CL_INTENSITY", "CL_LUMINANCE","CL_RG", "CL_RA","CL_RGB","CL_RGBA","CL_ARGB", "CL_BGRA" };
              cl_int types[]= {
                  CL_SNORM_INT8 , CL_SNORM_INT16, CL_UNORM_INT8, CL_UNORM_INT16, CL_UNORM_SHORT_565, CL_UNORM_SHORT_555, CL_UNORM_INT_101010,CL_SIGNED_INT8,
                  CL_SIGNED_INT16,  CL_SIGNED_INT32, CL_UNSIGNED_INT8, CL_UNSIGNED_INT16, CL_UNSIGNED_INT32, CL_HALF_FLOAT, CL_FLOAT};
              std::string tt[] = {"CL_SNORM_INT8" ,"CL_SNORM_INT16","CL_UNORM_INT8","CL_UNORM_INT16","CL_UNORM_SHORT_565","CL_UNORM_SHORT_555","CL_UNORM_INT_101010",
                  "CL_SIGNED_INT8","CL_SIGNED_INT16","CL_SIGNED_INT32","CL_UNSIGNED_INT8","CL_UNSIGNED_INT16","CL_UNSIGNED_INT32","CL_HALF_FLOAT","CL_FLOAT"};

              for(i=0; i<num_entries; i++)
                {
                  for(j=0; j<sizeof(orders)/sizeof(orders[0]); j++)
                    {
                      if(image_formats[i].image_channel_order==orders[j])
                        o=j;
                    }
                  for(j=0; j<sizeof(types)/sizeof(orders[0]); j++)
                    {
                      if(image_formats[i].image_channel_data_type==types[j])
                        t=j;
                    }
                  INFO_MARK("Format " << i << " : " << or2[o] << ", " << tt[t]);
                }
            }
          free(image_formats);
        }
    }
  }
}

#endif
