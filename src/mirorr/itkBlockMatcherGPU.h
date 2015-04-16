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
#ifndef __itkBlockMatcherGPU_h
#define __itkBlockMatcherGPU_h

#include <iostream>
#include <fstream>
#include <itkImage.h>
#include <itkSize.h>
#include <itkIndent.h>
#include <itkExceptionObject.h>
#include <vector>
#include <sstream>
#include <string>
#include <stdexcept>
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
#include "itkAbstractBlockMatcher.h"

///////////////////////////////
//
#define GPU_BM_USE_MASK
//
///////////////////////////////

namespace itk
{
  class BlockMatcherGPUException {
  public:
    BlockMatcherGPUException() {msg = "Undefined GPU BG exception";}
    BlockMatcherGPUException(const std::string _msg) {this->msg = _msg;}
    std::string GetMessage() const {return this->msg;}
  protected:
    std::string msg;
  };

  class BlockMatcherGPUExceptionAlgo : public BlockMatcherGPUException{
    // Exception related to the BM algorithm. If things fail here,
    // the CPU implementation should normally fail at the same point
  public:
    BlockMatcherGPUExceptionAlgo() {msg = "Undefined GPU BM algo exception";}
    BlockMatcherGPUExceptionAlgo(const std::string _msg) {this->msg = _msg;}
  };

  class BlockMatcherGPUExceptionPlatform : public BlockMatcherGPUException{
    // Exception related to the GPU platform implementation...
    // if things fail here, it might be worth using the CPU implementation
  public:
    BlockMatcherGPUExceptionPlatform() {msg = "Undefined GPU BM algo exception";}
    BlockMatcherGPUExceptionPlatform(const std::string _msg) {this->msg = _msg;}
  };

  template < typename ImageType >
  class BlockMatcherGPU : public AbstractBlockMatcher<ImageType >
  {
  public:

    ///////////////////////////////////////
    //        ITK standards
    /** Standard class typedefs. */
    typedef BlockMatcherGPU                       Self;
    typedef AbstractBlockMatcher<ImageType >      Superclass;
    typedef SmartPointer<Self>                    Pointer;
    typedef SmartPointer<const Self>              ConstPointer;
    itkStaticConstMacro(ImageDimension, unsigned int,
        ImageType::ImageDimension);

    //typedef itk::OrientedImage<unsigned char,3> MaskImageType;
    //typedef itk::Image<unsigned char,3> MaskImageType;
    //typedef typename MaskImageType::Pointer     MaskImagePointer;
    typedef typename Superclass::MaskType MaskType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ImageRegistrationMethod, ProcessObject);
    //
    ///////////////////////////////////////

    ///////////////////////////////////////
    //      Types
    enum 	devType 						// Types of device available
    { CPU, GPU };
    typedef typename ImageType::SizeType    SizeType;
    typedef Index<ImageDimension>           IndexType;
    typedef std::vector< IndexType >        IndexListType;
    typedef typename ImageType::PointType   PointType;
    typedef std::vector< PointType >        PointListType;
    typedef itk::ImageRegion<ImageType::ImageDimension>      RegionType;
    typedef itk::SparseImageRegionConstIterator< ImageType > SparseConstIterator;
    typedef itk::ImageRegionConstIterator< ImageType >       ConstIterator;
    //
    ///////////////////////////////////////

    ///////////////////////////////////////
    //      Images
    void SetBaseImage( const ImageType * movingImage );
    void SetSearchImage( const ImageType * fixedImage );
    //
    ///////////////////////////////////////

    ///////////////////////////////////////
    //      Algorithm Parameters
    void SetDeviceType(const BlockMatcherGPU<ImageType>::devType t);
    //
    ///////////////////////////////////////

    ///////////////////////////////////////
    // Allows to perform block matching in the other direction
    // Swaps m_BaseImage <-> m_SearchImage,
    //       m_BaseMask  <-> m_SearchMask,
    // Does NOT allocate OpenCL buffer: use with caution...
    virtual void SwapImageOrder();
    ///////////////////////////////////////

    ///////////////////////////////////////
    //      Match function
    void GetOffsets( PointListType &blockPositions,
        PointListType &blockMatchPositions );
    //
    ///////////////////////////////////////

  /**
   * Filter match lists by removing low-score or ambiguous matches
   * (copied from itkBlockMatcher.txx)
   */
  virtual void PostProcessOffsets(const IndexListType &blockPositions, const IndexListType &blockMatchPositions,
      const std::vector<double> &blockMatchScores, const std::vector<double> &blockMatchScores2ndBest,
      PointListType &blockPositionsOut, PointListType &blockMatchPositionsOut) const;

    ///////////////////////////////////////
    //      Displays
    void displayInfo();
    void PrintSelf(std::ostream& os, Indent indent) const;
    //
    ///////////////////////////////////////

    ///////////////////////////////////////
    //      Metrics
    virtual void SetBlockMetricType(std::string )
    {
      //Dummy function as the GPU implementation always compute normalized_correlation
    }
    //
    ///////////////////////////////////////

  private:
    // Algorithm parameters
    //int 			count;
    size_t region[3];
    size_t blockGrid[3];

    // OpenCL parameters
    int			m_DeviceType;		// Type of device to use (GPU/CPU/...)
    cl_int 		err;     		// error code returned from api calls
    cl_platform_id 	cpPlatform; 		// OpenCL platform
    cl_device_id 		device_id;    		// compute device id
    cl_context 		context;        	// compute context
    cl_command_queue 	commands; 		// compute command queue
    cl_program 		program;        	// compute program
/* The new version uses only one kernel called from host code
 * because means are calculated by it as they are needed
 */
    cl_kernel 		kernels[1];             // compute kernel
    cl_uint 		num_kernels;		// number of kernels
    // device memory used for the input array
    cl_mem 		input_L_;		//
    cl_mem 		input_R_;		//
#ifdef GPU_BM_USE_MASK
    cl_mem    input_ML_;
    cl_mem    input_MR_;
#endif
    // device memory used for the output array
    cl_mem 		output_L_kMean;		//
    cl_mem 		output_R_kMean;		//
    cl_mem 		output_L_kVariance;
    cl_mem 		output_R_kVariance;
    cl_mem		output_Matches;		//
    cl_mem		output_Scores;		//
    cl_mem    output_Scores2nd;       //

    bool      releasedAll_;
    bool      releasedBuffers_;
    bool      recycleImageBufferData_; //DRH: optimisation hack - do not release OpenCL buffer in next GetOffset() call

    ///////////////////////////////////////
    //      Constructors
    BlockMatcherGPU();
    ~BlockMatcherGPU();
    void Initialize();
    void release();
    void releaseBuffers();
    //
    ///////////////////////////////////////

    ///////////////////////////////////////
    //      Buffer managers
    void updateBuffers();
    void createBuffer(cl_mem * buffer, cl_mem_flags flags, size_t size,void *host_ptr);
    void createImage3D(cl_mem * buffer, cl_mem_flags flags, size_t w, size_t h, size_t d, cl_channel_type data_type, void *host_ptr);
    //
    ///////////////////////////////////////

    ///////////////////////////////////////
    //      Mean Variance computation
    void computeMeanVarianceImages();
    //
    ///////////////////////////////////////

    ///////////////////////////////////////
    //      GPU device tools
    void checkDevice();
    char ** loadKernels();
    void getImageFormatsInfo(cl_context context,cl_mem_flags m,cl_mem_object_type te);
    //
    ///////////////////////////////////////

  };

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatcherGPU.txx"
#endif

#endif
