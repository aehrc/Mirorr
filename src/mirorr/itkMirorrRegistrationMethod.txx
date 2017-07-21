/*=========================================================================
Program: mirorr
Module: itkMirorrRegistrationMethod.txx
Author: Nicholas Dowson
Created: 20 Feb 2009

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
/**
* A class for registering images using block matching methods
*/

#ifndef __itkMirorrRegistrationMethod_txx
#define __itkMirorrRegistrationMethod_txx

#include "itkMirorrRegistrationMethod.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkBlockMatcher.h"
#include "itkBlockMatcherThreaded.h"
#include "itkEuler3DTransform.h"

#ifdef USE_OPENCL
  #include "itkBlockMatcherGPU.h"
#endif
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_trace.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_quaternion.h>

#include <boost/timer/timer.hpp>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

namespace itk
{
/**
* Constructor
*/
template < typename TMovingImage, typename TFixedImage >
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::MirorrRegistrationMethod()
 {
  //imageFileName = "";

  this->SetNumberOfRequiredOutputs( 1 );  // for the Transform

  m_MovingImage   = 0; // has to be provided by the user.
  m_FixedImage  = 0; // has to be provided by the user.
  m_Transform    = 0; // has to be provided by the user.
  m_Interpolator = 0; // has to be provided by the user.
  m_MaskInterpolator = 0; // has to be provided by the user.
  //m_Metric       = MetricType::New(); // has to be provided by the user.
  //m_Optimizer    = 0; // has to be provided by the user.
  m_Verbosity    = 1; //Be somewhat verbose by default
  m_PortionMatchesKept = 0.5; //Default //Keep 25% of matches
  m_BlockMetricType = ""; // Use block matcher default

  m_MaxIterations = 5; //Maximum no. iterations

#ifdef USE_OPENCL
  m_UseGpuOn      = false;
#endif

  //m_UseSymResampling = true;
  //m_UseMaxResampling = false;
  m_ResamplingType = EResamplingMiddle;

  //m_InitialTransformParameters = ParametersType(1);
  m_LastTransformParameters = ParametersType(1);

  //m_InitialTransformParameters.Fill( 0.0f );
  m_LastTransformParameters.Fill( 0.0f );

  m_MovingImageRegionDefined = false;

  m_UseMultiThreading = true;

  //Default Block Matcher Parameters
  m_NhoodWidth        = 7; //Blocks
  m_NhoodGap          = 3; //Pixels
  m_BlockWidth        = 4; //Pixels
  m_BlockGap          = 1; //Pixels
#ifdef USE_NPW
  m_NPWbins           = 32;
#ifndef USE_OPENCL
  m_NPWscaleFactor    = 2;
  m_NPWshapeExpansion = false;
#endif
#endif

  TransformOutputPointer transformDecorator =
  static_cast< TransformOutputType * >(
  this->MakeOutput(0).GetPointer() );

  this->ProcessObject::SetNthOutput( 0, transformDecorator.GetPointer() );

#ifdef ITK_USE_OPTIMIZED_REGISTRATION_METHODS
  this->SetNumberOfThreads( this->GetMultiThreader()->GetNumberOfThreads() );
#else
  this->SetNumberOfThreads( 1 );
  this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
#endif

  this->m_NumberOfBlockMatcherThreads = this->GetMultiThreader()->GetGlobalDefaultNumberOfThreads();
 }

template < typename TMovingImage, typename TFixedImage >
unsigned long
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::GetMTime() const
 {
  unsigned long mtime = Superclass::GetMTime();
  unsigned long m;

  // Some of the following should be removed once ivars are put in the
  // input and output lists

  if (m_Transform)
  {
    m = m_Transform->GetMTime();
    mtime = (m > mtime ? m : mtime);
  }

  if (m_Interpolator)
  {
    m = m_Interpolator->GetMTime();
    mtime = (m > mtime ? m : mtime);
  }

  //FIXME: Check MaskInterpolator mtime

  if (m_MovingImage)
  {
    m = m_MovingImage->GetMTime();
    mtime = (m > mtime ? m : mtime);
  }

  if (m_FixedImage)
  {
    m = m_FixedImage->GetMTime();
    mtime = (m > mtime ? m : mtime);
  }

  return mtime;
 }

// Set the region of the moving image to be considered for registration
template < typename TMovingImage, typename TFixedImage >
void
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::SetMovingImageRegion( const MovingImageRegionType & region )
 {
  m_MovingImageRegion = region;
  m_MovingImageRegionDefined = true;
  this->Modified();
 }

//Initialize by setting the interconnects between components.
template < typename TMovingImage, typename TFixedImage >
void
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::Initialize() throw (ExceptionObject)
 {
  if( !m_MovingImage )
  {
    itkExceptionMacro(<<"MovingImage is not present");
  }

  if( !m_FixedImage )
  {
    itkExceptionMacro(<<"FixedImage is not present");
  }

  if( !m_Transform )
  {
    itkExceptionMacro(<<"Transform is not present");
  }

  //
  // Connect the transform to the Decorator.
  //
  TransformOutputType * transformOutput =
  static_cast< TransformOutputType * >( this->ProcessObject::GetOutput(0) );

  transformOutput->Set( m_Transform.GetPointer() );

  if( !m_Interpolator )
  {
    itkExceptionMacro(<<"Interpolator is not present");
  }
  //FIXME: Check for MaskInterpolator != 0

  //Set up the block matcher
  this->CreateBlockMatcher();
 }

template < typename TMovingImage, typename TFixedImage >
void
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::CreateBlockMatcher()
{
#ifdef USE_OPENCL
  if( m_UseGpuOn  && (m_BlockMetricType.find("npwmi") == std::string::npos)) {
    this->CreateBlockMatcherGPU();
  }
  else {
    this->CreateBlockMatcherCPU();
  }
#else
  this->CreateBlockMatcherCPU();
#endif
}

template < typename TMovingImage, typename TFixedImage >
void
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::CreateBlockMatcherCPU()
{
  if (m_BlockMatcher.IsNotNull()) {
    //Destroy any existing BlockMatcher, in case some settings have changed.
    m_BlockMatcher = 0;
  }

  if (m_UseMultiThreading && (m_BlockMetricType.find("npwmi") == std::string::npos) ) {
    //std::cout << "# Using MultiThreading" << std::endl;
    m_BlockMatcher = (BlockMatcherPointer)BlockMatcherThreaded<TMovingImage >::New();
    m_BlockMatcher->SetNumberOfThreads(this->m_NumberOfBlockMatcherThreads);
  } else {
    //std::cout << "# NOT Using Multi Threading" << std::endl;
    m_BlockMatcher = (BlockMatcherPointer)BlockMatcher<TMovingImage >::New();
  }
}

#ifdef USE_OPENCL
template < typename TMovingImage, typename TFixedImage >
void
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::CreateBlockMatcherGPU()
{
  if (m_BlockMatcher.IsNotNull()) {
    //Destroy any existing BlockMatcher, in case some settings have changed.
    m_BlockMatcher = 0;
  }

  m_BlockMatcher = (BlockMatcherPointer)BlockMatcherGPU<TMovingImage >::New();
}
#endif

template < typename TMovingImage, typename TFixedImage >
void
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::SetupBlockMatcher()
{
  m_BlockMatcher->SetBlockWidth( m_BlockWidth );  //Pixels
  m_BlockMatcher->SetBlockGap( m_BlockGap ); //Pixels
  m_BlockMatcher->SetNhoodGap( m_NhoodGap ); //Pixels

  m_BlockMatcher->SetNhoodWidth( m_NhoodWidth ); //Blocks
  m_BlockMatcher->SetVerbosity( m_Verbosity ); //Say something - same level
  m_BlockMatcher->SetPortionMatchesKept(m_PortionMatchesKept); //Set portion of matches kept

  #ifdef USE_NPW
  this->m_BlockMatcher->SetNPWbins( this->m_NPWbins );
  #ifndef USE_OPENCL
  this->m_BlockMatcher->SetNPWscaleFactor( this->m_NPWscaleFactor );
  this->m_BlockMatcher->SetNPWshapeExpansion( this->m_NPWshapeExpansion );
  #endif
  #endif

  try
  {
    m_BlockMatcher->SetBlockMetricType(m_BlockMetricType);
  }
  catch( std::exception &ee )
  {
    std::cerr << "WARNING: Problem setting block metric type to '"
        << m_BlockMetricType << "', because:"
        << ee.what()
        << "\nWARNING: Defaulting to 'normalized_correlation'."
        << std::endl;
    this->SetupBlockMatcher_tryNC();
  }
  catch( ... )
  {
    std::cerr << "WARNING: Problem setting block metric type to '"
        << m_BlockMetricType << "'. Defaulting to 'normalized_correlation'."
        << std::endl;
    this->SetupBlockMatcher_tryNC();
  }
}

template < typename TMovingImage, typename TFixedImage >
void
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::SetupBlockMatcher_tryNC()
{
  m_BlockMetricType = "normalized_correlation";

  try
  {
    m_BlockMatcher->SetBlockMetricType(m_BlockMetricType);
  }
  catch( std::exception &ee )
  {
    std::stringstream ss;
    ss << "ERROR: Unable to set BlockMetricType to '"
        << m_BlockMetricType << "', because " << ee.what()
        <<std::endl;
    throw std::runtime_error(ss.str());
  }
}

//Starts the Registration Process
template < typename TMovingImage, typename TFixedImage >
void
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::StartRegistration( void )
{

  // StartRegistration is an old API from before
  // ImageRegistrationMethod was a subclass of ProcessObject.
  // Historically, one could call StartRegistration() instead of
  // calling Update().  However, when called directly by the user, the
  // inputs to ImageRegistrationMethod may not be up to date.  This
  // may cause an unexpected behavior.
  //
  // Since we cannot eliminate StartRegistration for backward
  // compability reasons, we check whether StartRegistration was
  // called directly or whether Update() (which in turn called
  // StartRegistration()).
  if (!m_Updating)
  {
    this->Update();
  }
  else
  {
    ParametersType empty(1);
    empty.Fill( 0.0 );
    try
    {
      // initialize the interconnects between components
      this->Initialize();
    }
    catch( ExceptionObject& err )
    {
      m_LastTransformParameters = empty;

      // pass exception to caller
      throw err;
    }

    TransformPointer tBackup = this->m_Transform->Clone();

    try
    {
      this->StartOptimization();
    }
#ifdef USE_OPENCL
    catch (BlockMatcherGPUExceptionPlatform &exp)
    {
      // Try to fall back on CPU if GPU fail
      std::cerr
        << "Warning in MirorrRegistrationMethod: GPU BM implementation failed.\n"
        << "Message: " << exp.GetMessage() << std::endl
        << "Attempting to revert to a CPU implementation." << std::endl;

      this->m_Transform->SetParametersByValue(tBackup->GetParameters());
      this->m_Transform->SetFixedParameters(tBackup->GetFixedParameters());
      this->m_Transform->GetParameters(); //Force update

      this->CreateBlockMatcherCPU();
      this->StartOptimization();
    }
#endif
    catch (...)
    {
      throw;
    }
  }
}

template<typename TMovingImage, typename TFixedImage>
typename MirorrRegistrationMethod<TMovingImage, TFixedImage>::PointType
MirorrRegistrationMethod<TMovingImage, TFixedImage>::ComputeImageCenter(
    MovingImageConstPointer img)
{
  typedef ContinuousIndex< typename TMovingImage::PointType::ValueType,
      TFixedImage::ImageDimension >  ContinuousIndexType;
  typedef typename ContinuousIndexType::ValueType ContinuousIndexValueType;

  typename TMovingImage::RegionType::IndexType fixedIndex = img->GetLargestPossibleRegion().GetIndex();
  typename TMovingImage::RegionType::SizeType fixedSize = img->GetLargestPossibleRegion().GetSize();

  PointType centerPoint;
  ContinuousIndexType centerIndex;

  for ( unsigned int k = 0; k < 3; k++ ) {
    centerIndex[k] =
      static_cast< ContinuousIndexValueType >( fixedIndex[k] )
      + static_cast< ContinuousIndexValueType >( fixedSize[k] - 1 ) / 2.0;
  }

  img->TransformContinuousIndexToPhysicalPoint(centerIndex, centerPoint);

  return centerPoint;
}

/*
* Starts the Optimization process
*/
template < typename TMovingImage, typename TFixedImage >
void
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::StartOptimization( void )
 {
  if (this->m_Verbosity >= 1)
  {
    std::cout << "Using classic implementation" << std::endl;
  }

  typedef itk::ResampleImageFilter< TMovingImage,  TFixedImage > ResampleFilterType;
  typedef itk::ResampleImageFilter< MaskImageType, MaskImageType > MaskResampleFilterType;

  typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  typename MaskResampleFilterType::Pointer mask_resampler = MaskResampleFilterType::New();

  MovingImageConstPointer movingImage = 0; //m_MovingImage;
  MaskImagePointer        movingMask  = 0; //m_MovingMask;

  /**
   * DRH - Important!
   *
   * When "this->m_ResamplingType == EResamplingMiddle", the image resampling space
   * is computed by taking a mid point between the two input images. This is in opposition
   * to the previous scheme where the resampling space was the space of the fixed (or moving) image.
   *
   * It has been demonstrated on a CT-MR database of 35 prostate patients that using the
   * symmetrical resampling increase inverse-consistency by 1-2 orders of magnitude.
   * Citation to come...: D. Rivest-Henault et al. 2014.
   */
  switch(this->m_ResamplingType)
  {
  case EResamplingBasic: {
    std::cout << "# DRH - Mirorr - Using basic resampling" << std::endl;
    //Create an image resample filter to interpolate the FIXED image into the MOVING space
    resampler->SetInterpolator( m_Interpolator );
    resampler->SetTransform( m_Transform );
    resampler->SetDefaultPixelValue( 0 );
    resampler->SetUseReferenceImage(true);

    mask_resampler->SetInterpolator( m_MaskInterpolator );
    mask_resampler->SetTransform( m_Transform );
    mask_resampler->SetDefaultPixelValue( 0 );
    mask_resampler->SetUseReferenceImage(true);

    //Sets the convention of which image is to be resampled
    resampler->SetInput( m_FixedImage );
    resampler->SetReferenceImage(m_MovingImage);
    mask_resampler->SetInput( m_FixedMask );
    mask_resampler->SetReferenceImage( m_MovingMask );

    movingImage = m_MovingImage;
    movingMask  = m_MovingMask;
  } break;

  case EResamplingMiddle:
  case EResamplingMoving:
  case EResamplingFixed: {
    typename TMovingImage::PointType origin;
    typename TMovingImage::SpacingType spacing;
    typename TMovingImage::SizeType size;

    if (this->m_ResamplingType == EResamplingFixed) { //Use fixed image space
      std::cout << "# DRH - Mirorr - Using fixedImage resampling" << std::endl;
      origin  = this->m_FixedImage->GetOrigin();
      spacing = this->m_FixedImage->GetSpacing();
      size    = this->m_FixedImage->GetLargestPossibleRegion().GetSize();

      origin = m_Transform->GetInverseTransform()->TransformPoint(origin);
    } else if (this->m_ResamplingType == EResamplingMoving) {
      std::cout << "# DRH - Mirorr - Using movingImage resampling" << std::endl;
      origin  = this->m_MovingImage->GetOrigin();
      spacing = this->m_MovingImage->GetSpacing();
      size    = this->m_MovingImage->GetLargestPossibleRegion().GetSize();
    } else if (this->m_ResamplingType == EResamplingMiddle) {
      std::cout << "# DRH - Mirorr - Using middle space resampling" << std::endl;
      // Get the image centre coordinates
      PointType centerFixedPoint = this->ComputeImageCenter(this->m_FixedImage);
      PointType centerMovingPoint= this->ComputeImageCenter(this->m_MovingImage);

      //Compute the average (half-space) size, spacing, and origin
      typename TMovingImage::PointType::VectorType deltaF = this->m_FixedImage->GetOrigin() - centerFixedPoint;
      typename TMovingImage::PointType::VectorType deltaM = this->m_MovingImage->GetOrigin() - centerMovingPoint;

      origin.GetVnlVector() = ( centerMovingPoint + ((deltaF + deltaM) / 2.0) ).GetVnlVector();

      spacing = ( this->m_FixedImage->GetSpacing() + this->m_MovingImage->GetSpacing() ) / 2.0;

      size = ( this->m_FixedImage->GetLargestPossibleRegion().GetSize()
               + this->m_MovingImage->GetLargestPossibleRegion().GetSize() );
      size[0] = (size[0] + 0.5) / 2;
      size[1] = (size[1] + 0.5) / 2;
      size[2] = (size[2] + 0.5) / 2;
    } else {
      std::stringstream ss;
      ss
       << "ERROR: itk::MirorrRegistrationMethod::StartOptimization() "
       << " Unknown (sub-)resampling type. This should not happen..." << std::endl;
      throw std::runtime_error(ss.str());
    }

    //Create an image resample filter to interpolate the pseudo-FIXED image
    resampler->SetInterpolator(m_Interpolator);
    resampler->SetTransform(m_Transform);
    resampler->SetDefaultPixelValue(0);
    resampler->SetUseReferenceImage(false);
    resampler->SetInput(this->m_FixedImage);
    resampler->SetOutputOrigin(origin);
    resampler->SetOutputSpacing(spacing);
    resampler->SetOutputDirection(this->m_FixedImage->GetDirection());
    resampler->SetSize(size);

    //Create an image resample filter to interpolate the MOVING image
    typename ResampleFilterType::Pointer movingImageResampler = ResampleFilterType::New();
    movingImageResampler->SetInterpolator(m_Interpolator);
    movingImageResampler->SetDefaultPixelValue(0);
    movingImageResampler->SetUseReferenceImage(false);
    movingImageResampler->SetInput(this->m_MovingImage);
    movingImageResampler->SetOutputOrigin(origin);
    movingImageResampler->SetOutputSpacing(spacing);
    movingImageResampler->SetOutputDirection(this->m_MovingImage->GetDirection());
    movingImageResampler->SetSize(size);

    //typedef itk::ResampleImageFilter<MaskImageType, MaskImageType> MaskResampleFilterType;
    mask_resampler->SetInterpolator(m_MaskInterpolator);
    mask_resampler->SetTransform(m_Transform);
    mask_resampler->SetDefaultPixelValue(0);
    mask_resampler->SetUseReferenceImage(false);
    mask_resampler->SetInput(this->m_FixedMask);
    mask_resampler->SetOutputOrigin(origin);
    mask_resampler->SetOutputSpacing(spacing);
    mask_resampler->SetOutputDirection(this->m_FixedMask->GetDirection());
    mask_resampler->SetSize(size);

    typename MaskResampleFilterType::Pointer movingmask_resampler = MaskResampleFilterType::New();
    movingmask_resampler->SetInterpolator(m_MaskInterpolator);
    movingmask_resampler->SetDefaultPixelValue(0);
    movingmask_resampler->SetUseReferenceImage(false);
    movingmask_resampler->SetInput(this->m_MovingMask);
    movingmask_resampler->SetOutputOrigin(origin);
    movingmask_resampler->SetOutputSpacing(spacing);
    movingmask_resampler->SetOutputDirection(this->m_MovingMask->GetDirection());
    movingmask_resampler->SetSize(size);

    movingImageResampler->Update();
    movingmask_resampler->Update();

    movingImage = movingImageResampler->GetOutput();
    movingMask  = movingmask_resampler->GetOutput();
  } break;

  default:  {
    std::stringstream ss;
    ss
     << "ERROR: itk::MirorrRegistrationMethod::StartOptimization() "
     << " Unknown resampling type. This should not happen..." << std::endl;
    throw std::runtime_error(ss.str());
  } break;
  }

  this->SetupBlockMatcher();

  //Initial position
  m_LastTransformParameters = m_Transform->GetParameters();

  if( m_Verbosity >= 1 )
  {
    std::cout << "Initial Transform: " << m_LastTransformParameters << std::endl;
    m_BlockMatcher->displayInfo();
    std::cout << "N_Iterations: " << m_MaxIterations << std::endl;
  }

  boost::timer::cpu_timer iteration_timer;

  //For N iterationsm_MaxIterations
  const unsigned int maxIterations = static_cast<unsigned int>(std::abs(this->m_MaxIterations));
  for( unsigned int iteration=0; iteration<maxIterations; iteration++ )
  {
    iteration_timer.start();

    //1. Get interpolated fixed image for current transform
    resampler->Update();
    mask_resampler->Update();

    //2. Block Match the two images
    m_BlockMatcher->SetSearchImage( resampler->GetOutput() );
    m_BlockMatcher->SetBaseImage( movingImage );
    m_BlockMatcher->SetSearchMask( mask_resampler->GetOutput() );
    m_BlockMatcher->SetBaseMask( movingMask );

    try 
    {
      m_BlockMatcher->GetOffsets( m_InitialPositions, m_MovedPositions );
    } catch( std::exception &e ) {
      std::cerr << "Exception caught in BlockMatcher: " << e.what() << std::endl;
      exit(-1);
    }

    //Get the update to the transform
    OptimizeTransform();

    //4a. Check if parameters have changed
    double squareDifferenceInPosition = 0;
    for( unsigned int ii = 0; ii<m_LastTransformParameters.GetSize(); ++ii )
    {
      double diff = m_LastTransformParameters[ii] -
      m_Transform->GetParameters()[ii];
      squareDifferenceInPosition += diff*diff;
    }

    //Display
    if (this->m_Verbosity >= 1) {
      std::cout
      << std::fixed << std::setprecision(4)
      << "Iteration " << iteration + 1 << ",  " << "No. blocks: "
      << this->m_InitialPositions.size() << ",  ";
    }
    if (this->m_Verbosity >= 2) {
      std::cout << "\nPrev Pos: " << this->m_LastTransformParameters << ",\n";
    }
    if (this->m_Verbosity >= 1) {
      std::cout << "Opt Pos: " << std::setprecision(3)
      << PrettyPrint(this->m_Transform->GetParameters(), 7) << ",  "
      << "Delta: " << std::setw(9) << squareDifferenceInPosition << ",  "
      << std::setw(7) << boost::timer::format(iteration_timer.elapsed(), 2, "%ws") << std::endl;
    }

    //4. Update the transform
    m_LastTransformParameters = m_Transform->GetParameters();

    //4c. Terminate if movement less than 1/10 voxels or params hardly change
    if( squareDifferenceInPosition < 1e-3 )
    {
    	if( m_Verbosity >= 1 ) std::cout <<"Early breakout as movement small." << std::endl;
        break;
    }
  }
 }

/*
* PrintSelf
*/
template < typename TMovingImage, typename TFixedImage >
void
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::PrintSelf(std::ostream& os, Indent indent) const
 {
  Superclass::PrintSelf( os, indent );
  //  os << indent << "Metric: " << m_Metric.GetPointer() << std::endl;
  //  os << indent << "Optimizer: " << m_Optimizer.GetPointer() << std::endl;
  os << indent << "Transform: " << m_Transform.GetPointer() << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
  os << indent << "Moving Image: " << m_MovingImage.GetPointer() << std::endl;
  os << indent << "Fixed Image: " << m_FixedImage.GetPointer() << std::endl;
  os << indent << "Moving Image Region Defined: " << m_MovingImageRegionDefined << std::endl;
  os << indent << "Moving Image Region: " << m_MovingImageRegion << std::endl;
  //os << indent << "Initial Transform Parameters: " << m_InitialTransformParameters << std::endl;
  os << indent << "Last    Transform Parameters: " << m_LastTransformParameters << std::endl;
 }

// Generate Data
template < typename TMovingImage, typename TFixedImage >
void
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::GenerateData()
 {
  this->StartRegistration();
 }

// Get Output
template < typename TMovingImage, typename TFixedImage >
const typename MirorrRegistrationMethod<TMovingImage,TFixedImage>::TransformOutputType *
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::GetOutput() const
 {
  return static_cast< const TransformOutputType * >( this->ProcessObject::GetOutput(0) );
 }

template < typename TMovingImage, typename TFixedImage >
DataObject::Pointer
MirorrRegistrationMethod<TMovingImage,TFixedImage>
#if ITK_VERSION_MAJOR < 4
::MakeOutput(unsigned int idx)
#else
::MakeOutput(DataObjectPointerArraySizeType idx)
#endif
 {
  switch (idx)
  {
  case 0:
    return static_cast<DataObject*>(TransformOutputType::New().GetPointer());
    break;
  default:
    itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs");
    return 0;
  }
 }

template < typename TMovingImage, typename TFixedImage >
void
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::SetMovingImage( const MovingImageType * movingImage )
 {
  itkDebugMacro("setting Moving Image to " << movingImage );

  if (this->m_MovingImage.GetPointer() != movingImage )
  {
    this->m_MovingImage = movingImage;

    // Process object is not const-correct so the const_cast is required here
    this->ProcessObject::SetNthInput(0,
    const_cast< MovingImageType *>( movingImage ) );

    this->Modified();
  }
 }

template < typename TMovingImage, typename TFixedImage >
void
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::SetFixedImage( const FixedImageType * fixedImage )
 {
  itkDebugMacro("setting Fixed Image to " << fixedImage );

  if (this->m_FixedImage.GetPointer() != fixedImage )
  {
    this->m_FixedImage = fixedImage;

    // Process object is not const-correct so the const_cast is required here
    this->ProcessObject::SetNthInput(1,
    const_cast< FixedImageType *>( fixedImage ) );

    this->Modified();
  }
 }

template < typename TMovingImage, typename TFixedImage >
//typename MirorrRegistrationMethod<TMovingImage,TFixedImage>::ParametersType
void
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::OptimizeTransform( )
 {
  const unsigned int kDimension = 3;
  typedef typename itk::MatrixOffsetTransformBase<double,kDimension,kDimension> LocalTransformType;
  LocalTransformType::Pointer final_transform =
  dynamic_cast<LocalTransformType*>( m_Transform->CreateAnother().GetPointer() );

  PointType rotation_centre = m_Transform->GetCenter();
  final_transform->SetIdentity(); // This must be before SetCenter, otherwise center will be reset to 0,0,0
  final_transform->SetCenter( rotation_centre );
  ParametersType last_parameters( final_transform->GetParameters() );

  // LTS
  double fraction = 1.0; //This changes to 0.5 later - initially we use all points

  std::string transform_type = m_Transform->GetNameOfClass();

  // Affine transform gets performs a single pass
  //While change in transform is sufficiently low
  double last_delta_param_tol = 0;
  double delta_param_tol = 1;
  int iter = -1;
  //Check tolerance met, max no. iterations, flip flopping between solutions
  while( delta_param_tol > 1e-3 && iter < 100
      && fabs(delta_param_tol-last_delta_param_tol) > 1e-3 ) //Flip flopping
  {
    ++iter;

    if( m_InitialPositions.empty() )
    {
      break;
    }

    //Get residuals
    std::vector<double> residuals(m_InitialPositions.size(), 0.0);
    for( unsigned int point = 0; point<m_InitialPositions.size(); ++point)
    {
      PointType transformed_point = final_transform->TransformPoint( m_InitialPositions[point] );
      for( unsigned int col =0; col<kDimension; ++col )
      {
        //Transform m_InitialPositions
        double dist = transformed_point[col] - m_MovedPositions[point][col];
        residuals[point] += dist * dist;
      }
    }


    //Find nth residual
    double nth_residual = 0.0;
    if (fraction > 0.99999) {
      //David RH, 2012-11-09: 'if' necessary since std::nth_element does
      //nothing when fraction ~= 1.0.
      nth_residual = *std::max_element(residuals.begin(), residuals.end());
    } else {
      //Sort points in order of residuals and retain only the 'fraction' (~50%) of points
      //with the lowest residual
      std::vector<double> residuals_sorted( residuals );

      //if( residuals_sorted.size() > 0)
      {
        unsigned int len =
        static_cast<unsigned int>( fraction * residuals_sorted.size() );
        len = std::min<unsigned int>( len, residuals_sorted.size()-1 );
        std::nth_element( residuals_sorted.begin(),
        residuals_sorted.begin()+len, residuals_sorted.end() );
        nth_residual = residuals_sorted[len];
      }
    }

    //Extract points with lowest residual
    PointListType initialSubset; initialSubset.reserve(m_InitialPositions.size());
    PointListType movedSubset;   movedSubset.reserve(m_InitialPositions.size());
    for( unsigned int point = 0; point<m_InitialPositions.size(); ++point)
      if(residuals[point] <= nth_residual)
      {
        initialSubset.push_back( m_InitialPositions[point] );
        movedSubset.push_back( m_MovedPositions[point] );
      }

    //After first iteration, we reduce the fraction to 0.5
    fraction = 0.5;

    //Begin LS ************************
    //Get barycentres
    vnl_vector<double> initialMean = GetBaryCentre( initialSubset );
    vnl_vector<double> movedMean = GetBaryCentre( movedSubset );

    //Convert subset to barycentric coordinates in mm
    for( unsigned int point = 0; point<initialSubset.size(); ++point)
      for( unsigned int axis=0; axis<kDimension; ++axis )
      {
        initialSubset[point][axis] = (initialSubset[point][axis] - initialMean[axis]);// * spacing[axis];
        movedSubset[point][axis] = (movedSubset[point][axis] - movedMean[axis]);// * spacing[axis];
      }

    //Branch for Affine / rigid *****************************
    //Only a parameters typ object should leave this part
    if( m_Verbosity >= 3 )
    {
      std::cout << "\n  Inner iteration: " << iter << std::endl;
      std::cout << std::setprecision(5) << "  Init centre: " << initialMean << "   Moved centre: " << movedMean << std::endl;
      std::cout << "  Num points: " << initialSubset.size()
          << " of " << m_InitialPositions.size()
          << "  Nth residual: " << nth_residual << std::endl;
      std::cout << "  TransformType: " << transform_type << std::endl;
    }

    ParametersType new_parameters( last_parameters );

    if( 0 == transform_type.compare("AffineTransform") )
    {
      //Get summed outer product of x and u
      vnl_matrix<double> outer_xy = CovarianceOfPoints(movedSubset,initialSubset); //ORder is important here
      vnl_matrix<double> outer_xx = CovarianceOfPoints(initialSubset,initialSubset);

      //Get rotation matrix and translation
      vnl_matrix<double> rotation = outer_xy;
      if( fabs(vnl_determinant(outer_xx))<1e-4 ) //Check for low determinants
      std::cout << "\nERROR: Problem inverting outer_xx as its determinant is low:\n" << outer_xx
      << "Treating it as the identity matrix.\n" << std::endl;
      else
      rotation = rotation * vnl_inverse( outer_xx );

      vnl_vector<double> movedMeanFromCentre   = movedMean - rotation_centre.GetVnlVector();
      vnl_vector<double> initialMeanFromCentre = initialMean - rotation_centre.GetVnlVector();
      vnl_vector<double> translation = movedMeanFromCentre - rotation * initialMeanFromCentre;

      if( m_Verbosity >= 3 )
      {
        std::cout<<"Field"<<std::endl;
        for(unsigned int ii=0; ii<movedSubset.size(); ++ii)
        {
          std::cout
          << movedSubset[ii][0] << " "
          << movedSubset[ii][1] << " "
          << movedSubset[ii][2] << " "
          << initialSubset[ii][0] << " "
          << initialSubset[ii][1] << " "
          << initialSubset[ii][2] << " "
          <<std::endl;
        }
        std::cout<<"----"<<std::endl;

        std::cout << "outer_xy:\n" <<  outer_xy;
        std::cout << "outer_xx:\n" << outer_xx;
        std::cout << "Rotation:\n" << rotation;
        std::cout << "Translation: " << translation << std::endl << std::endl;
      }

      //TODO - check this assignment is correct
      //Extract parameters
      //Rotation
      new_parameters[0] = rotation[0][0]; new_parameters[1] = rotation[0][1]; new_parameters[2] = rotation[0][2];
      new_parameters[3] = rotation[1][0]; new_parameters[4] = rotation[1][1]; new_parameters[5] = rotation[1][2];
      new_parameters[6] = rotation[2][0]; new_parameters[7] = rotation[2][1]; new_parameters[8] = rotation[2][2];
      //Translation
      new_parameters[9] = translation[0]; new_parameters[10]= translation[1]; new_parameters[11]= translation[2];

      //Get summed change in parameters
      last_delta_param_tol = delta_param_tol;
      delta_param_tol = GetParametersDelta( last_parameters, new_parameters);
    }
    else if( 0 == transform_type.compare("Euler3DTransform") )
    {
      //Create quaternion matrix for each point pair and sum Q^T . Q
      vnl_quaternion<double> quaternion = CalculateRigidQuaternion( initialSubset, movedSubset );

      //Get rotation matrix and translation
      vnl_matrix<double> rotation = quaternion.rotation_matrix_transpose().transpose();
      vnl_vector<double> movedMeanFromCentre   = movedMean - rotation_centre.GetVnlVector();
      vnl_vector<double> initialMeanFromCentre = initialMean - rotation_centre.GetVnlVector();
      vnl_vector<double> translation = movedMeanFromCentre - rotation * initialMeanFromCentre;
      //vnl_vector<double> translation = rotation * ( movedMean - initialMean );
      vnl_vector<double> euler_angles = quaternion.rotation_euler_angles();

      //Modulus of Euler angles
      if( m_Verbosity >= 3 )
      {
        //std::cout << "  Rotation:\n" << rotation;
        std::cout << "  Translation: " << translation << std::endl;
        std::cout << "  Euler Angles: " << euler_angles << std::endl;
      }

      //Extract parameters
      for( int ii=0; ii<3; ++ii ) new_parameters[ii] = euler_angles[ii];
      for( int ii=0; ii<3; ++ii ) new_parameters[ii+3] = translation[ii];
      //for( int ii=0; ii<3; ++ii ) rotation_centre[ii] = initialMean[ii];

      //Get summed change in parameters
      last_delta_param_tol = delta_param_tol;
      delta_param_tol = GetParametersDelta(new_parameters,last_parameters);
    }
    else if( 0 == transform_type.compare("QuaternionRigid") )
    {
      //Create quaternion matrix for each point pair and sum Q^T . Q
      vnl_quaternion<double> quaternion = CalculateRigidQuaternion( initialSubset, movedSubset );

      //Get rotation matrix and translation
      //TODO: Test
      vnl_matrix<double> rotation = quaternion.rotation_matrix_transpose().transpose();
      vnl_vector<double> movedMeanFromCentre   = movedMean - rotation_centre.GetVnlVector();
      vnl_vector<double> initialMeanFromCentre = initialMean - rotation_centre.GetVnlVector();
      vnl_vector<double> translation = movedMeanFromCentre - rotation * initialMeanFromCentre;

      //Modulus of Euler angles
      if( m_Verbosity >= 3 )
      {
        std::cout << "Rotation:\n" << rotation;
        std::cout << "Translation: " << translation << std::endl;
      }

      //Extract parameters
      //Rotation
      for( int ii=0; ii<4; ++ii ) new_parameters[ii] = quaternion[ii];
      for( int ii=4; ii<7; ++ii ) new_parameters[ii] = translation[ii-4];

      //Get summed change in parameters
      last_delta_param_tol = delta_param_tol;
      delta_param_tol = GetParametersDelta(last_parameters,new_parameters);
    }// If Branch for Affine / rigid / QuaternionRigid*****************************

    // End LS
    // LTS ****************************************************

    //Display what is happening
    if( m_Verbosity >= 3 )
    {
      std::cout << "  Prev Parameters: " <<  last_parameters;
      std::cout << "\n  New  Parameters: " <<  new_parameters;
      std::cout << "   Shift: " << delta_param_tol << std::endl;
    }

    //Update by storing wap - do NOT compose as we do not shift the moving points
    //only select a different subset
    rotation_centre = initialMean.begin();
    final_transform->SetCenter( rotation_centre );
    final_transform->SetParameters( new_parameters );
  }

  //Now compose - Remember we need to pre multiply
  m_Transform->Compose( final_transform, true );
 }

template < typename TMovingImage, typename TFixedImage >
vnl_vector<double>
MirorrRegistrationMethod<TMovingImage,TFixedImage>::
GetBaryCentre( const PointListType &points ) const
{
  vnl_vector<double> mean(3,0.0);
  for( unsigned int point = 0; point<points.size(); ++point)
  for( unsigned int element =0; element<3; ++element )
  mean[element] += points[point][element];

  //Divide results to get mean
  mean /= points.size();

  return mean;
}

template < typename TMovingImage, typename TFixedImage >
double
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::GetParametersDelta(
const ParametersType & params1,
const ParametersType & params2
) const
{
  double delta_param = 0;
  for( unsigned int ii=0; ii<params1.GetNumberOfElements(); ++ii )
  {
    double difference = params1[ii] - params2[ii];
    delta_param += difference * difference;
  }
  return delta_param;
}

template < typename TMovingImage, typename TFixedImage >
vnl_matrix<double>
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::CovarianceOfPoints(
const PointListType &initialSubset,
const PointListType &movedSubset
) const
{
  unsigned int kDimension = TMovingImage::ImageDimension;

  //Create matrix
  vnl_matrix<double> covariance(kDimension,kDimension,0.0);

  //Calculate covariance
  for( unsigned int point = 0; point<initialSubset.size(); ++point)
  for( unsigned int row =0; row<kDimension; ++row )
  for( unsigned int col=0; col<kDimension; ++col )
  covariance[row][col] += initialSubset[point][row] * movedSubset[point][col];

  //Normalise matrix values
  covariance /= initialSubset.size();

  return covariance;
}

template < typename TMovingImage, typename TFixedImage >
vnl_quaternion<double>
MirorrRegistrationMethod<TMovingImage,TFixedImage>
::CalculateRigidQuaternion(
const PointListType &initialSubset,
const PointListType &movedSubset
) const
{
  //Create quaternion matrix for each point pair and sum Q^T . Q
  vnl_matrix<double> total_quat_matrix(4,4,0.0);

  for( unsigned int point = 0; point<initialSubset.size(); ++point)
  {
    //Get upper triangle of matrix
    vnl_matrix<double> quat_matrix(4,4,0.0);
    for( unsigned int ii = 0; ii<3; ++ii )
      quat_matrix[0][1+ii] = initialSubset[point][ii] - movedSubset[point][ii];
    quat_matrix[1][2] = -initialSubset[point][2] - movedSubset[point][2]; //z
    quat_matrix[1][3] =  initialSubset[point][1] + movedSubset[point][1]; //y
    quat_matrix[2][3] = -initialSubset[point][0] - movedSubset[point][0]; //x

    //Copy negative to lower triangle
    for( unsigned int ii = 0; ii<4; ++ii )
    for( unsigned int jj = ii+1; jj<4; ++jj )
    quat_matrix[jj][ii] = -quat_matrix[ii][jj];

    //Add Q^T Q to total
    total_quat_matrix += quat_matrix.transpose() * quat_matrix;
  }

  //Get primary eigenvector of matrix for use as a quaternion
  vnl_symmetric_eigensystem<double> quat_eigs( total_quat_matrix );
  vnl_vector<double> evec = quat_eigs.get_eigenvector(0);
  vnl_quaternion<double> quaternion(
  (evec[1]),(evec[2]),(evec[3]),(evec[0]) );

  //Modulus of Euler angles
  if( m_Verbosity >= 3 )
  {
    std::cout << "  Total quat matrix:\n" <<  total_quat_matrix;
    std::cout << "  Eigvec 0: " << quaternion << std::endl;
  }
  return quaternion;
}

} // end namespace itk
#endif
