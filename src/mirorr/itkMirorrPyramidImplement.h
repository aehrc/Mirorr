/*=========================================================================
Program: mirorr
Module: itkMirorrPyramidImplement.h
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
An implementation of the Mirorr registration class using multiple scales
*/

#ifndef __itkMirorrPyramidImplement_h
#define __itkMirorrPyramidImplement_h

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include <itkLinearInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkMultiResolutionImageRegistrationMethod.h>

#include <itkPowellOptimizer.h>
#include <itkImageRegistrationMethod.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkGradientDifferenceImageToImageMetric.h>
#include <itkMutualInformationImageToImageMetric.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>

#if ITK_VERSION_MAJOR < 4
#include <itkOrientedImage.h>
#else
#include <itkImageMaskSpatialObject.h>
#include <itkImage.h>
#endif

#include "itkMirorrRegistrationMethod.h"
#include "itkSymMirorrRegistrationMethod.h"
#include "itkPyramidScheduleTuner.h"

/** Implements the Mirorr registration algorithm in a
 *  a multi-scale scheme, taking away some design decisions:
 *  Interpolation: Linear
 *  Image Type: OrientedImage<short, ?>
 *  Only the image dimension and the transformation type need to be specified.
 *  Some useful helper functions are also implemented:
 *  - Calculate image moments to automatically select initial transform
 *  - Output results at key parts of registration
 *
 * Used by: MirorrPyramidWrapper
 * Uses: itkMirorrRegistrationMethod
 */
//template <unsigned int DIMENSION >
class MirorrPyramidImplement
{
public:
  static const unsigned int DIMENSION = 3;
  //! Standard class typedefs. */
  typedef float                                        PixelType;
#if ITK_VERSION_MAJOR < 4
  typedef itk::OrientedImage<PixelType, DIMENSION>     ImageType;
  typedef itk::OrientedImage<unsigned char, DIMENSION> MaskType;
#else
  typedef itk::Image<PixelType, DIMENSION>     ImageType;
  typedef itk::Image<unsigned char, DIMENSION> MaskType;
  typedef itk::ImageMaskSpatialObject<DIMENSION>       MaskObjectType;
#endif
  
  typedef /*typename*/ ImageType::Pointer              ImagePointer;
  typedef MaskType::Pointer                            MaskPointer;

  //! Typedefs defining the particular registration algorithm */
  typedef itk::LinearInterpolateImageFunction< ImageType, double > InterpolatorType;
  typedef itk::LinearInterpolateImageFunction< MaskType, double > MaskInterpolatorType;
  typedef itk::MirorrRegistrationMethod< ImageType, ImageType >   RegistrationType;
  typedef itk::SymMirorrRegistrationMethod< ImageType, ImageType >   SymRegistrationType;

  typedef /*typename*/ RegistrationType::ParametersType            ParametersType;
  typedef /*typename*/ RegistrationType::TransformType             TransformType;
  typedef /*typename*/ itk::PyramidScheduleTuner<DIMENSION>        PyramidScheduleTunerType;

  //! We also use ITK for a final fine tuning registration */
  typedef itk::PowellOptimizer                                                   SecondOptimizerType;
  //typedef itk::MattesMutualInformationImageToImageMetric< ImageType, ImageType > SecondMetricType;
  typedef itk::ImageToImageMetric< ImageType, ImageType >                        SecondMetricType;
  typedef itk::MattesMutualInformationImageToImageMetric< ImageType, ImageType > SecondMetricTypeMMI;
  typedef itk::NormalizedCorrelationImageToImageMetric< ImageType, ImageType >   SecondMetricTypeNCC;
  typedef itk::GradientDifferenceImageToImageMetric< ImageType, ImageType >      SecondMetricTypeGD;
  typedef itk::MutualInformationImageToImageMetric< ImageType, ImageType >       SecondMetricTypeMI;
  typedef itk::ImageRegistrationMethod< ImageType, ImageType >                   SecondRegistrationType;

  //! Define a multi-resolution image filter*/
  typedef itk::RecursiveMultiResolutionPyramidImageFilter<ImageType, ImageType > ImagePyramidFilterType;

  enum ERegistrationType {ERegClassic=0, ERegSymmetrical};
  enum ESecondRegistrationType {ERegItkMMI=0, ERegItkNCC, ERegItkGD, ERegItkMI}; //MI requires normalized images, MMI doesn't
  enum EResamplingType {EResamplingBasic=0, EResamplingMiddle, EResamplingMaxResolution, EResamplingMaxSize, EResamplingFixed, EResamplingMoving}; //Mirorr inner resampling mode

  //! Constructor. A fixed and moving image must be defined along with additional inputs
  MirorrPyramidImplement();

  //! Run the hierarchical registration.
  void RunMirorr( /*typename*/ TransformType::Pointer iTransform );

  void SetRegistrationType(ERegistrationType rtype);
  void SetRegistrationTypeClassic();
  void SetRegistrationTypeSymmetrical();
  ERegistrationType GetRegistrationType();

  void SetResamplingType(EResamplingType rtype);
  EResamplingType GetResamplingType();

  /** Set and Get methods */
  void SetVerbosity( int in ) { verbosity = in; registration->SetVerbosity(  verbosity  ); }
  void SetMaxIterations( int in )
  {
    if( in != 0 ) registration->SetMaxIterations( std::abs( in ) );
    halveIterations = false;
    if( in < 0 )
      halveIterations = true;
  }
  void SetFixedImage( ImagePointer in ) { fixedImage = in; }
  void SetMovingImage( ImagePointer in ) { movingImage = in; }
  void SetFixedMask( MaskPointer in ) { fixedMask = in; }
  void SetMovingMask( MaskPointer in ) { movingMask = in; }
  void SetLevelToChangeMethod( int in ) { m_LevelToChangeMethod = in; }
  void SetUseBlockMatchingAlgorithm( bool in ) { m_UseBlockMatchingAlgorithm = in; }
  void SetUseMutualInformationAlgorithm( bool in ) { m_UseMutualInformationAlgorithm = in; }

  RegistrationType::Pointer GetRegistrationObject() { return registration; }

  //! Access the methods to tune the pyramid schedule
  PyramidScheduleTunerType::Pointer GetPyramidScheduleTuner() { return m_PyramidScheduleTuner; }

  //! Resample the fixed image into another space defined by the input transform
  ImagePointer GetResampledImage( /*typename*/ TransformType::Pointer transform,
      bool resampling_moving_image = false );

  //! Set the metric for the secondary ITK registration
  void SetSecondaryRegistrationMetricToMMI(); //Mattes mutual information
  void SetSecondaryRegistrationMetricToNCC(); //Normalized cross-correlation
  void SetSecondaryRegistrationMetricToGD();  //Gradient difference
  void SetSecondaryRegistrationMetricToMI();  //Viola and Well mutual information (requires normalized input images)

  //! Request a particular sample rate
  void SetRequestedSampleRate( double in )
  { m_RequestedSampleRate = std::max(0.0,std::min(1.0,in)); }
  double GetRequestedSampleRate()
  { return m_RequestedSampleRate; }

private:
  //! The moving and fixed images */
  ImagePointer movingImage;
  ImagePointer fixedImage;
  MaskPointer movingMask;
  MaskPointer fixedMask;

  ERegistrationType registrationType;
  ESecondRegistrationType secondRegistrationType;
  EResamplingType resamplingType;

  //! Define how much is displayed */
  int verbosity;
  //! Should we halve the number of iterations each time we go up a pyramid level
  bool halveIterations;
  //! The minimum length of the image in any one dimension
  const double minLength;
  //! The pyramid level when we switch from aladin registration to ITK
  int m_LevelToChangeMethod;
  bool m_UseBlockMatchingAlgorithm;
  bool m_UseMutualInformationAlgorithm;

  //! Pointers to classes that do actual registration
  /*typename*/ InterpolatorType::Pointer         interpolator;
  /*typename*/ MaskInterpolatorType::Pointer     mask_interpolator;
  /*typename*/ RegistrationType::Pointer         registration;

  //------------------------------------------------------
    //DRH: Test SymMirorr
  /*typename*/ InterpolatorType::Pointer         interpolatorFixed;
  /*typename*/ MaskInterpolatorType::Pointer     fixedmask_interpolator;
  //------------------------------------------------------

  //! Image pyramids
  /*typename*/ ImagePyramidFilterType::Pointer   movingPyramid;
  /*typename*/ ImagePyramidFilterType::Pointer   fixedPyramid;

  //! Object for deciding what schedule to use
  PyramidScheduleTunerType::Pointer m_PyramidScheduleTuner;

  //! Extra objects used for final ITK registration
  SecondOptimizerType::Pointer      m_SecondOptimizer;
  SecondMetricType::Pointer         m_SecondMetric;
  SecondRegistrationType::Pointer   m_SecondRegistration;
  double m_RequestedSampleRate;

  void InitializeRegistrationObjet();

  void SetupRegistrationObjetResamplingMode();

  //! Create the image-to-image metric object for the secondary ITK registration
  void CreateSecondaryRegistrationMetric();

  //! Helper function to generate the image pyramids
  void CreateImagePyramids();

  //! Display useful image information for debugging */
  void displayImgStats( std::string prefix,
      /*typename*/ ImageType::Pointer image,
      bool displayDir = false )
  {
    std::cout << prefix
    << " Size: " << image->GetLargestPossibleRegion().GetSize()
    << " Origin: " << image->GetOrigin()
    << " Spacing: " << image->GetSpacing();

    if( displayDir )
      std::cout << " Direction: " << image->GetDirection();
    else
      std::cout << std::endl;
  }
  void DisplayCurrentLevel( int level );
  void DisplayLevelTimings( int level, boost::timer::cpu_timer &level_timer, boost::timer::cpu_timer &total_timer );

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMirorrPyramidImplement.txx"
#endif


#endif
