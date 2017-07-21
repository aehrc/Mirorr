/*=========================================================================
Program: mirorr
Module: itkMirorrPyramidImplement.txx
Author: Nicholas Dowson
Created: 20 Feb 2009

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
//History
// 24 Apr 2009: Moving small error in the way the parameters are updated

#ifndef __itkMirorrPyramidImplement_txx
#define __itkMirorrPyramidImplement_txx

#include "itkMirorrPyramidImplement.h"
#include <itkCommand.h>
#include <itkResampleImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkWindowedSincInterpolateImageFunction.h>
#include <itkImageFileWriter.h>
#include <boost/timer/timer.hpp>
#include <string>
#include <iomanip>
#include <exception>
#include "itkIOUtils.h"
//#include "../../../../../../usr/local/cuda-7.0/targets/x86_64-linux/include/nppdefs.h"

#include <itkChangeInformationImageFilter.h>


MirorrPyramidImplement
::MirorrPyramidImplement()
:
minLength( 32 ) //This ensures that minimum length along
//one side of image is not smaller than 32 voxels
{
  interpolator = InterpolatorType::New();
  mask_interpolator = MaskInterpolatorType::New();

  // DRH: this is for symmetrical registration:
  fixedmask_interpolator = MaskInterpolatorType::New();
  interpolatorFixed = InterpolatorType::New();

  verbosity          = 1;

  this->registrationType = ERegSymmetrical;
  this->InitializeRegistrationObjet(); //Initialise this->registration
  this->resamplingType = EResamplingMiddle;

  halveIterations = false;

  m_PyramidScheduleTuner = PyramidScheduleTunerType::New();

#ifdef USE_NPW
  registration->SetNPWbins( 32 ); //inputs.NPWbins );
#ifndef USE_OPENCL
  registration->SetNPWscaleFactor( 2 ); //inputs.NPWscaleFactor );
  registration->SetNPWshapeExpansion( false ); //inputs.NPWshapeExpansion );
#endif
#endif

}

void
MirorrPyramidImplement
::RunMirorr( /*typename*/ TransformType::Pointer iTransform )
{
  //Check moving and fixed image have been assigned
  if( movingImage.IsNull() )
    throw std::runtime_error( "RunMirorr Exception: Moving Image has not been assigned." );
  //itkExceptionMacro( "RunMirorr Exception: Moving Image has not been assigned." );
  if( fixedImage.IsNull() )
    throw std::runtime_error( "RunMirorr Exception: Fixed Image has not been assigned." );
  //itkExceptionMacro( "RunMirorr Exception: Fixed Image has not been assigned." );

  //Create the image pyramid
  CreateImagePyramids();

  SetupRegistrationObjetResamplingMode();

  //Initialize with identity transform. How can we set this properly?
  registration->SetTransform( iTransform );
  ParametersType currentParameters( iTransform->GetParameters() );

  //Check change-over level is >= min and <= max
  unsigned int t_min = m_PyramidScheduleTuner->GetCroppedLevel(m_PyramidScheduleTuner->GetLevelMin());
  unsigned int t_max = m_PyramidScheduleTuner->GetCroppedLevel(m_PyramidScheduleTuner->GetLevelMax(),movingPyramid->GetNumberOfLevels());

  boost::timer::cpu_timer level_timer, total_timer;

  //Now run the registration algorithm
  int first_max_iterations = registration->GetMaxIterations();
  double portionMatchesKept = registration->GetPortionMatchesKept();

  if( verbosity >= 1 )
    std::cout << "\nBlock Matching registration start ====== \n" << std::endl;
  int max_iterations = registration->GetMaxIterations();

  //EXCLUSIVE pyramid level
  for( unsigned int level=0; level< t_max; ++level )
  {
    //Get the current images
    registration->SetMovingImage( movingPyramid->GetOutput(level) );
    registration->SetFixedImage( fixedPyramid->GetOutput(level) );
    registration->SetMovingMask(movingMask);
    registration->SetFixedMask(fixedMask);

    registration->SetMovingImageRegion( movingImage->GetLargestPossibleRegion() );
    if( level <= 1 ) //For small images, use all available data
      registration->SetPortionMatchesKept( 1.0 );
    else
      registration->SetPortionMatchesKept( portionMatchesKept );

    //Set the parameters
    registration->GetTransform()->SetParameters( currentParameters );

    //Set the number of iterations
    if( max_iterations != 0 )
      registration->SetMaxIterations( max_iterations );
    if( halveIterations )
      max_iterations = std::max( 1, max_iterations/2 );

    //Tell user what is going on
    if( verbosity >= 1 )
      DisplayCurrentLevel( level + t_min -1);

    //Run registration
    registration->Update();

    //Update the current parameters
    //ParametersType lastParameters(currentParameters);
    currentParameters = registration->GetLastTransformParameters();

    if( verbosity >= 1 ) {
      DisplayLevelTimings( level, level_timer, total_timer );
      std::cout << std::endl;
    }
  }

  if( verbosity >= 1)
    std::cout << "Block Match Registration DONE" << std::endl;

  iTransform->SetParameters( currentParameters );

  //Reset iterations
  registration->SetMaxIterations( first_max_iterations );

  //Extract the optimized parameters
  if(verbosity >= 1)
    std::cout << "Best Pos: " << std::setprecision(5) << itk::PrettyPrint(iTransform->GetParameters()) << std::endl;
}

void
MirorrPyramidImplement
::SetRegistrationType(ERegistrationType rtype)
{
  if (this->registrationType != rtype) {
    this->registrationType = rtype;
    this->InitializeRegistrationObjet();
  }
}

void
MirorrPyramidImplement
::SetRegistrationTypeClassic()
{
  if (this->registrationType != ERegClassic) {
      this->registrationType = ERegClassic;
      this->InitializeRegistrationObjet();
    }
}

void
MirorrPyramidImplement
::SetRegistrationTypeSymmetrical()
{
  if (this->registrationType != ERegSymmetrical) {
    this->registrationType = ERegSymmetrical;
    this->InitializeRegistrationObjet();
  }
}

MirorrPyramidImplement::ERegistrationType
MirorrPyramidImplement
::GetRegistrationType()
{
  return this->registrationType;
}

void
MirorrPyramidImplement
::SetResamplingType(EResamplingType rtype)
{
  this->resamplingType = rtype;
}

MirorrPyramidImplement::EResamplingType
MirorrPyramidImplement
::GetResamplingType()
{
  return this->resamplingType;
}


//template <unsigned int DIMENSION>
void
MirorrPyramidImplement/*<DIMENSION>*/
::DisplayCurrentLevel( int level )
{
  std::cout << "Pyramid level: " << level+1 << std::endl;

  if( verbosity >= 2 ) {
    //Moving image dimensions
    displayImgStats("Moving ", movingPyramid->GetOutput(level) );
    displayImgStats("Fixed ", fixedPyramid->GetOutput(level) );

    //Subsampling
    std::cout << "Moving Shrink Factor: ";
    for( unsigned int ii=0; ii<DIMENSION; ++ii )
      std::cout << movingPyramid->GetSchedule()[level][ii] << " ";
    std::cout << "  Fixed Shrink Factor: ";
    for( unsigned int ii=0; ii<fixedImage->GetImageDimension(); ++ii )
      std::cout << fixedPyramid->GetSchedule()[level][ii] << " ";
    std::cout << std::endl;
  }
}

//template <unsigned int DIMENSION>
void
MirorrPyramidImplement/*<DIMENSION>*/
::DisplayLevelTimings( int level, boost::timer::cpu_timer &level_timer, boost::timer::cpu_timer &total_timer )
{
  std::cout << "Level " << level+1 << " completed in "
  << std::setprecision(2) << level_timer.elapsed().wall/1.E9 << "s of "
  << std::setprecision(2) << total_timer.elapsed().wall/1.E9
  << "s <-^------"
  << std::endl;
  level_timer.start();
}

/*typename*/ MirorrPyramidImplement/*<DIMENSION>*/::ImagePointer
MirorrPyramidImplement/*<DIMENSION>*/
::GetReorientedImage(
    /*typename*/ TransformType::Pointer itransform,
    //bool resample_moving
                 ImageType::Pointer image
)
{
  typedef itk::ChangeInformationImageFilter< ImageType > ChangeInfoType;
  TransformType::Pointer transform = itransform;

  ChangeInfoType::Pointer imageChanger = ChangeInfoType::New();
  //ImageType::Pointer image;
  //
  //  if( resample_moving )
  //  {
  //    image = movingImage;
  //    itransform = dynamic_cast<TransformType*>(
  //        transform->GetInverseTransform().GetPointer() );
  //  }
  //  else
  //    image = fixedImage;
  imageChanger->SetInput( image );
  imageChanger->ChangeDirectionOn();
  imageChanger->ChangeOriginOn();

  //Get relevant information
  typedef ImageType::DirectionType DirectionType;
  typedef ImageType::PointType PointType;
  typedef ImageType::SpacingType SpacingType;
  DirectionType in_direction = image->GetDirection();
  PointType in_origin = image->GetOrigin();
  //SpacingType in_spacing = image->GetSpacing();

  DirectionType out_direction;
  PointType out_origin;
  SpacingType out_spacing;

  //Get the matrix manually
  {
    for(int ii=0;ii<3; ++ii)
    {
      itk::Point<double, 3> ipts_4d, ipts_4d_zero;
      itk::Point<double, 3> ipts_3d, ipts_3d_zero, ipts_3d_diff;
      itk::Point<double, 3> opts_3d, opts_3d_zero, opts_3d_diff;

      itk::ContinuousIndex<double,3> index;
      index.Fill(0);
      image->TransformContinuousIndexToPhysicalPoint( index, ipts_4d_zero );
      index[ii] = 1;
      image->TransformContinuousIndexToPhysicalPoint( index, ipts_4d );

      for(int jj=0; jj<3; ++jj ) ipts_3d[jj] = ipts_4d[jj];
      for(int jj=0; jj<3; ++jj ) ipts_3d_zero[jj] = ipts_4d_zero[jj];

      double mag_3d_diff = 0;
      for(int jj=0; jj<3; ++jj )
      {
        ipts_3d_diff[jj] = ipts_4d[jj]-ipts_4d_zero[jj];
        mag_3d_diff += ipts_3d_diff[jj]*ipts_3d_diff[jj];
      }
      //Normalise
      mag_3d_diff = sqrt(mag_3d_diff);
      for(int jj=0; jj<3; ++jj )
        ipts_3d_diff[jj] /= mag_3d_diff;

      //std::cout<<"Point "<<ii<<": "<<ipts_3d_diff<<" ---> ";

      opts_3d = transform->TransformPoint( ipts_3d );
      opts_3d_zero = transform->TransformPoint( ipts_3d_zero );
      opts_3d_diff = transform->TransformPoint( ipts_3d_diff );

      mag_3d_diff = 0;
      for(int jj=0; jj<3; ++jj )
      {
        opts_3d_diff[jj] = opts_3d[jj] - opts_3d_zero[jj];
        mag_3d_diff += opts_3d_diff[jj]*opts_3d_diff[jj];
      }
      //Normalise!
      mag_3d_diff = sqrt(mag_3d_diff);
      for(int jj=0; jj<3; ++jj )
        opts_3d_diff[jj] /= mag_3d_diff;

      //std::cout<<opts_3d_diff<<std::endl;

      for(int jj=0; jj<3; ++jj )
        out_direction[jj][ii] = opts_3d_diff[jj];

    }
    //out_direction[3][3] = 1;

    //BEWARE of negative determinants - they will flip your axes!
    double det = fabs(vnl_det( out_direction.GetVnlMatrix() ));
    for( int ii=0; ii<3; ++ii )
      for( int jj=0; jj<3; ++jj )
        out_direction[ii][jj] /= det;
  }

  itk::Point<double, 3> in_origin3d;
  for( int ii=0; ii<3; ++ii ) in_origin3d[ii] = in_origin[ii];
  itk::Point<double, 3> out_origin3d = transform->TransformPoint( in_origin3d );
  for( int ii=0; ii<3; ++ii ) out_origin[ii] = out_origin3d[ii];

//  std::cout << "   Tfm Offset: " << transform->GetOffset() << std::endl;
//  std::cout << "   Tfm Centre: " << transform->GetCenter() << std::endl;
//  std::cout << "Tfm Translate: " << transform->GetTranslation() << std::endl;
//  std::cout << "   Tfm Rotate: " << transform->GetMatrix() << std::endl;

  //Apply the transform
  imageChanger->SetOutputDirection( out_direction );
  imageChanger->SetOutputOrigin( out_origin );
  imageChanger->Update();

  return imageChanger->GetOutput();
}


itk::InterpolateImageFunction<MirorrPyramidImplement::ImageType>::Pointer
MirorrPyramidImplement::GetInterpolatorFromString(std::string interpolator_name)
{
  /**
   * The interpolators provided by this function are meant to provide the "best" (YMMV) interpolation
   * method for a given neighborhood size (thus related to computational cost). Based on:
   *
   * Erik H. W. Meijering, Wiro J. Niessen, Josien P. W. Pluim, Max A. Viergever: Quantitative Comparison
   * of Sinc-Approximating Kernels for Medical Image Interpolation. MICCAI 1999, pp. 210-217
   */
  enum EInterpolatorType {NN, LINEAR, BSPLINE, SINC, _ERROR_};
  EInterpolatorType eInterpolatorType = _ERROR_;

  if (interpolator_name.find("nn") != std::string::npos) {
    eInterpolatorType = NN;
  } else if (interpolator_name.find("linear") != std::string::npos) {
    eInterpolatorType = LINEAR;
  } else if(interpolator_name.find("bspline") != std::string::npos) {
    eInterpolatorType = BSPLINE;
  } else if(interpolator_name.find("sinc") != std::string::npos) {
    eInterpolatorType = SINC;
  } else {
    std::cout << "GetInterpolatorFromString(...) Unspecified or unknown interpolator name. "
    << "Defaulting to bspline." << std::cout;
    eInterpolatorType = BSPLINE;
  }

  itk::InterpolateImageFunction<ImageType>::Pointer final_interpolator;

  switch (eInterpolatorType) {
    case NN:
    {
      final_interpolator = dynamic_cast<itk::InterpolateImageFunction<ImageType> *>(
          itk::NearestNeighborInterpolateImageFunction<ImageType>::New().GetPointer());

    }
    break;

    case LINEAR:
    {
      final_interpolator = dynamic_cast<itk::InterpolateImageFunction<ImageType> *>(
          itk::LinearInterpolateImageFunction<ImageType>::New().GetPointer());
    }
    break;

    case BSPLINE:
    {
      itk::BSplineInterpolateImageFunction<ImageType>::Pointer bspline_interpolator = itk::BSplineInterpolateImageFunction<ImageType>::New();
      bspline_interpolator->SetSplineOrder(3); //default is 3
      final_interpolator = dynamic_cast<itk::InterpolateImageFunction<ImageType> *>(bspline_interpolator.GetPointer());
    }
    break;

    case SINC:
    {
      const unsigned int WindowRadius = 5;
      typedef itk::Function::WelchWindowFunction<WindowRadius> WindowFunctionType;
      itk::WindowedSincInterpolateImageFunction<ImageType, WindowRadius, WindowFunctionType>::Pointer sinc_interpolator
              = itk::WindowedSincInterpolateImageFunction<ImageType, WindowRadius, WindowFunctionType>::New();
      final_interpolator = dynamic_cast<itk::InterpolateImageFunction<ImageType> *>(sinc_interpolator.GetPointer());
    }
    break;

    default:
      std::cout << "GetInterpolatorFromString(...) Implementation error!" << std::endl;
      break;
  }

  return final_interpolator;
}


//template <unsigned int DIMENSION>
/*typename*/ MirorrPyramidImplement/*<DIMENSION>*/::ImagePointer
MirorrPyramidImplement/*<DIMENSION>*/
::GetResampledImage(
    /*typename*/ TransformType::Pointer transform,
    //bool resample_moving,
                 ImageType::Pointer image,
                 ImageType::Pointer ref_image,
                 std::string interpolator_name
)
{
  //Create the interpolator used for the final resampling
  itk::InterpolateImageFunction<ImageType>::Pointer final_interpolator = GetInterpolatorFromString(interpolator_name);

  //Create an image resampler
  typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;
  /*typename*/ ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  resampler->SetDefaultPixelValue( 0 );
  resampler->SetUseReferenceImage(true);
  resampler->SetInterpolator(final_interpolator);

//  if( resample_moving )
//  {
//    resampler->SetInput( movingImage );
//    resampler->SetReferenceImage(fixedImage);
//    resampler->SetTransform( dynamic_cast<TransformType*>( transform->GetInverseTransform().GetPointer() ) );
//  }
//  else
//  {
//    resampler->SetInput( fixedImage );
//    resampler->SetReferenceImage(movingImage);
//    resampler->SetTransform( transform );
//  }
  resampler->SetInput( image );
  resampler->SetReferenceImage(ref_image);
  resampler->SetTransform( transform );

  resampler->Update();
  ImagePointer resampledFixedImage = resampler->GetOutput();

  return resampledFixedImage;
}


void
MirorrPyramidImplement
::InitializeRegistrationObjet()
{
  switch (this->registrationType)
  {
    case ERegClassic : {
      this->registration = RegistrationType::New();
      this->registration->SetInterpolator( interpolator );
      this->registration->SetMaskInterpolator( mask_interpolator );
    }
    break;
    case ERegSymmetrical : {
      SymRegistrationType::Pointer reg = SymRegistrationType::New();
      reg->SetInterpolator( interpolator );
      reg->SetMaskInterpolator( mask_interpolator );
      reg->SetFixedImageInterpolator(interpolatorFixed);
      reg->SetFixedMaskInterpolator(fixedmask_interpolator);
      this->registration = reg;
    }
    break;
    default: {
      std::stringstream ss;
      ss
       << "ERROR: MirorrPyramidImplement::InitializeRegistrationObjet() "
       << " Unknown registration type. This should not happen..." << std::endl;
      throw std::runtime_error(ss.str());
    }
  }

  registration->SetVerbosity(verbosity);
}

void
MirorrPyramidImplement
::SetupRegistrationObjetResamplingMode()
{
  // Assuming: not(movingImage.IsNull()) and not(fixedImage.IsNull())
  //=0, EResamplingMiddle, EResamplingMaxResolution, EResamplingMaxSize

  switch(this->resamplingType) {
  case EResamplingBasic:
    this->registration->SetResamplingType(RegistrationType::EResamplingBasic);
    break;
  case EResamplingMiddle:
    this->registration->SetResamplingType(RegistrationType::EResamplingMiddle);
    break;
  case EResamplingMaxResolution: {
      ImageType::SpacingType sF = this->fixedImage->GetSpacing();
      ImageType::SpacingType sM = this->movingImage->GetSpacing();
      double dF = sF.GetSquaredNorm();
      double dM = sM.GetSquaredNorm();

      if (dF < dM) {
        this->registration->SetResamplingType(RegistrationType::EResamplingFixed);
      } else {
        this->registration->SetResamplingType(RegistrationType::EResamplingMoving);
      }
    }
    break;
  case EResamplingMaxSize: {
      ImageType::RegionType::SizeType sF = this->fixedImage->GetLargestPossibleRegion().GetSize();
      ImageType::RegionType::SizeType sM = this->movingImage->GetLargestPossibleRegion().GetSize();
      double vF = 1., vM = 1.;
      for (unsigned int i=0; i<ImageType::RegionType::SizeType::Dimension; i++) {
        vF *= sF[i];
        vM *= sM[i];
      }

      if (vF > vM) {
        this->registration->SetResamplingType(RegistrationType::EResamplingFixed);
      } else {
        this->registration->SetResamplingType(RegistrationType::EResamplingMoving);
      }
    }
    break;
  case EResamplingFixed:
    this->registration->SetResamplingType(RegistrationType::EResamplingFixed);
    break;
  case EResamplingMoving:
    this->registration->SetResamplingType(RegistrationType::EResamplingMoving);
    break;
  default:  {
    std::stringstream ss;
    ss
     << "ERROR: MirorrPyramidImplement::SetupRegistrationObjetResamplingMode() "
     << " Unknown resampling type. This should not happen..." << std::endl;
    throw std::runtime_error(ss.str());
    }
    break;
  }
}

//template <unsigned int DIMENSION>
void
MirorrPyramidImplement
::CreateImagePyramids()
{
  //Generate a decent pyramid schedule
  m_PyramidScheduleTuner->SetFixedSpacing( fixedImage->GetSpacing() );
  m_PyramidScheduleTuner->SetMovingSpacing( movingImage->GetSpacing() );
  m_PyramidScheduleTuner->SetFixedSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  m_PyramidScheduleTuner->SetMovingSize( movingImage->GetLargestPossibleRegion().GetSize() );

  m_PyramidScheduleTuner->Update();

  if( verbosity >= 1 ) {
    std::cout << "\nPyramid Schedule: \n";
    m_PyramidScheduleTuner->PrintPyramidInfo(std::cout, verbosity);
    std::cout << "Building Image pyramid..." << std::endl;
  }

  movingPyramid = ImagePyramidFilterType::New();
  fixedPyramid = ImagePyramidFilterType::New();

  movingPyramid->SetInput( movingImage );
  fixedPyramid->SetInput( fixedImage );

  movingPyramid->SetNumberOfLevels( m_PyramidScheduleTuner->GetMovingSchedule().rows() );
  fixedPyramid->SetNumberOfLevels( m_PyramidScheduleTuner->GetFixedSchedule().rows() );

  movingPyramid->SetSchedule( m_PyramidScheduleTuner->GetMovingSchedule() );
  fixedPyramid->SetSchedule( m_PyramidScheduleTuner->GetFixedSchedule() );

  //Create the image pyramids
  movingPyramid->UpdateLargestPossibleRegion();
  fixedPyramid->UpdateLargestPossibleRegion();
}



#endif
