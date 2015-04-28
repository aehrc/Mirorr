/*=========================================================================
Program: mirorr
Module: itkSymMirorrRegistrationMethod.txx
Author: David Rivest-Henault
Created: 26 Oct 2012

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#ifndef __symMirorrRegistrationMethod_txx_
#define __symMirorrRegistrationMethod_txx_


#include "itkSymMirorrRegistrationMethod.h"
#include "itkEuler3DTransform.h"
#include "vnl_utils.h"
#include <numeric>
#include <iomanip>
#include "itkTransformUtils.h"
#include "itkIOUtils.h"

#include "itkImageRegionIterator.h"

#ifdef VTK_MIRORR
#include "vtkMirorrUtils.h"
#endif

namespace itk
{

template <typename TMovingImage, typename TFixedImage>
SymMirorrRegistrationMethod<TMovingImage,TFixedImage>::SymMirorrRegistrationMethod() :
  m_MovingImageInterpolator(Superclass::m_Interpolator),
  m_MovingMaskInterpolator(Superclass::m_MaskInterpolator),
  m_FixedImageInterpolator(0),
  m_FixedMaskInterpolator(0)
{
}

template <typename TMovingImage, typename TFixedImage>
SymMirorrRegistrationMethod<TMovingImage,TFixedImage>::~SymMirorrRegistrationMethod() {
}

template < typename TMovingImage, typename TFixedImage >
unsigned long
SymMirorrRegistrationMethod<TMovingImage,TFixedImage>
::GetMTime() const
{
  unsigned long mtime = Superclass::GetMTime();
  unsigned long m;

  if (m_MovingMaskInterpolator) //DRH: should be in parent class
  {
    m = m_MovingMaskInterpolator->GetMTime();
    mtime = (m > mtime ? m : mtime);
  }

  if (m_FixedImageInterpolator)
  {
    m = m_FixedImageInterpolator->GetMTime();
    mtime = (m > mtime ? m : mtime);
  }

  if (m_FixedMaskInterpolator)
  {
    m = m_FixedMaskInterpolator->GetMTime();
    mtime = (m > mtime ? m : mtime);
  }

  return mtime;
}

/*
* PrintSelf
*/
template < typename TMovingImage, typename TFixedImage >
void
SymMirorrRegistrationMethod<TMovingImage,TFixedImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "### SymMirorrRegistrationMethod ###" << std::endl;
}

template < typename TMovingImage, typename TFixedImage >
void
SymMirorrRegistrationMethod<TMovingImage,TFixedImage>
::Initialize() throw (ExceptionObject) {
  Superclass::Initialize();

  if( !m_MovingMaskInterpolator ) //DRH: should be in parent class
  {
    itkExceptionMacro(<<"Moving mask interpolator is not present");
  }
  if( !m_FixedImageInterpolator )
  {
    itkExceptionMacro(<<"Fixed image interpolator is not present");
  }
  if( !m_FixedMaskInterpolator )
  {
    itkExceptionMacro(<<"Fixed mask interpolator is not present");
  }
}


template<typename TMovingImage, typename TFixedImage>
typename SymMirorrRegistrationMethod<TMovingImage, TFixedImage>::PointType
SymMirorrRegistrationMethod<TMovingImage, TFixedImage>::ComputeImageCenter(
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

template<typename TMovingImage, typename TFixedImage>
void SymMirorrRegistrationMethod<TMovingImage, TFixedImage>
::ComputeHalfSpaceTransforms(
    const TransformPointer t, TransformPointer t12F, TransformPointer t12B) const
{
  vnl_vector<double> trans = t->GetTranslation().GetVnlVector();

  //-----------------

  vnl_matrix<double> L = t->GetMatrix().GetVnlMatrix();
  vnl_matrix<double> m (4,4,0);
  m(0,0) = L(0,0); m(0,1) = L(0,1); m(0,2) = L(0,2); m(0,3) = trans(0);
  m(1,0) = L(1,0); m(1,1) = L(1,1); m(1,2) = L(1,2); m(1,3) = trans(1);
  m(2,0) = L(2,0); m(2,1) = L(2,1); m(2,2) = L(2,2); m(2,3) = trans(2);
  m(3,3) = 1.;

  vnl_matrix<double> m12, im12, tmp;
  vnl_matrix_sqrt_and_inverse_sqrt(m, m12, tmp, 1e-14, 50);

  //-----------------

  TransformPointer iT = dynamic_cast<TransformType*>(t->CreateAnother().GetPointer());
  t->GetInverse(iT);

  vnl_matrix<double> K = iT->GetMatrix().GetVnlMatrix();
  trans = K * (-trans);

  vnl_matrix<double> mm (4,4,0);
  mm(0,0) = K(0,0); mm(0,1) = K(0,1); mm(0,2) = K(0,2); mm(0,3) = trans(0);
  mm(1,0) = K(1,0); mm(1,1) = K(1,1); mm(1,2) = K(1,2); mm(1,3) = trans(1);
  mm(2,0) = K(2,0); mm(2,1) = K(2,1); mm(2,2) = K(2,2); mm(2,3) = trans(2);
  mm(3,3) = 1.;

  mm = vnl_matrix_inverse<double>(m);
  vnl_matrix_sqrt_and_inverse_sqrt(mm, im12, tmp, 1e-14, 50);

  //-----------------

  if (this->m_Verbosity >=3 )
  {
    std::cout << "ComputeHalfSpaceTransforms m:\n" << m << std::endl;
    std::cout << "ComputeHalfSpaceTransforms m12\n:  " << m12 << std::endl;
    std::cout << "ComputeHalfSpaceTransforms im12\n: " << im12 << std::endl;
  }

  K(0,0) = m12(0,0); K(0,1) = m12(0,1); K(0,2) = m12(0,2);
  K(1,0) = m12(1,0); K(1,1) = m12(1,1); K(1,2) = m12(1,2);
  K(2,0) = m12(2,0); K(2,1) = m12(2,1); K(2,2) = m12(2,2);

  L(0,0) = im12(0,0); L(0,1) = im12(0,1); L(0,2) = im12(0,2);
  L(1,0) = im12(1,0); L(1,1) = im12(1,1); L(1,2) = im12(1,2);
  L(2,0) = im12(2,0); L(2,1) = im12(2,1); L(2,2) = im12(2,2);

  typename TransformType::OutputVectorType _trans12F, _trans12B;
  _trans12F[0] =  m12(0,3); _trans12F[1] =  m12(1,3); _trans12F[2] =  m12(2,3);
  _trans12B[0] = im12(0,3); _trans12B[1] = im12(1,3); _trans12B[2] = im12(2,3);

  t12F->SetMatrix(K);
  t12B->SetMatrix(L);

  t12F->SetTranslation(_trans12F);
  t12B->SetTranslation(_trans12B);
}

template<typename TMovingImage, typename TFixedImage>
void SymMirorrRegistrationMethod<TMovingImage, TFixedImage>::PrintMeanDisplacement(
    PointType p1, PointType p2) const
{
  std::cout << "Displacement is from " << p1 << " to " << p2 << std::endl;
}


template<typename TMovingImage, typename TFixedImage>
void SymMirorrRegistrationMethod<TMovingImage, TFixedImage>::PrintMeanDisplacement(
    PointListType pl1, PointListType pl2) const
{
  if (pl1.size() != pl2.size()) {
    std::cout << "WARNING: cannot print the mean displacement, vectors are of different size." << std::endl;
    return;
  }

  itk::Vector<double,3> d;
  d.Fill(0.);
  const size_t N = pl1.size();
  for (size_t i=0; i<N; i++)  {
    d += pl2[i] - pl1[i];
  }
  d /= N;

  std::cout << "--> INFO: Mean displacement is " << d << ", N=" << N << std::endl;
}

/**
 * The Symmetrical Mirorr registration method -- main loop
 *
 * The current implementation assume that at the beginning this->m_Transform
 * represent a crude alignment of the image (e.g.: translation between the
 * images' barycentres).
 */
template<typename TMovingImage, typename TFixedImage>
void SymMirorrRegistrationMethod<TMovingImage, TFixedImage>::StartOptimization(void)
{
  if (this->m_Verbosity >=3 )
  {
    std::cout << "StartOptimization" << std::endl;
    std::cout.precision(16);

    std::cout << "#DRH - Initial optimised m_Transform (centred at "
        << this->m_Transform->GetCenter() << "):\n"
        << this->m_Transform << std::endl;
  }
  if (this->m_Verbosity >= 1)
  {
    std::cout << "Using symmetric implementation" << std::endl;
  }

  // Get the image centre coordinates
  PointType centerFixedPoint = this->ComputeImageCenter(this->m_FixedImage);
  PointType centerMovingPoint= this->ComputeImageCenter(this->m_MovingImage);

  typename TMovingImage::PointType origin;
  typename TMovingImage::SpacingType spacing;
  typename TMovingImage::SizeType size;
  typedef typename TMovingImage::PointType::VectorType LocalVecType;

  switch(this->m_ResamplingType)
  {
  case Superclass::EResamplingBasic:
  case Superclass::EResamplingMiddle: {
    if (this->m_Verbosity >= 1) {
      std::cout << "Using middle space resampling" << std::endl;
    }
    //Compute the average (half-space) size, spacing, and origin
    typename TMovingImage::PointType::VectorType deltaF = this->m_FixedImage->GetOrigin() - centerFixedPoint;
    typename TMovingImage::PointType::VectorType deltaM = this->m_MovingImage->GetOrigin() - centerMovingPoint;

    origin.GetVnlVector() = ( (deltaF + deltaM) / 2.0 ).GetVnlVector();

    spacing = ( this->m_FixedImage->GetSpacing() + this->m_MovingImage->GetSpacing() ) / 2.0;

    size = ( this->m_FixedImage->GetLargestPossibleRegion().GetSize()
             + this->m_MovingImage->GetLargestPossibleRegion().GetSize() );
    size[0] = (size[0] + 0.5) / 2;
    size[1] = (size[1] + 0.5) / 2;
    size[2] = (size[2] + 0.5) / 2;
  } break;

  case Superclass::EResamplingMoving: {
    if (this->m_Verbosity >= 1) {
      std::cout << "Using movingImage resampling" << std::endl;
    }
    origin.GetVnlVector() = LocalVecType(this->m_MovingImage->GetOrigin() - centerMovingPoint).GetVnlVector();
    spacing               = this->m_MovingImage->GetSpacing();
    size                  = this->m_MovingImage->GetLargestPossibleRegion().GetSize();
  } break;

  case Superclass::EResamplingFixed: {
    if (this->m_Verbosity >= 1) {
      std::cout << "Using fixedImage resampling" << std::endl;
    }
    origin.GetVnlVector() = LocalVecType(this->m_FixedImage->GetOrigin() - centerFixedPoint).GetVnlVector();
    spacing               = this->m_FixedImage->GetSpacing();
    size                  = this->m_FixedImage->GetLargestPossibleRegion().GetSize();
  } break;

  default: {
    std::stringstream ss;
    ss
     << "ERROR: itk::SymMirorrRegistrationMethod::StartOptimization() "
     << " Unknown resampling type. This should not happen..." << std::endl;
    throw std::runtime_error(ss.str());
  } break;
  }

  if (this->m_Verbosity >=2 )
  {
    std::cout << "Size fixed:  " <<  this->m_FixedImage->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "Size moving: " <<  this->m_MovingImage->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "New size:    " <<  size << std::endl;

    std::cout << "Origin fixed:  " <<  this->m_FixedImage->GetOrigin() << std::endl;
    std::cout << "Origin moving: " <<  this->m_MovingImage->GetOrigin() << std::endl;
    std::cout << "New origin:    " <<  origin << std::endl;

    std::cout << "Spacing fixed:  " <<  this->m_FixedImage->GetSpacing() << std::endl;
    std::cout << "Spacing moving: " <<  this->m_MovingImage->GetSpacing() << std::endl;
    std::cout << "New spacing:    " <<  spacing << std::endl;
  }

  // Re-centering m_Transform around zero
  typename TransformType::InputPointType initialCenter = this->m_Transform->GetCenter();
  setCenterToOrigin(this->m_Transform);

  if (this->m_Verbosity >=3 )
  {
    std::cout << "#DRH - Initial optimised m_Transform (centred at "
        << this->m_Transform->GetCenter() << "):\n"
        << this->m_Transform << std::endl;
  }

  // Translation of the centroid of the image at the origin
  typedef Euler3DTransform< double > TranslationType;
  typename TranslationType::Pointer fixedTranslation = TranslationType::New();
  typename TranslationType::Pointer movingTranslation = TranslationType::New();

  typename TranslationType::ParametersType fixedTranslationParam(6);
  typename TranslationType::ParametersType movingTranslationParam(6);
  for (unsigned long k=0; k<3; k++) {
    fixedTranslationParam[k]  = 0.0;
    movingTranslationParam[k] = 0.0;
    fixedTranslationParam[k+3]  = centerFixedPoint[k];
    movingTranslationParam[k+3] = centerMovingPoint[k];
  }

  fixedTranslation->SetParameters(fixedTranslationParam);
  movingTranslation->SetParameters(movingTranslationParam);

  // Modified m_Transform to account for the translation of the
  // centroid of the images at the origin.
  typename TranslationType::Pointer ifixedTranslation  =
      dynamic_cast<TranslationType*>(fixedTranslation->CreateAnother().GetPointer());
  fixedTranslation->GetInverse(ifixedTranslation);

  this->m_Transform->Compose(ifixedTranslation, false);
  this->m_Transform->Compose(movingTranslation, true);

  // Compute the half-space transform (Forward and Backward)
  // At this stage, the centre of rotation are not taken into account
  TransformPointer t12F = dynamic_cast<TransformType*>(this->m_Transform->CreateAnother().GetPointer());
  TransformPointer t12B = dynamic_cast<TransformType*>(this->m_Transform->CreateAnother().GetPointer());
  this->ComputeHalfSpaceTransforms(this->m_Transform, t12F, t12B);

  // Declare two new transform to hold the composed transformation (recentering + half space)
  //  -> centeredT12F, centeredT12B will be set to their correct values later...
  TransformPointer centeredT12F = dynamic_cast<TransformType*>(this->m_Transform->CreateAnother().GetPointer());
  TransformPointer centeredT12B = dynamic_cast<TransformType*>(this->m_Transform->CreateAnother().GetPointer());


  if (this->m_Verbosity >=2 )
  {
    std::cout << "Modif'd translation: " << this->m_Transform->GetTranslation() << std::endl;
    std::cout << " t12F   translation: " << t12F->GetTranslation() << std::endl;
    std::cout << " t12B   translation: " << t12B->GetTranslation() << std::endl;

    std::cout << "\nInitial matrix: \n" << this->m_Transform->GetMatrix() << std::endl;
    std::cout << " t12F   matrix: \n" << t12F->GetMatrix() << std::endl;
    std::cout << " t12B   matrix: \n" << t12B->GetMatrix() << std::endl;
  }

  typedef itk::ResampleImageFilter<TMovingImage, TFixedImage> ResampleFilterType;

  //Create an image resample filter to interpolate the pseudo-FIXED image
  typename ResampleFilterType::Pointer fixedImageResampler = ResampleFilterType::New();
  fixedImageResampler->SetInterpolator(m_FixedImageInterpolator);
  fixedImageResampler->SetTransform(centeredT12F);
  fixedImageResampler->SetDefaultPixelValue(0);
  fixedImageResampler->SetUseReferenceImage(false);
  fixedImageResampler->SetInput(this->m_FixedImage);
  fixedImageResampler->SetOutputOrigin(origin);
  fixedImageResampler->SetOutputSpacing(spacing);
  fixedImageResampler->SetOutputDirection(this->m_FixedImage->GetDirection());
  fixedImageResampler->SetSize(size);

  //Create an image resample filter to interpolate the MOVING image
  typename ResampleFilterType::Pointer movingImageResampler = ResampleFilterType::New();
  movingImageResampler->SetInterpolator(m_MovingImageInterpolator);
  movingImageResampler->SetTransform(centeredT12B);
  movingImageResampler->SetDefaultPixelValue(0);
  movingImageResampler->SetUseReferenceImage(false);
  movingImageResampler->SetInput(this->m_MovingImage);
  movingImageResampler->SetOutputOrigin(origin);
  movingImageResampler->SetOutputSpacing(spacing);
  movingImageResampler->SetOutputDirection(this->m_MovingImage->GetDirection());
  movingImageResampler->SetSize(size);

  typedef itk::ResampleImageFilter<MaskImageType, MaskImageType> MaskResampleFilterType;
  typename MaskResampleFilterType::Pointer fixedmask_resampler =
      MaskResampleFilterType::New();
  fixedmask_resampler->SetInterpolator(m_FixedMaskInterpolator);
  fixedmask_resampler->SetTransform(centeredT12F);
  fixedmask_resampler->SetDefaultPixelValue(0);
  fixedmask_resampler->SetUseReferenceImage(false);
  fixedmask_resampler->SetInput(this->m_FixedMask);
  fixedmask_resampler->SetOutputOrigin(origin);
  fixedmask_resampler->SetOutputSpacing(spacing);
  fixedmask_resampler->SetOutputDirection(this->m_FixedMask->GetDirection());
  fixedmask_resampler->SetSize(size);

  typename MaskResampleFilterType::Pointer movingmask_resampler =
      MaskResampleFilterType::New();
  movingmask_resampler->SetInterpolator(m_FixedMaskInterpolator);
  movingmask_resampler->SetTransform(centeredT12B);
  movingmask_resampler->SetDefaultPixelValue(0);
  movingmask_resampler->SetUseReferenceImage(false);
  movingmask_resampler->SetInput(this->m_MovingMask);
  movingmask_resampler->SetOutputOrigin(origin);
  movingmask_resampler->SetOutputSpacing(spacing);
  movingmask_resampler->SetOutputDirection(this->m_MovingMask->GetDirection());
  movingmask_resampler->SetSize(size);

  this->SetupBlockMatcher();

  //Initial position
  this->m_LastTransformParameters = this->m_Transform->GetParameters(); //m_InitialTransformParameters;

  if (this->m_Verbosity >= 2)
  {
    std::cout << "InitialTransform = " << this->m_LastTransformParameters
        << " FixedParam = "<<this->m_Transform->GetFixedParameters()
        << std::endl;
    this->m_BlockMatcher->displayInfo();
    std::cout << "N_Iterations: " << this->m_MaxIterations << std::endl;
  }

  boost::timer::cpu_timer iteration_timer;
  PointListType fixedInitialPosition, fixedMovedPosition;

  //For N iterationsm_MaxIterations
  const unsigned int maxIterations = static_cast<unsigned int>(std::abs(this->m_MaxIterations));
  for (unsigned int iteration = 0; iteration < maxIterations; iteration++)
  {
    iteration_timer.start();

    //Update the half-space transforms
    centeredT12F->SetParameters(t12F->GetParameters());
    centeredT12B->SetParameters(t12B->GetParameters());
    centeredT12F->Compose(fixedTranslation, false);
    centeredT12B->Compose(movingTranslation, false);

    if (this->m_Verbosity >=3 )
    {
      std::cout << "# DRH - centeredT12F:\n    translation:"
          << centeredT12F->GetTranslation() << "\n    center:"
          << centeredT12F->GetCenter()
          << "\n    Matrix:\n" << centeredT12F->GetMatrix() << std::endl;
      std::cout << "# DRH - centeredT12B:\n    translation:"
          << centeredT12B->GetTranslation() << "\n    center:"
          << centeredT12B->GetCenter()
          << "\n    Matrix:\n" << centeredT12B->GetMatrix() << std::endl;
    }

    //1. Get interpolated fixed image for current transform
    movingImageResampler->Update();
    fixedImageResampler->Update();
    movingmask_resampler->Update();
    fixedmask_resampler->Update();

    if (this->m_Verbosity >=2 )
    {
      typename TMovingImage::Pointer mi = movingImageResampler->GetOutput();
      typename TMovingImage::Pointer fi = fixedImageResampler->GetOutput();
      typedef typename itk::ImageFileWriter<TMovingImage> WriterType;
      typename WriterType::Pointer writer = WriterType::New();

      std::string movFN, fixFN;
      if (iteration==0) {
        movFN = "debugMoving0.nii.gz";
        fixFN = "debugFixed0.nii.gz";
      } else {
        movFN = "debugMoving.nii.gz";
        fixFN = "debugFixed.nii.gz";
      }
      writer->SetFileName(movFN.c_str());
      writer->SetInput(mi);
      writer->Update();
      writer->SetFileName(fixFN.c_str());
      writer->SetInput(fi);
      writer->Update();
      std::cout << "---- Debug images written ----" << std::endl;
      std::cout << "./" << movFN << std::endl;
      std::cout << "./" << fixFN << std::endl;
      std::cout << "------------------------------" << std::endl;
    }

    //2. Block Match the two images
    try
    {
      this->m_BlockMatcher->SetBaseImage(movingImageResampler->GetOutput());
      this->m_BlockMatcher->SetSearchImage(fixedImageResampler->GetOutput());
      this->m_BlockMatcher->SetSearchMask(fixedmask_resampler->GetOutput());
      this->m_BlockMatcher->SetBaseMask(movingmask_resampler->GetOutput());
      this->m_BlockMatcher->GetOffsets(this->m_InitialPositions, this->m_MovedPositions);

      if (this->m_Verbosity >= 2)
      {
        std::cout << "--> INFO: @idx=0, "; PrintMeanDisplacement(this->m_InitialPositions[0], this->m_MovedPositions[0]);
        PrintMeanDisplacement(this->m_InitialPositions, this->m_MovedPositions);
      }

      //---------------------------------------
      //---------------------------------------
#ifdef DEBUG_GPU_CPU_CONSISTENCY
      unsigned int nDiff = 0;
      std::cout << std::endl;
      std::cout << "#######################" << std::endl;
      std::cout << "####  Scores diff: ####" << std::endl;
      for (unsigned int k=0; k< this->m_MovedPositions.size(); k++)
      {
        typename TMovingImage::PointType basePt   = this->m_InitialPositions[k];
        typename TMovingImage::IndexType baseIdx  = this->m_BlockMatcher->pointToBlockPosition(basePt);

        typename TFixedImage::PointType searchPt = this->m_MovedPositions[k];
        typename TFixedImage::IndexType searchIdx= this->m_BlockMatcher->pointToBlockPosition(searchPt);

        try {
          float score = this->m_BlockMatcher->testNormalizedCorrelation(baseIdx, searchIdx);

          if (score != this->m_BlockMatcher->scoresOut[k]) {
            nDiff++;
            float diff = score - this->m_BlockMatcher->scoresOut[k];
            std::cout << "### Diff at index " << k << std::endl;
            std::cout << "###   From " << basePt << " to " << searchPt << "(mm)" << std::endl;
            std::cout << "###   From " << baseIdx << " to " << searchIdx << "(idx)" << std::endl;
            std::cout << std::setprecision(20);
            std::cout << "###   Score: " << score << " (Post), " << this->m_BlockMatcher->scoresOut[k] << " (BM), diff = " << diff << std::endl;
            std::cout << std::setprecision(4);
          }
        } catch (...) {
          std::cout << "###### 666 ######" << std::endl;
        }
      }
      std::cout << "### N diff: " << nDiff << "/" << this->m_BlockMatcher->scoresOut.size() << std::endl;
      std::cout << "#######################" << std::endl;

      const unsigned int K_MAX = 100;
      const unsigned int K = (this->m_MovedPositions.size() > K_MAX) ? K_MAX : this->m_MovedPositions.size();

#if 1
      std::cout << "\n#######################" << std::endl;
      std::cout << "####    Scores:    ####" << std::endl;
      std::cout << std::setprecision(20);
      for (unsigned int k=0; k<K; k++)
      {
        std::cout << "K=" << std::setw(4) << std::fixed << k << ", " << this->m_BlockMatcher->scoresOut[k] << std::endl;
      }
      std::cout << std::setprecision(4);
      std::cout << "#######################" << std::endl;
#endif

#if 1
      std::cout << "\n#######################" << std::endl;
      std::cout << "####  2nd Scores:  ####" << std::endl;
      std::cout << std::setprecision(20);
      for (unsigned int k=0; k<K; k++)
      {
        std::cout << "K=" << std::setw(4) << std::fixed << k << ", " << this->m_BlockMatcher->scoresOut2nd[k] << std::endl;
      }
      std::cout << std::setprecision(4);
      std::cout << "#######################" << std::endl;
#endif
#endif //DEBUG_GPU_CPU_CONSISTENCY
      //---------------------------------------
      //---------------------------------------

      fixedInitialPosition.clear();
      fixedMovedPosition.clear();
      this->m_BlockMatcher->SwapImageOrder();
      this->m_BlockMatcher->GetOffsets(fixedInitialPosition, fixedMovedPosition);

      this->m_InitialPositions.insert(this->m_InitialPositions.end(),
          fixedMovedPosition.begin(), fixedMovedPosition.end());

      this->m_MovedPositions.insert(this->m_MovedPositions.end(),
                fixedInitialPosition.begin(), fixedInitialPosition.end());

      if (this->m_Verbosity >= 2)
      {
        std::cout << "--> INFO: "; PrintMeanDisplacement(fixedInitialPosition[0], fixedMovedPosition[0]);
        PrintMeanDisplacement(fixedInitialPosition, fixedMovedPosition);
        std::cout << "--> INFO: "; PrintMeanDisplacement(this->m_InitialPositions[0], this->m_MovedPositions[0]);
        PrintMeanDisplacement(this->m_InitialPositions, this->m_MovedPositions);
      }
    }
    catch (std::exception &e)
    {
      std::cerr << "Exception caught in SymMirorrRegistrationMethod::StartOptimization() while updating BlockMatcher: " << e.what()
          << std::endl;
      exit(-1);
    }

    //Get the update to the transform
    this->OptimizeTransform(t12F, t12B);

    //4a. Check if parameters have changed
    double squareDifferenceInPosition = 0;
    for (unsigned int ii = 0; ii < this->m_LastTransformParameters.GetSize(); ++ii)
    {
      double diff = this->m_LastTransformParameters[ii]
          - this->m_Transform->GetParameters()[ii];
      squareDifferenceInPosition += diff * diff;
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
    this->m_LastTransformParameters = this->m_Transform->GetParameters();

    //4c. Terminate if movement less than 1/10 voxels or params hardly change
    if (squareDifferenceInPosition < 1e-3)
    {
      if (this->m_Verbosity >= 1)
        std::cout << "Early breakout as movement small." << std::endl;
      break;
    }
  }

  // Update the centred transform since they are used to compute the global transform
  // WARNING: must be as above!!! @see_comment: "Update the half-space transforms"
  centeredT12F->SetParameters(t12F->GetParameters());
  centeredT12B->SetParameters(t12B->GetParameters());
  centeredT12F->Compose(fixedTranslation, false);
  centeredT12B->Compose(movingTranslation, false);

  typename TransformType::Pointer invT12B = TransformType::New();
  centeredT12B->GetInverse(invT12B);

  // Compose the final transform
  this->m_Transform->SetIdentity();
  this->m_Transform->Compose(invT12B, true);
  this->m_Transform->Compose(centeredT12F, false);
  this->m_Transform->GetParameters();

  if (this->m_Verbosity >= 3) {
    std::cout << "#DRH - Final t12F:\n"
                << t12F << std::endl;
    std::cout << "#DRH - Final t12B:\n"
              << t12B << std::endl;
    std::cout << "#DRH - Final centeredT12F:\n"
              << centeredT12F << std::endl;
    std::cout << "#DRH - Final centeredT12B:\n"
              << centeredT12B << std::endl;
    std::cout << "#DRH - Final optimised m_Transform (centred at 0):\n"
              << this->m_Transform << std::endl;
  }

  // Re-center the transform around the standard center
  itk::updateCenter(this->m_Transform, initialCenter);
  this->m_LastTransformParameters = this->m_Transform->GetParameters();

  if (this->m_Verbosity >= 3) {
    std::cout << "#DRH - Final optimised m_Transform (centred at "
        << this->m_Transform->GetCenter() << "):\n"
        << this->m_Transform << std::endl;
  }

#ifdef VTK_MIRORR
#if 1
  typename TransformType::Pointer invT = TransformType::New();
  this->m_Transform->GetInverse(invT);

  const size_t N =  this->m_InitialPositions.size();
  for (size_t i=0; i<N; i++) {
    //centerFixedPoint centerMovingPoint
    this->m_InitialPositions[i][0] -= centerFixedPoint[0];
    this->m_InitialPositions[i][1] -= centerFixedPoint[1];
    this->m_InitialPositions[i][2] -= centerFixedPoint[2];
    this->m_MovedPositions[i][0] -= centerFixedPoint[0];
    this->m_MovedPositions[i][1] -= centerFixedPoint[1];
    this->m_MovedPositions[i][2] -= centerFixedPoint[2];

    this->m_InitialPositions[i] = this->m_Transform->TransformPoint(this->m_InitialPositions[i]);
    this->m_MovedPositions[i]   = this->m_Transform->TransformPoint(this->m_MovedPositions[i]);
  }
#endif
#endif
}


/**
 * DrH, 2012-11-09
  *
 * WARNING: At the moment (2012-11-09) this implementation is only
 * completed for the 'Euler3DTransform' and 'affine' transformation type
 */
template<typename TMovingImage, typename TFixedImage>
void
SymMirorrRegistrationMethod<TMovingImage, TFixedImage>
::OptimizeTransform(TransformPointer t12F, TransformPointer t12B)
{
  if (this->m_Verbosity >=3 ) {
    std::cout << "OptimizeTransform" << std::endl;
    std::cout.precision(16);
  }

  const unsigned int kDimension = 3;
  typedef typename itk::MatrixOffsetTransformBase<double, kDimension,
      kDimension> LocalTransformType;
  LocalTransformType::Pointer final_transform =
      dynamic_cast<LocalTransformType*>(this->m_Transform->CreateAnother().GetPointer());

  LocalTransformType::Pointer final_transformF =
        dynamic_cast<LocalTransformType*>(this->m_Transform->CreateAnother().GetPointer());

  LocalTransformType::Pointer final_transformB =
        dynamic_cast<LocalTransformType*>(this->m_Transform->CreateAnother().GetPointer());

  PointType rotation_centre = this->m_Transform->GetCenter();
  final_transform->SetIdentity(); // This must be before SetCenter, otherwise center will be reset to 0,0,0
  final_transform->SetCenter(rotation_centre);
  ParametersType last_parameters(final_transform->GetParameters());

  final_transformF->SetIdentity();
  final_transformF->SetCenter(rotation_centre);

  final_transformB->SetIdentity();
  final_transformB->SetCenter(rotation_centre);

  // LTS
  double fraction = 1.0; //This changes to 0.5 later - initially we use all points

  std::string transform_type = this->m_Transform->GetNameOfClass();

  // Affine transform gets performs a single pass
  //While change in transform is sufficiently low
  double last_delta_param_tol = 0;
  double delta_param_tol = 1;
  int iter = -1;
  //Check tolerance met, max no. iterations, flip flopping between solutions
  while (delta_param_tol > 1e-3 && iter < 100
      && fabs(delta_param_tol - last_delta_param_tol) > 1e-3) //Flip flopping
  {
    ++iter;

    if (this->m_InitialPositions.empty())
    {
      break;
    }

    if (this->m_Verbosity >=3 )
    {
      std::cout << "Transform for iteration no" << iter
          << ":\n" << final_transform << std::endl;
    }

    //Get residuals
    std::vector<double> residuals(this->m_InitialPositions.size(), 0.0);
    for (unsigned int point = 0; point < this->m_InitialPositions.size(); ++point)
    {
      PointType transformed_pointF = final_transformF->TransformPoint(this->m_InitialPositions[point]);
      PointType transformed_pointB = final_transformB->TransformPoint(this->m_MovedPositions[point]);
      for (unsigned int col = 0; col < kDimension; ++col)
      {
        //Transform m_InitialPositions
        double dist = transformed_pointF[col] - transformed_pointB[col];
        residuals[point] += dist * dist;
      }
    }

    //Find nth residual
    /**
     *
     * DRH -- Note: stable-sorting the point vectors according to the corresponding
     * residual values would most likely improve the symmetry of the algorithm.
     *
     */
    double nth_residual = 0.0;
    if (fraction > 0.99999) {
      //David RH, 2012-11-09: 'if' necessary since std::nth_element does
      //nothing when fraction ~= 1.0.
      nth_residual = *std::max_element(residuals.begin(), residuals.end());

      if (this->m_Verbosity >=2 )
      {
        std::streamsize p = std::cout.precision();
        std::cout.precision(16);
        std::cout << "sum residual: " << std::accumulate(residuals.begin(), residuals.end(), 0.) << std::endl;
        std::cout.precision(p);
      }
    }
    else
    {
      //Sort points in order of residuals and retain only the 'fraction' (~50%) of points
      //with the lowest residual
      std::vector<double> residuals_sorted(residuals);

      //if( residuals_sorted.size() > 0)
      {
        unsigned int len = static_cast<unsigned int>(fraction
            * residuals_sorted.size());
        len = std::min<unsigned int>( len, residuals_sorted.size()-1 );
        std::nth_element(residuals_sorted.begin(),
            residuals_sorted.begin() + len, residuals_sorted.end());
        //nth_residual = residuals_sorted[len - 1];
        nth_residual = residuals_sorted[len];

        if (this->m_Verbosity >=2 )
        {
          std::streamsize p = std::cout.precision();
          std::cout.precision(16);
          std::cout << "partial sum residual: " << std::accumulate(residuals_sorted.begin(), residuals_sorted.begin()+len, 0.) << std::endl;
          std::cout.precision(p);
        }
      }
    }


    if (this->m_Verbosity >=2 )
    {
      std::cout << "nth_residual: " << nth_residual << std::endl;
    }

    //Extract points with lowest residual
    PointListType initialSubset;
    initialSubset.reserve(this->m_InitialPositions.size());
    PointListType movedSubset;
    movedSubset.reserve(this->m_InitialPositions.size());
    for (unsigned int point = 0; point < this->m_InitialPositions.size(); ++point)
      if (residuals[point] <= nth_residual)
      {
        initialSubset.push_back(this->m_InitialPositions[point]);
        movedSubset.push_back(this->m_MovedPositions[point]);
      }

    //After first iteration, we reduce the fraction to 0.5
    fraction = 0.5;

    //Begin LS ************************
    //Get barycentres
    vnl_vector<double> initialMean = this->GetBaryCentre(initialSubset);
    vnl_vector<double> movedMean = this->GetBaryCentre(movedSubset);

    if (this->m_Verbosity >=2 )
    {
      std::cout << "initialMean: " << initialMean << std::endl;
      std::cout << "movedMean: " << movedMean << std::endl;
    }

    //Convert subset to barycentric coordinates in mm
    for (unsigned int point = 0; point < initialSubset.size(); ++point)
      for (unsigned int axis = 0; axis < kDimension; ++axis)
      {
        initialSubset[point][axis] = (initialSubset[point][axis]
            - initialMean[axis]);    // * spacing[axis];
        movedSubset[point][axis] =
            (movedSubset[point][axis] - movedMean[axis]);  // * spacing[axis];
      }

    //Branch for Affine / rigid *****************************
    //Only a parameters typ object should leave this part
    if (this->m_Verbosity >= 3)
    {
      std::cout << "\n  Inner iteration: " << iter << std::endl;
      std::cout << std::setprecision(5) << "  Init centre: " << initialMean
          << "   Moved centre: " << movedMean << std::endl;
      std::cout << "  Num points: " << initialSubset.size() << " of "
          << this->m_InitialPositions.size() << "  Nth residual: " << nth_residual
          << std::endl;
      std::cout << "  TransformType: " << transform_type << std::endl;
    }

    ParametersType new_parameters(last_parameters);
    ParametersType new_parametersF(last_parameters);
    ParametersType new_parametersB(last_parameters);

    if (0 == transform_type.compare("AffineTransform"))
    {
      /**
       * Compute the forward transformation
       */

      //Get summed outer product of x and u
      vnl_matrix<double> outer_xy = this->CovarianceOfPoints(movedSubset, initialSubset); //ORder is important here
      vnl_matrix<double> outer_xx = this->CovarianceOfPoints(initialSubset, initialSubset);

      //Get rotationF matrix and translation
      vnl_matrix<double> rotationF = outer_xy;
      if (fabs(vnl_determinant(outer_xx)) < 1e-4) //Check for low determinants
        std::cout
            << "\nERROR: Problem inverting outer_xx as its determinant is low:\n"
            << outer_xx << "Treating it as the identity matrix.\n"
            << std::endl;
      else
        rotationF = rotationF * vnl_matrix_inverse<double>(outer_xx);

      if (this->m_Verbosity >= 3)
      {
        std::cout << "Field" << std::endl;
        for (unsigned int ii = 0; ii < movedSubset.size(); ++ii)
        {
          std::cout << movedSubset[ii][0] << " " << movedSubset[ii][1] << " "
              << movedSubset[ii][2] << " " << initialSubset[ii][0] << " "
              << initialSubset[ii][1] << " " << initialSubset[ii][2] << " "
              << std::endl;
        }
        std::cout << "----" << std::endl;

        std::cout << "outer_xy:\n" << outer_xy;
        std::cout << "outer_xx:\n" << outer_xx;
        std::cout << "Rotation:\n" << rotationF;
      }

      /**
       * Compute the backward transformation
       */

      //Get summed outer product of x and u
      outer_xy = this->CovarianceOfPoints(initialSubset, movedSubset); //ORder is important here
      outer_xx = this->CovarianceOfPoints(movedSubset, movedSubset);

      //Get rotationF matrix and translation
      vnl_matrix<double> rotationB = outer_xy;
      if (fabs(vnl_determinant(outer_xx)) < 1e-4) //Check for low determinants
        std::cout
            << "\nERROR: Problem inverting outer_xx as its determinant is low:\n"
            << outer_xx << "Treating it as the identity matrix.\n"
            << std::endl;
      else
        rotationB = rotationB * vnl_matrix_inverse<double>(outer_xx);

      vnl_vector<double> translation  = movedMean - initialMean;
      vnl_vector<double> translationF =  translation;
      vnl_vector<double> translationB = -translation;

      if (this->m_Verbosity >=3 ) {
        std::cout << "translation F:\n" << translationF << "\n\n";
        std::cout << "outer_xy B   :\n" << outer_xy;
        std::cout << "outer_xx B   :\n" << outer_xx;
        std::cout << "Rotation B   :\n" << rotationB;
        std::cout << "translation B:\n" << translationB << "\n\n";
      }

      //TODO - check this assignment is correct
      //Extract parameters
      //Rotation
      new_parameters[0] = rotationF[0][0];
      new_parameters[1] = rotationF[0][1];
      new_parameters[2] = rotationF[0][2];
      new_parameters[3] = rotationF[1][0];
      new_parameters[4] = rotationF[1][1];
      new_parameters[5] = rotationF[1][2];
      new_parameters[6] = rotationF[2][0];
      new_parameters[7] = rotationF[2][1];
      new_parameters[8] = rotationF[2][2];
      //Translation
      new_parameters[9] = translation[0];
      new_parameters[10] = translation[1];
      new_parameters[11] = translation[2];

      //Extract parameters
      //Rotation
      new_parametersF[0] = rotationF[0][0];
      new_parametersF[1] = rotationF[0][1];
      new_parametersF[2] = rotationF[0][2];
      new_parametersF[3] = rotationF[1][0];
      new_parametersF[4] = rotationF[1][1];
      new_parametersF[5] = rotationF[1][2];
      new_parametersF[6] = rotationF[2][0];
      new_parametersF[7] = rotationF[2][1];
      new_parametersF[8] = rotationF[2][2];
      //Translation
      new_parametersF[9] = translationF[0];
      new_parametersF[10] = translationF[1];
      new_parametersF[11] = translationF[2];

      //Extract parameters
      //Rotation
      new_parametersB[0] = rotationB[0][0];
      new_parametersB[1] = rotationB[0][1];
      new_parametersB[2] = rotationB[0][2];
      new_parametersB[3] = rotationB[1][0];
      new_parametersB[4] = rotationB[1][1];
      new_parametersB[5] = rotationB[1][2];
      new_parametersB[6] = rotationB[2][0];
      new_parametersB[7] = rotationB[2][1];
      new_parametersB[8] = rotationB[2][2];
      //Translation
      new_parametersB[9] = translationB[0];
      new_parametersB[10] = translationB[1];
      new_parametersB[11] = translationB[2];

      rotation_centre = initialMean.begin();
      final_transformF->SetCenter(rotation_centre);
      final_transformF->SetParameters(new_parametersF);

      rotation_centre = movedMean.begin();
      final_transformB->SetCenter(rotation_centre);
      final_transformB->SetParameters(new_parametersB);

      if (this->m_Verbosity >=3 ) {
        std::cout << "# DRH - final_transformF:\n    translation:"
            << final_transformF->GetTranslation() << "\n    center:"
            << final_transformF->GetCenter()
            << "\n    Matrix:\n" << final_transformF->GetMatrix() << std::endl;
        std::cout << "# DRH - final_transformB:\n    translation:"
            << final_transformB->GetTranslation() << "\n    center:"
            << final_transformB->GetCenter()
            << "\n    Matrix:\n" << final_transformB->GetMatrix() << std::endl;
      }

      /**
       * Average forward and backward matrices to improve numerical stability & symmetry
       */
      vnl_matrix<double> K = final_transformF->GetMatrix().GetVnlMatrix();
      vnl_matrix<double> L =  vnl_matrix_inverse<double>(final_transformB->GetMatrix().GetVnlMatrix());

      if (this->m_Verbosity >=3 ) {
        std::cout << "K F:\n" << K << std::endl;
        std::cout << "L F:\n" << L << std::endl;
      }

      K = (K + L) / 2.0;

      vnl_matrix<double> M = vnl_matrix_inverse<double>(final_transformF->GetMatrix().GetVnlMatrix());
      vnl_matrix<double> N = final_transformB->GetMatrix().GetVnlMatrix();

      if (this->m_Verbosity >=3 ) {
        std::cout << "M B:\n" << M << std::endl;
        std::cout << "N B:\n" << N << std::endl;
      }

      M = (M + N) / 2.0;

      final_transformF->SetMatrix(K);
      final_transformB->SetMatrix(M);

      /**
       * Re-centre the transformation about the origin
       */
      typename TransformType::OutputVectorType offF = final_transformF->GetOffset();
      typename TransformType::OutputVectorType offB = final_transformB->GetOffset();

      typename TransformType::InputPointType zerosCenter;
      for(unsigned int i=0; i<3; ++i)
        zerosCenter[i] = 0.0;
      final_transformF->SetCenter(zerosCenter);
      final_transformB->SetCenter(zerosCenter);

      final_transformF->SetTranslation(offF);
      final_transformB->SetTranslation(offB);

      /**
       * Compute the half space transforms
       * The forward and backward^-1 matrices are averaged for numerical stability
       */
      LocalTransformType::Pointer iF =
              dynamic_cast<LocalTransformType*>(this->m_Transform->CreateAnother().GetPointer());
      LocalTransformType::Pointer iB =
                    dynamic_cast<LocalTransformType*>(this->m_Transform->CreateAnother().GetPointer());

      final_transformF->GetInverse(iF);
      final_transformB->GetInverse(iB);

      if (this->m_Verbosity >=3 ) {
        std::cout << "# DRH - final_transformF:\n    translation:"
              << final_transformF->GetTranslation()
              << "\n    Matrix:\n" << final_transformF->GetMatrix() << std::endl;
        std::cout << "# DRH - final_transformB:\n    translation:"
            << final_transformB->GetTranslation()
            << "\n    Matrix:\n" << final_transformB->GetMatrix() << std::endl;
      }

      // Forward half space transform
      K = final_transformF->GetMatrix().GetVnlMatrix();
      L = iB->GetMatrix().GetVnlMatrix();
      K = (K + L) / 2.0;

      vnl_vector<double> transF = final_transformF->GetTranslation().GetVnlVector();
      vnl_vector<double> transB = iB->GetTranslation().GetVnlVector();
      vnl_vector<double> trans  = (transF + transB) / 2.0;

      vnl_matrix<double> m (4,4,0);
      m(0,0) = K(0,0); m(0,1) = K(0,1); m(0,2) = K(0,2); m(0,3) = trans(0);
      m(1,0) = K(1,0); m(1,1) = K(1,1); m(1,2) = K(1,2); m(1,3) = trans(1);
      m(2,0) = K(2,0); m(2,1) = K(2,1); m(2,2) = K(2,2); m(2,3) = trans(2);
      m(3,3) = 1.;

      vnl_matrix<double> m12, im12, tmp;
      vnl_matrix_sqrt_and_inverse_sqrt(m, m12, tmp, 1e-14, 50);

      // Backward half space transform
      K = iF->GetMatrix().GetVnlMatrix();
      L = final_transformB->GetMatrix().GetVnlMatrix();
      K = (K + L) / 2.0;

      transF = iF->GetTranslation().GetVnlVector();
      transB = final_transformB->GetTranslation().GetVnlVector();
      trans  = (transF + transB) / 2.0;

      vnl_matrix<double> mm (4,4,0);
      mm(0,0) = K(0,0); mm(0,1) = K(0,1); mm(0,2) = K(0,2); mm(0,3) = trans(0);
      mm(1,0) = K(1,0); mm(1,1) = K(1,1); mm(1,2) = K(1,2); mm(1,3) = trans(1);
      mm(2,0) = K(2,0); mm(2,1) = K(2,1); mm(2,2) = K(2,2); mm(2,3) = trans(2);
      mm(3,3) = 1.;

      vnl_matrix_sqrt_and_inverse_sqrt(mm, im12, tmp, 1e-14, 50);

      //Final assignment
      K(0,0) = m12(0,0); K(0,1) = m12(0,1); K(0,2) = m12(0,2);
      K(1,0) = m12(1,0); K(1,1) = m12(1,1); K(1,2) = m12(1,2);
      K(2,0) = m12(2,0); K(2,1) = m12(2,1); K(2,2) = m12(2,2);

      L(0,0) = im12(0,0); L(0,1) = im12(0,1); L(0,2) = im12(0,2);
      L(1,0) = im12(1,0); L(1,1) = im12(1,1); L(1,2) = im12(1,2);
      L(2,0) = im12(2,0); L(2,1) = im12(2,1); L(2,2) = im12(2,2);

      typename TransformType::OutputVectorType _trans12F, _trans12B;
      _trans12F[0] =  m12(0,3); _trans12F[1] =  m12(1,3); _trans12F[2] =  m12(2,3);
      _trans12B[0] = im12(0,3); _trans12B[1] = im12(1,3); _trans12B[2] = im12(2,3);

      final_transformF->SetMatrix(K);
      final_transformB->SetMatrix(L);

      final_transformF->SetTranslation(_trans12F);
      final_transformB->SetTranslation(_trans12B);

      if (this->m_Verbosity >=3 ) {
        std::cout << "m:\n" << m << std::endl;
        std::cout << "mm:\n" << mm << std::endl;
        std::cout << "m12:\n" << m12 << std::endl;
        std::cout << "im12:\n" << im12 << std::endl;

        //vnl_matrix_inverse<T>(m12)
        //vnl_inverse(outer_xx);
        vnl_matrix<double> testF = m12 * vnl_inverse(im12);
        vnl_matrix<double> testB = im12 * vnl_inverse(m12);
        std::cout << "testF:\n" << testF << std::endl;
        std::cout << "testB:\n" << testB << std::endl;
      }

      //Get summed change in parameters
      last_delta_param_tol = delta_param_tol;
      delta_param_tol = this->GetParametersDelta(last_parameters, new_parameters);
    }
    else if (0 == transform_type.compare("Euler3DTransform"))
    {
      vnl_vector<double> euler_angles, euler_anglesF, euler_anglesB;
      vnl_vector_fixed<double,3> axis;

      //Create quaternion matrix for each point pair and sum Q^T . Q
      vnl_quaternion<double> quaternion = this->CalculateRigidQuaternion(initialSubset, movedSubset);

      euler_angles = quaternion.rotation_euler_angles();

      axis = quaternion.axis();
      double angle = quaternion.angle();

      /**
       * Compute forward and backward transform
       */
      vnl_quaternion<double> quaternion12F(axis,  angle);
      vnl_quaternion<double> quaternion12B(axis, -angle);

      euler_anglesF = quaternion12F.rotation_euler_angles();
      euler_anglesB = quaternion12B.rotation_euler_angles();

      vnl_vector<double> translation  = movedMean - initialMean;
      vnl_vector<double> translationF =  translation;
      vnl_vector<double> translationB = -translation;

      //Extract parameters
      for (int ii = 0; ii < 3; ++ii) {
        new_parameters[ii]  = euler_angles[ii];
        new_parametersF[ii] = euler_anglesF[ii];
        new_parametersB[ii] = euler_anglesB[ii];

        new_parameters[ii + 3]  = translation[ii];
        new_parametersF[ii + 3] = translationF[ii];
        new_parametersB[ii + 3] = translationB[ii];
      }

      rotation_centre = initialMean.begin();
      final_transformF->SetCenter(rotation_centre);
      final_transformF->SetParameters(new_parametersF);

      rotation_centre = movedMean.begin();
      final_transformB->SetCenter(rotation_centre);
      final_transformB->SetParameters(new_parametersB);

      /**
       * Re-centre the transformation about the origin
       */
      vnl_matrix<double> K = final_transformF->GetMatrix().GetVnlMatrix();
      vnl_matrix<double> L = final_transformB->GetInverseMatrix().GetVnlMatrix();
      K = (K + L) / 2.0;
      K = orthogonalize<double,3>(K);

      L = K.transpose();

      final_transformF->SetMatrix(K);
      final_transformB->SetMatrix(L);
      //-----------------

      typename TransformType::OutputVectorType offF = final_transformF->GetOffset();
      typename TransformType::OutputVectorType offB = final_transformB->GetOffset();

      typename TransformType::InputPointType zerosCenter;
      for(unsigned int i=0; i<3; ++i)
        zerosCenter[i] = 0.0;

      final_transformF->SetCenter(zerosCenter);
      final_transformB->SetCenter(zerosCenter);

      final_transformF->SetTranslation(offF);
      final_transformB->SetTranslation(offB);

      /**
       * Compute the half space transforms
       * The forward and backward^-1 matrices are averaged for numerical stability
       */
      LocalTransformType::Pointer iF =
              dynamic_cast<LocalTransformType*>(this->m_Transform->CreateAnother().GetPointer());
      LocalTransformType::Pointer iB =
                    dynamic_cast<LocalTransformType*>(this->m_Transform->CreateAnother().GetPointer());

      final_transformF->GetInverse(iF);
      final_transformB->GetInverse(iB);

      if (this->m_Verbosity >=3 ) {
        std::cout << "# DRH - final_transformF:\n    translation:"
              << final_transformF->GetTranslation()
              << "\n    Matrix:\n" << final_transformF->GetMatrix() << std::endl;
        std::cout << "# DRH - inv final_transformF:\n    translation:"
              << iF->GetTranslation()
              << "\n    Matrix:\n" << iF->GetMatrix() << std::endl;
        std::cout << "# DRH - final_transformB:\n    translation:"
            << final_transformB->GetTranslation()
            << "\n    Matrix:\n" << final_transformB->GetMatrix() << std::endl;
        std::cout << "# DRH - inv final_transformB:\n    translation:"
            << iB->GetTranslation()
            << "\n    Matrix:\n" << iB->GetMatrix() << std::endl;
      }

      // Forward half space transform
      K = final_transformF->GetMatrix().GetVnlMatrix();
      L = iB->GetMatrix().GetVnlMatrix();
      K = (K + L) / 2.0;
      K = orthogonalize<double,3>(K);

      vnl_vector<double> transF = final_transformF->GetTranslation().GetVnlVector();
      vnl_vector<double> transB = iB->GetTranslation().GetVnlVector();
      vnl_vector<double> trans  = (transF + transB) / 2.0;

      vnl_matrix<double> m (4,4,0);
      m(0,0) = K(0,0); m(0,1) = K(0,1); m(0,2) = K(0,2); m(0,3) = trans(0);
      m(1,0) = K(1,0); m(1,1) = K(1,1); m(1,2) = K(1,2); m(1,3) = trans(1);
      m(2,0) = K(2,0); m(2,1) = K(2,1); m(2,2) = K(2,2); m(2,3) = trans(2);
      m(3,3) = 1.;

      vnl_matrix<double> m12, im12, tmp;
      vnl_matrix_sqrt_and_inverse_sqrt(m, m12, tmp, 1e-14, 50);

      // Backward half space transform
      K = iF->GetMatrix().GetVnlMatrix();
      L = final_transformB->GetMatrix().GetVnlMatrix();
      K = (K + L) / 2.0;
      K = orthogonalize<double,3>(K);

      transF = iF->GetTranslation().GetVnlVector();
      transB = final_transformB->GetTranslation().GetVnlVector();
      trans  = (transF + transB) / 2.0;

      vnl_matrix<double> mm (4,4,0);
      mm(0,0) = K(0,0); mm(0,1) = K(0,1); mm(0,2) = K(0,2); mm(0,3) = trans(0);
      mm(1,0) = K(1,0); mm(1,1) = K(1,1); mm(1,2) = K(1,2); mm(1,3) = trans(1);
      mm(2,0) = K(2,0); mm(2,1) = K(2,1); mm(2,2) = K(2,2); mm(2,3) = trans(2);
      mm(3,3) = 1.;

      vnl_matrix_sqrt_and_inverse_sqrt(mm, im12, tmp, 1e-14, 50);

      //Final assignment
      K(0,0) = m12(0,0); K(0,1) = m12(0,1); K(0,2) = m12(0,2);
      K(1,0) = m12(1,0); K(1,1) = m12(1,1); K(1,2) = m12(1,2);
      K(2,0) = m12(2,0); K(2,1) = m12(2,1); K(2,2) = m12(2,2);

      L(0,0) = im12(0,0); L(0,1) = im12(0,1); L(0,2) = im12(0,2);
      L(1,0) = im12(1,0); L(1,1) = im12(1,1); L(1,2) = im12(1,2);
      L(2,0) = im12(2,0); L(2,1) = im12(2,1); L(2,2) = im12(2,2);

      typename TransformType::OutputVectorType _trans12F, _trans12B;
      _trans12F[0] =  m12(0,3); _trans12F[1] =  m12(1,3); _trans12F[2] =  m12(2,3);
      _trans12B[0] = im12(0,3); _trans12B[1] = im12(1,3); _trans12B[2] = im12(2,3);

      try {
        final_transformF->SetMatrix(K);
        final_transformB->SetMatrix(L);
      } catch (itk::ExceptionObject &) {
        // DRH: In rare occasions, the K and L matrices can be too far from
        // being orthogonal to be accepted by ITK, so we orthogonalize them.
        // The following operations reduce the symmetry of the algorithm, that's why
        // they are not applied systematically.
        std::cout << "WARNING: Questionable accuracy in "
            "SymMirorrRegistrationMethod::OptimizeTransform(...). "
            "Re-orthogonalizing." << std::endl;
        K = orthogonalize<double,3>(K);
        L = orthogonalize<double,3>(L);
        final_transformF->SetMatrix(K);
        final_transformB->SetMatrix(L);
      }

      final_transformF->SetTranslation(_trans12F);
      final_transformB->SetTranslation(_trans12B);

      if (this->m_Verbosity >=3 ) {
        std::cout << "m:\n" << m << std::endl;
        std::cout << "mm:\n" << mm << std::endl;
        std::cout << "m12\n:  " << m12 << std::endl;
        std::cout << "im12\n: " << im12 << std::endl;
      }

      //Get summed change in parameters
      last_delta_param_tol = delta_param_tol;
      delta_param_tol = this->GetParametersDelta(new_parameters, last_parameters);
    }
    else if (0 == transform_type.compare("QuaternionRigid"))
    {
      //Create quaternion matrix for each point pair and sum Q^T . Q
      vnl_quaternion<double> quaternion = this->CalculateRigidQuaternion(
          initialSubset, movedSubset);

      //Get rotation matrix and translation
      //TODO: Test
      vnl_matrix<double> rotation =
          quaternion.rotation_matrix_transpose().transpose();
      vnl_vector<double> movedMeanFromCentre = movedMean
          - rotation_centre.GetVnlVector();
      vnl_vector<double> initialMeanFromCentre = initialMean
          - rotation_centre.GetVnlVector();
      vnl_vector<double> translation = movedMeanFromCentre
          - rotation * initialMeanFromCentre;

      //Modulus of Euler angles
      if (this->m_Verbosity >= 3)
      {
        std::cout << "Rotation:\n" << rotation;
        std::cout << "Translation: " << translation << std::endl;
      }

      //Extract parameters
      //Rotation
      for (int ii = 0; ii < 4; ++ii)
        new_parameters[ii] = quaternion[ii];
      for (int ii = 4; ii < 7; ++ii)
        new_parameters[ii] = translation[ii - 4];

      //Get summed change in parameters
      last_delta_param_tol = delta_param_tol;
      delta_param_tol = this->GetParametersDelta(last_parameters, new_parameters);
    } // If Branch for Affine / rigid / QuaternionRigid*****************************

    // End LS
    // LTS ****************************************************

    //Display what is happening
    if (this->m_Verbosity >= 3)
    {
      std::cout << "  Prev Parameters: " << last_parameters;
      std::cout << "\n  New  Parameters: " << new_parameters;
      std::cout << "   Shift: " << delta_param_tol << std::endl;
    }

    for (size_t i=0; i<last_parameters.Size(); i++) {
      last_parameters[i] = new_parameters[i];
    }
  }

  if (this->m_Verbosity >= 3)
  {
    std::cout << "Final optimised transform:\n" << final_transform
        << std::endl;
    std::cout << "Initial transform:\n" << this->m_Transform << std::endl;
  }

  //Now compose - Remember we need to pre multiply
  t12F->Compose(final_transformF, true);
  t12B->Compose(final_transformB, true);

  //Update main transform (for the stopping criterion) -- ignoring
  //any centering operation at this point
  typename TransformType::Pointer invT12B = TransformType::New();
  t12B->GetInverse(invT12B);
  this->m_Transform->SetIdentity();
  this->m_Transform->Compose(invT12B, true);
  this->m_Transform->Compose(t12F, false);
  this->m_Transform->GetParameters();


  if (this->m_Verbosity >= 3)
  {
    std::cout << "# DRH - t12F:\n    translation:"
        << t12F->GetTranslation()
        << "\n    Matrix:\n" << t12F->GetMatrix() << std::endl;
    std::cout << "# DRH - t12B:\n    translation:"
        << t12B->GetTranslation()
        << "\n    Matrix:\n" << t12B->GetMatrix() << std::endl;
    std::cout << "#DRH - Final optimised transform (composed):\n"
        << this->m_Transform << std::endl;
  }
}

} //namespace itk

#endif
