/*=========================================================================
Program: mirorr
Module: itkSymMirorrRegistrationMethod.h
Author: David Rivest-Henault (DRH, riv019)
Created: 26 Oct 2012

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#ifndef __symMirorrRegistrationMethod_h_
#define __symMirorrRegistrationMethod_h_

#include "itkMirorrRegistrationMethod.h"

namespace itk
{

template <typename TMovingImage, typename TFixedImage>
class ITK_EXPORT SymMirorrRegistrationMethod :
  public MirorrRegistrationMethod<TMovingImage, TFixedImage>
{
public:
  /** Standard class typedefs. */
  typedef SymMirorrRegistrationMethod Self;
  typedef MirorrRegistrationMethod<TMovingImage, TFixedImage> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SymMirorrRegistrationMethod, MirorrRegistrationMethod);

  /**  Type of the Moving image. */
  typedef TMovingImage MovingImageType;
  typedef typename MovingImageType::ConstPointer MovingImageConstPointer;
  typedef typename MovingImageType::RegionType MovingImageRegionType;

  /**  Type of the Fixed image. */
  typedef TFixedImage FixedImageType;
  typedef typename FixedImageType::ConstPointer                 FixedImageConstPointer;

  /** Typedef from superclass */
  typedef typename Superclass::MaskImageType                    MaskImageType;
  typedef typename Superclass::MaskImagePointer                 MaskImagePointer;

  typedef typename Superclass::TransformType                    TransformType;
  typedef typename Superclass::TransformPointer                 TransformPointer;

  typedef typename Superclass::TransformOutputType              TransformOutputType;
  typedef typename Superclass::TransformOutputPointer           ransformOutputPointer;
  typedef typename Superclass::TransformOutputConstPointer      TransformOutputConstPointer;

  typedef typename Superclass::InterpolatorType                 InterpolatorType;
  typedef typename Superclass::InterpolatorPointer              InterpolatorPointer;
  typedef typename Superclass::MaskInterpolatorType             MaskInterpolatorType;
  typedef typename Superclass::MaskInterpolatorPointer          MaskInterpolatorPointer;

  typedef typename Superclass::ParametersType                   ParametersType;

  typedef typename Superclass::BlockMatcherType                 BlockMatcherType;
  typedef typename Superclass::BlockMatcherPointer              BlockMatcherPointer;

  typedef typename Superclass::PointListType                    PointListType;
  typedef typename Superclass::PointType                        PointType;
  typedef typename Superclass::VectorType                       VectorType;

  /** Set/Get the Interpolator. */
  itkSetObjectMacro( MovingImageInterpolator, InterpolatorType );
  itkGetObjectMacro( MovingImageInterpolator, InterpolatorType );
  itkSetObjectMacro( MovingMaskInterpolator, MaskInterpolatorType );
  itkGetObjectMacro( MovingMaskInterpolator, MaskInterpolatorType );
  itkSetObjectMacro( FixedImageInterpolator, InterpolatorType );
  itkGetObjectMacro( FixedImageInterpolator, InterpolatorType );
  itkSetObjectMacro( FixedMaskInterpolator, MaskInterpolatorType );
  itkGetObjectMacro( FixedMaskInterpolator, MaskInterpolatorType );

  /** Method to return the latest modified time of this object or
    * any of its cached ivars */
  unsigned long GetMTime() const;

  /** Initialize by setting the interconnects between the components. */
  virtual void Initialize() throw (ExceptionObject);

protected:
  InterpolatorPointer      &m_MovingImageInterpolator;
  MaskInterpolatorPointer  &m_MovingMaskInterpolator;

  InterpolatorPointer      m_FixedImageInterpolator;
  MaskInterpolatorPointer  m_FixedMaskInterpolator;

  SymMirorrRegistrationMethod();
  virtual ~SymMirorrRegistrationMethod();
  virtual void PrintSelf(std::ostream& os, Indent indent) const;

  virtual void PrintMeanDisplacement(PointType p1, PointType p2) const;
  virtual void PrintMeanDisplacement(PointListType pl1, PointListType pl2) const;

  virtual PointType ComputeImageCenter(MovingImageConstPointer img);
  virtual void ComputeHalfSpaceTransforms(const TransformPointer t,
      TransformPointer t12, TransformPointer it12) const;

  virtual void StartOptimization(void);

  virtual void OptimizeTransform() { /* Not to be implemented */ };
  virtual void OptimizeTransform(TransformPointer t12F, TransformPointer t12B);

private:
  SymMirorrRegistrationMethod(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSymMirorrRegistrationMethod.txx"
#endif

#endif /* SYMMIRORRREGISTRATIONMETHOD_H_ */
