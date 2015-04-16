/*=========================================================================
Program: mirorr
Module: itkMirorrRegistrationMethod.h
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
 * Header file for aligning two images using the Mirorr algorithm
 */

#ifndef __itkMirorrRegistrationMethod_h
#define __itkMirorrRegistrationMethod_h

#include <itkProcessObject.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_quaternion.h>
#include "itkAbstractBlockMatcher.h"

namespace itk
{

/**
 * TODO: Create Special TranslationTransform type with SetIdentity()
 */


/**
 * \class Class implementing the Mirorr registration method. It has been
 * copied directly from itkImageRegistrationMethod, with some tweaks,
 * because it does not use a metric based on the similarity between
 * two images. Instead it locates nearby blocks in the fixed image
 * (after transformation and interpolation) that match the
 * corresponding block in the moving image. The 75% of block matches
 * with low scores are thrown away. Then the metric is used to
 * minimise the discrepancy between the global transformation and the
 * transform suggested by each block match (voting).  As with
 * itkImageRegistrationMethod, the Optimizer, Transform, Interpolator,
 * Metric and Images need to be set before the algorithm is run.
 *
 * Used by: MirorrImplement
 * Uses: itkBlockMatcher
 *
 * To do:
 * - Ideally the class would inherit directly from
 *   itkImageRegistrationMethod, but this proved to be problematic. This
 *   on the to do list.
 * - Allow the block matching metric to be selected.
 * - Stop iterating when the largest voxcel movement is less than 0.2
 *   voxels.
 */
template <typename TMovingImage, typename TFixedImage>
class ITK_EXPORT MirorrRegistrationMethod : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef MirorrRegistrationMethod  Self;
  typedef ProcessObject              Superclass;
  typedef SmartPointer<Self>         Pointer;
  typedef SmartPointer<const Self>   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageRegistrationMethod, ProcessObject);

  /**  Type of the Moving image. */
  typedef          TMovingImage                  MovingImageType;
  typedef typename MovingImageType::ConstPointer MovingImageConstPointer;

  /**  Type of the Fixed image. */
  typedef          TFixedImage                   FixedImageType;
  typedef typename FixedImageType::ConstPointer  FixedImageConstPointer;

#if ITK_VERSION_MAJOR < 4
  typedef itk::OrientedImage<unsigned char,3>    MaskImageType;
  typedef int ThreadIdType;
#else
  typedef itk::Image<unsigned char,3>    MaskImageType;
#endif
  
  typedef typename MaskImageType::Pointer        MaskImagePointer;

  typedef typename TMovingImage::RegionType MovingImageRegionType;

  /**  Type of the Transform . */
  //typedef  typename MetricType::TransformType      TransformType;
  typedef MatrixOffsetTransformBase< double, 3, 3 >  TransformType;
  typedef  typename TransformType::Pointer           TransformPointer;

  /** Type for the output: Using Decorator pattern for enabling
   *  the Transform to be passed in the data pipeline */
  typedef  DataObjectDecorator< TransformType >       TransformOutputType;
  typedef typename TransformOutputType::Pointer       TransformOutputPointer;
  typedef typename TransformOutputType::ConstPointer  TransformOutputConstPointer;

  /**  Type of the Interpolator. */
  //typedef  typename MetricType::InterpolatorType   InterpolatorType;
  typedef InterpolateImageFunction< FixedImageType, double > InterpolatorType;
  typedef typename InterpolatorType::Pointer                 InterpolatorPointer;
  typedef InterpolateImageFunction< MaskImageType, double >  MaskInterpolatorType;
  typedef typename MaskInterpolatorType::Pointer             MaskInterpolatorPointer;

  /** Type of the Transformation parameters This is the same type used to
   *  represent the search space of the optimization algorithm */
  typedef typename TransformType::ParametersType    ParametersType;

  /** Smart Pointer type to a DataObject. */
  typedef typename DataObject::Pointer DataObjectPointer;

  /** Type for block matching */
  typedef AbstractBlockMatcher<TMovingImage >           BlockMatcherType;
  typedef typename BlockMatcherType::Pointer            BlockMatcherPointer;
  
  typedef typename BlockMatcherType::PointListType      PointListType;
  typedef typename BlockMatcherType::PointType          PointType;
  typedef Vector<double,FixedImageType::ImageDimension> VectorType;

  // Define the inner resampling space
  enum EResamplingType {EResamplingBasic=0, EResamplingMiddle, EResamplingMoving, EResamplingFixed};

  /** Method that initiates the registration. This will Initialize and ensure
   * that all inputs the registration needs are in place, via a call to
   * Initialize() will then start the optimization process via a call to
   * StartOptimization()  */
  void StartRegistration(void);

  /** Method that initiates the optimization process. This method should not be
   * called directly by the users. Instead, this method is intended to be
   * invoked internally by the StartRegistration() which is in turn invoked by
   * the Update() method.
   * FIXME: This method should be declared protected. */
  virtual void StartOptimization(void);

  /** Set/Get the Moving image. */
  void SetMovingImage( const MovingImageType * movingImage );
  itkGetConstObjectMacro( MovingImage, MovingImageType );

  /** Set/Get the Fixed image. */
  void SetFixedImage( const FixedImageType * fixedImage );
  itkGetConstObjectMacro( FixedImage, FixedImageType );

#ifdef USE_OPENCL
  /** Set/Get UseGpuOn option */
  itkSetMacro( UseGpuOn, bool );
  itkGetMacro( UseGpuOn, bool );
#endif

  itkSetMacro( ResamplingType, EResamplingType );
  itkGetMacro( ResamplingType, EResamplingType );

  itkSetMacro( UseMultiThreading, bool );
  itkBooleanMacro( UseMultiThreading );

  /** Set/Get the Transfrom. */
  itkSetObjectMacro( Transform, TransformType );
  itkGetObjectMacro( Transform, TransformType );

  /** Set/Get the Interpolator. */
  itkSetObjectMacro( Interpolator, InterpolatorType );
  itkGetObjectMacro( Interpolator, InterpolatorType );
  itkSetObjectMacro( MaskInterpolator, MaskInterpolatorType );
  itkGetObjectMacro( MaskInterpolator, MaskInterpolatorType );

  /** Set/Get the Masks. */
  itkSetObjectMacro( FixedMask, MaskImageType );
  itkGetObjectMacro( FixedMask, MaskImageType );
  itkSetObjectMacro( MovingMask, MaskImageType );
  itkGetObjectMacro( MovingMask, MaskImageType );

  /** Set/Get the Metric. */
  itkSetMacro( Verbosity, int );
  itkGetMacro( Verbosity, int );

  /** Set/Get the Max No. iterations. */
  itkSetMacro( MaxIterations, int );
  itkGetMacro( MaxIterations, int );

  //How many block matches to keep
  itkSetMacro( PortionMatchesKept, double );
  itkGetMacro( PortionMatchesKept, double );

  /** Set/Get the string specifying the block matcher metric */
  itkGetMacro( BlockMetricType, std::string );
  itkSetMacro( BlockMetricType, std::string );

  /** Set / Get "deep" algorithm parameters */
  itkGetMacro( NhoodWidth, int );
  itkGetMacro( NhoodGap, int );
  itkGetMacro( BlockWidth, int );
  itkGetMacro( BlockGap, int );
#ifdef USE_NPW
  itkGetMacro( NPWbins, int );
#ifndef USE_OPENCL
  itkGetMacro( NPWscaleFactor, int );
  itkGetMacro( NPWshapeExpansion, bool );
#endif
#endif

  itkSetMacro( NhoodWidth, int );
  itkSetMacro( NhoodGap, int );
  itkSetMacro( BlockWidth, int );
  itkSetMacro( BlockGap, int );
#ifdef USE_NPW
  itkSetMacro( NPWbins, int );
#ifndef USE_OPENCL
  itkSetMacro( NPWscaleFactor, int );
  itkSetMacro( NPWshapeExpansion, bool );
#endif
#endif

  /** Get the last transformation parameters visited by
   * the optimizer. */
  itkGetConstReferenceMacro( LastTransformParameters, ParametersType );

  /** Set the region of the moving image to be considered as region of
     interest during the registration. This region will be passed to
     the ImageMetric in order to restrict the metric computation to
     consider only this region.
     \warning The same region can also be set directly into the metric.
     please avoid to set the region in both places since this can lead
     to inconsistent configurations.  */
  void SetMovingImageRegion( const  MovingImageRegionType & region );
  /** Get the region of the moving image to be considered as region of
     interest during the registration. This region will be passed to
     the ImageMetric in order to restrict the metric computation to
     consider only this region.  */
  itkGetConstReferenceMacro( MovingImageRegion, MovingImageRegionType );
  /** True if a region has been defined for the moving image to which
     the ImageMetric will limit its computation */
  itkGetMacro( MovingImageRegionDefined, bool );
  /** Turn on/off the use of a moving image region to which
     the ImageMetric will limit its computation.
     \warning The region must have been previously defined using the
     SetMovingImageRegion member function */
  itkSetMacro( MovingImageRegionDefined, bool );

  /** Initialize by setting the interconnects between the components. */
  virtual void Initialize() throw (ExceptionObject);

  /** Block-matcher object management */
  virtual void CreateBlockMatcher();
  virtual void CreateBlockMatcherCPU();
#ifdef USE_OPENCL
  virtual void CreateBlockMatcherGPU();
#endif
  virtual void SetupBlockMatcher();
  virtual void SetupBlockMatcher_tryNC();

  /** Returns the transform resulting from the registration process  */
  const TransformOutputType * GetOutput() const;

  /** Make a DataObject of the correct type to be used as the specified
   * output. */
#if ITK_VERSION_MAJOR < 4
  virtual DataObjectPointer MakeOutput(unsigned int idx);
#else
  virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx);
  virtual DataObjectPointer MakeOutput( const DataObjectIdentifierType & ) { /* Place holder */ return 0;}
#endif

  /** Method to return the latest modified time of this object or
   * any of its cached ivars */
  virtual unsigned long GetMTime() const;

  itkSetMacro( NumberOfBlockMatcherThreads, ThreadIdType );
  itkGetMacro( NumberOfBlockMatcherThreads, ThreadIdType );

  //NumberOfBlockMatcherThreads

protected:
  MirorrRegistrationMethod();
  virtual ~MirorrRegistrationMethod() {};
  virtual void PrintSelf(std::ostream& os, Indent indent) const;

  /** Method invoked by the pipeline in order to trigger the computation of
   * the registration. */
  virtual void  GenerateData ();

  PointType ComputeImageCenter(MovingImageConstPointer img);

  /** Provides derived classes with the ability to set this private var */
  itkSetMacro( LastTransformParameters, ParametersType );

private:
  MirorrRegistrationMethod(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  //MetricPointer                    m_Metric; //Used
  //OptimizerType::Pointer           m_Optimizer; //Used

protected:
  FixedImageConstPointer    m_FixedImage; //Used
  MovingImageConstPointer   m_MovingImage; //Used

  MaskImagePointer          m_FixedMask;
  MaskImagePointer          m_MovingMask;

  TransformPointer          m_Transform; //Used
  InterpolatorPointer       m_Interpolator; //Used
  MaskInterpolatorPointer   m_MaskInterpolator; //Used
  //  ParametersType         m_InitialTransformParameters; //Used
  ParametersType            m_LastTransformParameters; //Used

  bool                      m_MovingImageRegionDefined;
  MovingImageRegionType     m_MovingImageRegion;

  bool                      m_UseMultiThreading;

  double                    m_PortionMatchesKept;

  std::string               m_BlockMetricType;

  bool                      m_UseGpuOn;

  //bool                      m_UseSymResampling;
  //bool                      m_UseMaxResampling;
  EResamplingType           m_ResamplingType;

  int m_Verbosity;

  int m_MaxIterations; //The number of iterations to apply

  //"Deep" algorithm parameters
  int m_NhoodWidth;
  int m_NhoodGap;
  int m_BlockWidth;
  int m_BlockGap;
#ifdef USE_NPW
  int  m_NPWbins;
#ifndef USE_OPENCL
  int  m_NPWscaleFactor;
  bool m_NPWshapeExpansion;
#endif
#endif

  //This is NOT set by the user!
  BlockMatcherPointer m_BlockMatcher;

  //Number of thread to be used by the block matcher (if the threaded block matcher is used)
  ThreadIdType m_NumberOfBlockMatcherThreads;

  PointListType       m_InitialPositions;
  PointListType       m_MovedPositions;

  ///Functions for optimisation which use quaternions
  /** Given a set of positions and offsets compute the best affine
   *  transform to make them align
   */
  //ParametersType OptimizeTransform( double & delta_parameters);
  virtual void OptimizeTransform( );

  //Utility functions
  vnl_vector<double> GetBaryCentre( const PointListType &points ) const;
  vnl_matrix<double> CovarianceOfPoints( const PointListType &initialSubset,
      const PointListType &movedSubset ) const;
  vnl_quaternion<double> CalculateRigidQuaternion( const PointListType &initialSubset,
      const PointListType &movedSubset ) const;
  double GetParametersDelta( const ParametersType & last_parameters,
      const ParametersType & final_transform ) const;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMirorrRegistrationMethod.txx"
#endif

#endif
