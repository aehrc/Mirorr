/*=========================================================================
Program: mirorr
Module: itkBlockMatcher.h
Author: Nicholas Dowson
Created: 20 Feb 2009

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
/*
Header file for block matching two images
Todo:
* Make blockgap etc multidimensional
* Generalise to use other metrics
*/

#ifndef __itkBlockMatcher_h
#define __itkBlockMatcher_h

#include "itkAbstractBlockMatcher.h"
#include <itkImage.h>
#include <itkSize.h>
#include <itkImageRegionConstIterator.h>
#include <itkExceptionObject.h>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include "itkMultiThreader.h"
#include "itkImageSource.h"
#include "itkSparseImageRegionConstIterator.h"

#ifdef USE_NPW
  #include "NPWopengl/NPWGLImage.h"
  #ifdef USE_OPENCL
    #include "NPWcuda/npw/NPWCudaHistogramCreator.h"
  #else
    #include "NPWopengl/NPWJointHistogramCreator.h"
  #endif
//#define DB_TIME //Uncomment for timing statistics
#endif


namespace itk
{
/** \class Class for block matching two images.
 * The images must the same type and dimension.
 *
 * Uses: itkSparseImageRegionConstIterator
 */
template < typename ImageType >
class /*ITK_EXPORT*/ BlockMatcher : public AbstractBlockMatcher<ImageType >
{
public:
  /** Standard class typedefs. */
  typedef BlockMatcher                                  Self;
  typedef AbstractBlockMatcher<ImageType >              Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;
  itkStaticConstMacro(ImageDimension, unsigned int,
      ImageType::ImageDimension);

  /** Define the primary outputs of the algorithm */
  typedef Index<ImageDimension>     IndexType;
  typedef std::vector< IndexType >  IndexListType;

  typedef typename ImageType::PointType  PointType;
  typedef std::vector< PointType >       PointListType;

  typedef typename ImageType::SpacingType  SpacingType;
  typedef typename ImageType::SizeType     SizeType;
  typedef typename ImageType::PixelType    PixelType;

  typedef typename Superclass::MaskType    MaskType;
  typedef typename Superclass::MaskPointer MaskPointer;

  /** Define Region Types*/
  typedef itk::ImageRegion<ImageType::ImageDimension>      RegionType;
  typedef itk::SparseImageRegionConstIterator< ImageType > SparseConstIterator;
  typedef itk::ImageRegionConstIterator< ImageType >       ConstIterator;
  typedef itk::ImageRegionConstIterator< MaskType >        MaskConstIterator;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageRegistrationMethod, ProcessObject);

  /** Set/Get the Moving image. */
  void SetBaseImage( const ImageType * movingImage )
    {
      Superclass::SetBaseImage(movingImage);
#ifdef USE_NPW
    if ( this->metricFunc == &itk::BlockMatcher<ImageType>::npwMutualInformation )
      {
          this->setNPWBaseBinSize( movingImage, this->m_NPWbins );
      }
#endif
    }

  /** Set/Get the Fixed image. */
  void SetSearchImage( const ImageType * fixedImage )
    {
      Superclass::SetSearchImage(fixedImage);
#ifdef USE_NPW
    if ( this->metricFunc == &itk::BlockMatcher<ImageType>::npwMutualInformation )
      {
          this->setNPWSearchBinSize( fixedImage, this->m_NPWbins );
      }
#endif
    }

#ifdef USE_NPW
private:
  /*
   * The NPW openGL histogram creator
   */
#ifdef USE_OPENCL
  NPWCudaHistogramCreator * npwHistogramCreator;
#else
  NPWJointHistogramCreator * npwHistogramCreator;
#endif

#endif

public:
  //void SetBlockGap( unsigned int in ) { m_BlockGap = in; };
  /** Set and Get*/
  itkGetMacro( BlockMetricType, std::string );

  virtual void SetBlockMetricType(std::string _arg)
  {
    if( _arg.empty() )
      {
      _arg = "normalized_correlation";
      itkDebugMacro("Empty input. So defaulting normalized_correlation.");
      }
    itkDebugMacro("setting BlockMetricType to " << _arg);
    if (this->m_BlockMetricType!= _arg)
    {
      this->m_BlockMetricType= _arg;
      MetricMapTypeIterator it = metricMap.find( this->m_BlockMetricType );
      if( it == metricMap.end() )
      {
        std::ostringstream ss;
        ss << "ERROR: " << this->m_BlockMetricType
        << " is not a valid blockmetric type. Exiting!" << std::endl;
        ss << "possible options are:" << std::endl;
        for( it = metricMap.begin(); it != metricMap.end(); ++it )
        {
          ss << "  " << it->second << std::endl;
        }
        std::cerr << ss.str() << std::endl;
        throw std::runtime_error(ss.str());
      }
      metricFunc = it->second;
      /*
       * create the histogram creator if it is needed and it doesn't exist
       */
#ifdef USE_NPW
      if ( metricFunc == &itk::BlockMatcher<ImageType>::npwMutualInformation &&
           0 == this->npwHistogramCreator )
      {
#ifdef USE_OPENCL
        npwHistogramCreator = new NPWCudaHistogramCreator();
#else
        npwHistogramCreator = new NPWJointHistogramCreator();
        npwHistogramCreator->setScaleFactor ( this->m_NPWscaleFactor );
        if ( this->m_NPWshapeExpansion )
        {
            npwHistogramCreator->turnShapeExpansionOn();
        }
        else
        {
            npwHistogramCreator->turnShapeExpansionOff();
        }
#endif
      }
#endif
      this->Modified();
    } else {
    }
  }

  virtual void GetOffsets( RegionType region, PointListType &blockPositions,
      PointListType &blockMatchPositions ) const;

  /**This takes a set of block positions in a floating image,
   * computes the match with surrounding blocks in a moving image,
   * returns the offset of the block giving the best match*/
  virtual void GetOffsets( PointListType &blockPositions,
      PointListType &blockMatchPositions );

  /**WRAPPER FOR ALADIN2000
   * This takes a set of block positions in a floating image,
   * computes the match with surrounding blocks in a moving image,
   * returns the offset of the block giving the best match*/
  void GetOffsets_aladin2k( PointListType &blockPositions,
      PointListType &blockMatchPositions );


  //Display key info
  void displayInfo()
  {
    if(m_BlockMetricType.empty())
      {
      std::cout << "No block metric!" << std::endl;
      throw std::runtime_error( "m_BlockMetricType is empty" );
      }

    Superclass::displayInfo();

    std::cout << m_BlockMetricType << std::endl;
  }

protected: ///Attributes

  //The blockmetric
  typedef itk::ImageRegionConstIterator<ImageType> iteratorType;
  typedef double (itk::BlockMatcher<ImageType>::*MetricFun)(iteratorType,
      iteratorType,MaskConstIterator, MaskConstIterator) const;

  /** A map that links blockmetric type strings to functions */
  typedef std::map<std::string, MetricFun> MetricMapType;
  typedef typename MetricMapType::iterator MetricMapTypeIterator;

  MetricMapType metricMap;
  MetricFun metricFunc;
  std::string m_BlockMetricType;

#ifdef USE_NPW
  /*

  //npw parameters
  int  m_NPWbins;
#ifndef USE_OPENCL
  int  m_NPWscaleFactor;
  bool m_NPWshapeExpansion;
#endif
*/
#ifdef DB_TIME
  /*
   * counters to get average time to compute MI
   */
  double imageConvertTime;
  double histogramAcquisitionTime;
  double MIAcquisitionTime;
  int metricCounter;
#endif


  /*
   * set the bin sizes
   */
  void setNPWBaseBinSize( const ImageType * baseImage, const int bins )
  {
      const PixelType * buffer = baseImage->GetBufferPointer();
      int voxelCount  = baseImage->GetLargestPossibleRegion().GetNumberOfPixels();

#ifdef USE_OPENCL
      npwHistogramCreator->setBaseBinSize( buffer, voxelCount, bins );
#else
      float minVal = *std::min_element( buffer, buffer + voxelCount );
      float maxVal = *std::max_element( buffer, buffer + voxelCount );

      npwHistogramCreator->setBaseBinSize( minVal, maxVal, bins );
#endif
  }

  void setNPWSearchBinSize( const ImageType * searchImage, const int bins )
  {
      const PixelType * buffer = searchImage->GetBufferPointer();
      int voxelCount  = searchImage->GetLargestPossibleRegion().GetNumberOfPixels();

#ifdef USE_OPENCL
      npwHistogramCreator->setSearchBinSize( buffer, voxelCount, bins );
#else
      float minVal = *std::min_element( buffer, buffer + voxelCount );
      float maxVal = *std::max_element( buffer, buffer + voxelCount );

      npwHistogramCreator->setSearchBinSize( minVal, maxVal, bins );
#endif
  }

#endif

  //Constructor is private because classes are generated by object factory
  BlockMatcher() :
    Superclass(),
#ifdef USE_NPW
    npwHistogramCreator ( 0 ),
#endif
    metricFunc( NULL ),
    m_BlockMetricType ( "" ), //Default is below
#ifdef DB_TIME
    imageConvertTime ( 0.0 ),
    histogramAcquisitionTime ( 0.0 ),
    MIAcquisitionTime ( 0.0 ),
    metricCounter ( 0 ),
#endif
/*#ifdef USE_NPW
    m_NPWbins( 32 ),
#ifndef USE_OPENCL
    m_NPWscaleFactor( 2 ),
    m_NPWshapeExpansion( false ),
#endif
#endif*/
    scale_factor( 0.1 ), //This deals with high dynamic range images
    block_size( 343 ) //Expected block size
    {
      //Set a similarity measure for the blockmatching
      metricMap["normalized_correlation"]     = &itk::BlockMatcher<ImageType>::normalizedCorrelation;
      metricMap["correlation_ratio"]          = &itk::BlockMatcher<ImageType>::correlationRatio;
      metricMap["sum_of_squared_differences"] = &itk::BlockMatcher<ImageType>::sumOfSquareDifference;
#ifdef USE_NPW
      metricMap["mutual_information"]         = &itk::BlockMatcher<ImageType>::mutualInformation;
      metricMap["npw_mutual_information"]     = &itk::BlockMatcher<ImageType>::npwMutualInformation;
#endif

      m_BlockMetricType = "normalized_correlation";
      metricFunc = metricMap[m_BlockMetricType]; //DEFAULT
    }

  virtual ~BlockMatcher() {;}


  virtual void BeforeGetOffsets();

  /**
   * Filter match lists by removing low-score or ambiguous matches
   */
#ifdef USE_FLOAT_SCORE_BUFFER
  virtual void PostProcessOffsets(const IndexListType &blockPositions, const IndexListType &blockMatchPositions,
      const std::vector<float> &blockMatchScores, const std::vector<float> &blockMatchScores2ndBest,
      PointListType &blockPositionsOut, PointListType &blockMatchPositionsOut) const;
#else
  virtual void PostProcessOffsets(const IndexListType &blockPositions, const IndexListType &blockMatchPositions,
      const std::vector<double> &blockMatchScores, const std::vector<double> &blockMatchScores2ndBest,
      PointListType &blockPositionsOut, PointListType &blockMatchPositionsOut) const;
#endif

  virtual RegionType ComputeRegionForMatching(RegionType baseRegion) const;
  virtual int ComputeNhoodTotal(RegionType matchingRegion) const;

  /**Calculate normalised correlation between two blocks
   * NCC(f,g) = 1/N Sum \frac{f-mean(f)}{g-mean(g)}{ std(f) * std(g) }
   */
  double normalizedCorrelation(
      itk::ImageRegionConstIterator<ImageType> movingIterator,
      itk::ImageRegionConstIterator<ImageType> fixedIterator,
      MaskConstIterator mx_iter, MaskConstIterator my_iter ) const;

  /**
   * A Measure of the statistical dispersion within individual categories
   * 1 - \frac{ sum_{y category} variance(x_{within y category}) }{ variance(x)}
   */
  double correlationRatio(
      itk::ImageRegionConstIterator<ImageType> x_iter,
      itk::ImageRegionConstIterator<ImageType> y_iter,
      MaskConstIterator mx_iter, MaskConstIterator my_iter ) const;
  double scale_factor; //This deals with high dynamic range images
  unsigned int block_size; //Expected block size

  /**
   * Sum of square differences
   */
  double sumOfSquareDifference(
      itk::ImageRegionConstIterator<ImageType> movingIterator,
      itk::ImageRegionConstIterator<ImageType> fixedIterator,
      MaskConstIterator mx_iter, MaskConstIterator my_iter ) const;

    double mutualInformation(
        itk::ImageRegionConstIterator<ImageType> movingIterator,
        itk::ImageRegionConstIterator<ImageType> fixedIterator,
        MaskConstIterator mx_iter, MaskConstIterator my_iter
        ) const;

#ifdef USE_NPW
    /**
     * Conventional Mutual Information
     */
    double getMutualInformationFromHistogram( const NPWGLImage & hist, const bool verbose = false );

    /**
     * Daan: Non-Parametric Windows
     */
    double npwMutualInformation(
        itk::ImageRegionConstIterator<ImageType> movingIterator,
        itk::ImageRegionConstIterator<ImageType> fixedIterator,
        MaskConstIterator mx_iter, MaskConstIterator my_iter
        ) const;
#endif

};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatcher.txx"
#endif

#endif
