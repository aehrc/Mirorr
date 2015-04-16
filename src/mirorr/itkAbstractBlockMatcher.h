/*=========================================================================
Program: mirorr
Module: itkAbstractBlockMatcher.h
Author: David Rivest-Henault
Created: 02 Nov 2012

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
#ifndef __itkAbstractBlockMatcher_h
#define __itkAbstractBlockMatcher_h

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <itkImage.h>
#include <itkSize.h>
#include <itkIndent.h>

/////////////////////////////////////
//
#define ENFORCE_GPU_CPU_CONSISTENCY
//#define DEBUG_GPU_CPU_CONSISTENCY
/////////////////////////////////////

#ifdef ENFORCE_GPU_CPU_CONSISTENCY
#define USE_FLOAT_SCORE_BUFFER
#endif
#define USE_MASK_IN_METRIC

namespace itk
{

  template < typename ImageType >
  class AbstractBlockMatcher : public ProcessObject
  {
  public:

    ///////////////////////////////////////
    //        ITK standards
    /** Standard class typedefs. */
    typedef AbstractBlockMatcher      Self;
    typedef ProcessObject             Superclass;
    typedef SmartPointer<Self>        Pointer;
    typedef SmartPointer<const Self>  ConstPointer;
    itkStaticConstMacro(ImageDimension, unsigned int,
        ImageType::ImageDimension);

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ImageRegistrationMethod, ProcessObject);
    //
    ///////////////////////////////////////

    ///////////////////////////////////////
    //      Types
    typedef typename ImageType::SizeType    SizeType;
    typedef Index<ImageDimension>           IndexType;
    typedef typename ImageType::SpacingType SpacingType;
    typedef std::vector< IndexType >        IndexListType;
    typedef typename ImageType::PointType   PointType;
    typedef std::vector< PointType >        PointListType;
#if ITK_VERSION_MAJOR < 4
    typedef itk::OrientedImage<unsigned char,3> MaskType;
#else
    typedef itk::Image<unsigned char,3> MaskType;
#endif

    typedef typename MaskType::Pointer MaskPointer;
    //
    ///////////////////////////////////////

    ///////////////////////////////////////
    //      Images
    /** Set/Get the Moving image. */
    virtual void SetBaseImage( const ImageType * movingImage );
    itkGetConstObjectMacro( BaseImage, ImageType );
    itkSetObjectMacro(BaseMask,MaskType);
    itkGetObjectMacro(BaseMask,MaskType);
    itkGetConstObjectMacro(BaseMask,MaskType);
    /** Set/Get the Fixed image. */
    virtual void SetSearchImage( const ImageType * fixedImage );
    itkGetConstObjectMacro( SearchImage, ImageType );
    itkSetObjectMacro(SearchMask,MaskType);
    itkGetObjectMacro(SearchMask,MaskType);
    itkGetConstObjectMacro(SearchMask,MaskType);
    //
    ///////////////////////////////////////

    ///////////////////////////////////////
    //      Algorithm Parameters
    itkGetMacro( BlockGap,   	unsigned int );
    itkSetMacro( BlockGap,   	unsigned int );
    itkGetMacro( BlockWidth, 	unsigned int );
    itkSetMacro( BlockWidth, 	unsigned int );
    itkGetMacro( NhoodGap,   	unsigned int );
    itkSetMacro( NhoodGap,   	unsigned int );
    itkGetMacro( NhoodWidth, 	unsigned int );
    itkSetMacro( NhoodWidth, 	unsigned int );
    itkGetMacro( Padding, 	  unsigned int );
    void SetPadding( unsigned int in );
    itkGetMacro( Verbosity,  	int );
    itkSetMacro( Verbosity,  	int );
    itkGetMacro( PortionMatchesKept, double);
#ifdef USE_NPW
    itkGetMacro( NPWbins,  int );
    itkSetMacro( NPWbins,  int );
#ifndef USE_OPENCL
    itkGetMacro( NPWscaleFactor,  int );
    itkSetMacro( NPWscaleFactor,  int );
    itkGetMacro( NPWshapeExpansion,  bool );
    itkSetMacro( NPWshapeExpansion,  bool );
#endif
#endif

    /** Custom Set method for PortionMatchesKept, which ensures values
     * are kept between 0 and 1*/
    void SetPortionMatchesKept( double in );
    //
    ///////////////////////////////////////

    ///////////////////////////////////////
    //      Allows to perform block matching in the other direction
    // Swaps m_BaseImage <-> m_SearchImage,
    //   and m_BaseMask  <-> m_SearchMask.
    virtual void SwapImageOrder();

    ///////////////////////////////////////
    //      Match function
    virtual void GetOffsets(    PointListType &blockPositions,
                                PointListType &blockMatchPositions )
    {
      if(m_Verbosity>=1)
        std::cout<<"Usage of the empty GetOffsets method. Overlay it.";
      // to avoid warnings :
      blockPositions.empty(); blockMatchPositions.empty();
    }
    //
    ///////////////////////////////////////

    ///////////////////////////////////////
    //      Displays
    virtual void displayInfo();
    virtual void PrintSelf(std::ostream& os, Indent indent) const;
    //
    ///////////////////////////////////////

    ///////////////////////////////////////
    //      Metrics
    virtual void SetBlockMetricType(std::string _arg)
    {
      if(m_Verbosity>=1)
        std::cout<<"cannot set "<<_arg<<". Usage of the empty SetBlockMetricType method. Overlay it.";
    }
    //
    ///////////////////////////////////////

    double testNormalizedCorrelation(
        itk::ImageRegionConstIterator<ImageType> baseImageIterator,
        itk::ImageRegionConstIterator<ImageType> searchImageIterator) const
    {
      //Some variables to perform normalised correlation between the two blocks
      double count = 0, meanSearch = 0, meanBase = 0;

      //Get local means
      for( baseImageIterator.GoToBegin(), searchImageIterator.GoToBegin();
          !baseImageIterator.IsAtEnd();
          ++baseImageIterator, ++searchImageIterator, ++count)
      {
        meanBase += baseImageIterator.Get();
        meanSearch += searchImageIterator.Get();
      }
      meanBase /= count;
      meanSearch /= count;

      double baseSearchProd = 0.0;
      double baseSumSqr = 0.0;
      double searchSumSqr = 0.0;
      for( baseImageIterator.GoToBegin(), searchImageIterator.GoToBegin();
          !baseImageIterator.IsAtEnd();
          ++baseImageIterator, ++searchImageIterator)
      {
        double baseBary = baseImageIterator.Get() - meanBase;
        double searchBary = searchImageIterator.Get() - meanSearch;

        baseSearchProd += baseBary * searchBary;
        baseSumSqr += baseBary * baseBary;
        searchSumSqr += searchBary * searchBary;
      }

      if( baseSumSqr<=0.0 || searchSumSqr<=0.0 )
        return -2.0;

      return baseSearchProd*baseSearchProd / baseSumSqr / searchSumSqr;
    }

    double testNormalizedCorrelation(IndexType baseIdx, IndexType searchIdx)
    {
      typename ImageType::RegionType::SizeType size; size[0] = m_BlockWidth; size[1] = m_BlockWidth; size[2] = m_BlockWidth;
      typename ImageType::RegionType region;
      region.SetSize(size);

      region.SetIndex(baseIdx);
      typename itk::ImageRegionConstIterator<ImageType> baseImageIterator(this->m_BaseImage, region);

      region.SetIndex(searchIdx);
      typename itk::ImageRegionConstIterator<ImageType> searchImageIterator(this->m_SearchImage, region);


      return testNormalizedCorrelation(baseImageIterator, searchImageIterator);
    }

    double testNormalizedCorrelation(PointType basePt, PointType searchPt)
    {
      IndexType baseIdx, searchIdx;
      this->m_BaseImage->TransformPhysicalPointToIndex(basePt, baseIdx);
      this->m_SearchImage->TransformPhysicalPointToIndex(searchPt, searchIdx);
      return testNormalizedCorrelation(baseIdx, searchIdx);
    }
  protected:

    // Algorithm parameters
    //The gap between blocks within a nhood
    unsigned int m_BlockGap;   //In pixels
    //The size of blocks
    unsigned int m_BlockWidth; //In Pixels
    //The gap between nhoods
    unsigned int m_NhoodGap;   //In Pixels
    //How many blocks around an nhood should be matched
    unsigned int m_NhoodWidth; //In Blocks.
    //How much info should the algorithm display as it runs
    int           m_Verbosity;
    //Each nhood has one block which gives the best match.
    //How many of these should be kept?
    double        m_PortionMatchesKept; //How many best matches should be kept
    unsigned int  m_Padding;

#ifdef USE_NPW
    //npw parameters
    int  m_NPWbins;
#ifndef USE_OPENCL
    int  m_NPWscaleFactor;
    bool m_NPWshapeExpansion;
#endif
#endif

    typename ImageType::ConstPointer m_BaseImage;
    typename ImageType::ConstPointer m_SearchImage;
    typename MaskType::Pointer m_BaseMask;
    typename MaskType::Pointer m_SearchMask;


    ///////////////////////////////////////
    //      Constructors
    AbstractBlockMatcher();
    virtual ~AbstractBlockMatcher() {}
    //
    ///////////////////////////////////////

  public:
    //!Convert a block position to a point in mm
    PointType blockPositionToPoint( IndexType index ) const
    {
      PointType point;
      //Use Base image - not search image

      //Get position of corner
      this->m_BaseImage->TransformIndexToPhysicalPoint( index, point );

      //Calculate block centre
      double halfBlockSizePixels = 0.5 * this->m_BlockWidth;
      SpacingType halfBlockSizeMm =
          this->m_BaseImage->GetSpacing() * halfBlockSizePixels;

      point += halfBlockSizeMm;

      return point;
    }

    //!Inverse of the above. Convert a point in mm to block position
    IndexType pointToBlockPosition( PointType point ) const
    {
      IndexType index;
      //Use Base image - not search image

      //Calculate block centre
      double halfBlockSizePixels = 0.5 * this->m_BlockWidth;
      SpacingType halfBlockSizeMm =
          this->m_BaseImage->GetSpacing() * halfBlockSizePixels;

      point -= halfBlockSizeMm;

      //Get position of corner
      //this->m_BaseImage->TransformIndexToPhysicalPoint( index, point );
      this->m_BaseImage->TransformPhysicalPointToIndex( point, index );

      return index;
    }

#ifdef DEBUG_GPU_CPU_CONSISTENCY
    mutable std::vector<float> scoresOut;
    mutable std::vector<float> scoresOut2nd;
#endif
  };

}


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAbstractBlockMatcher.txx"
#endif


#endif
