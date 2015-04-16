/*=========================================================================
Program: mirorr
Module: itkBlockMatcher.txx
Author: Nicholas Dowson
Created: 20 Feb 2009

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
#ifndef __itkBlockMatcher_txx
#define __itkBlockMatcher_txx

#include "itkBlockMatcher.h"
#include <vector>
#include <algorithm>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <stdexcept>
#include <sstream>
#include <iomanip>

#ifdef USE_NPW
#include "NPWcuda/npw/MutualInformation.h"
#include "NPWopengl/NPWStopwatch.h"
#endif

bool go = false;


namespace itk
{

typedef boost::tuple< double, double, int > ScoreIndexType;

bool InvalidScoreIndex( const ScoreIndexType & in ){ return in.get<0>() <= 1e-4; };

//----------------------------------------------------------------------------------------------------

template<typename ImageType>
double
BlockMatcher<ImageType>
::normalizedCorrelation(
    itk::ImageRegionConstIterator<ImageType> baseImageIterator,
    itk::ImageRegionConstIterator<ImageType> searchImageIterator,
    MaskConstIterator mx_iter,
    MaskConstIterator my_iter
) const
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
      ++baseImageIterator, ++searchImageIterator, ++mx_iter, ++my_iter)
  {
#ifdef USE_MASK_IN_METRIC
    if(!mx_iter.Get() || !my_iter.Get()) {
      //continue;
      baseSearchProd += 0.;
      baseSumSqr += 0.;
      searchSumSqr += 0.;
    } else {
#endif

    double baseBary = baseImageIterator.Get() - meanBase;
    double searchBary = searchImageIterator.Get() - meanSearch;

    baseSearchProd += baseBary * searchBary;
    baseSumSqr += baseBary * baseBary;
    searchSumSqr += searchBary * searchBary;
#ifdef USE_MASK_IN_METRIC
    }
#endif
  }

  if( baseSumSqr<=0.0 || searchSumSqr<=0.0 )
    return -2.0;

  return baseSearchProd*baseSearchProd / baseSumSqr / searchSumSqr;
}

//----------------------------------------------------------------------------------------------------

template<typename ImageType>
double
BlockMatcher<ImageType>::correlationRatio(
    itk::ImageRegionConstIterator<ImageType> x_iter,
    itk::ImageRegionConstIterator<ImageType> y_iter,
    MaskConstIterator mx_iter,
    MaskConstIterator my_iter
) const
{
  typedef std::pair<int,double>  ZRT;
  typedef std::vector< ZRT >     VZRT;

  //Copy out x values
  VZRT xy_values;
  xy_values.reserve(block_size);
  for( x_iter.GoToBegin(), y_iter.GoToBegin(); !x_iter.IsAtEnd();
      ++x_iter, ++y_iter, ++mx_iter, ++my_iter )
  {
    if(!mx_iter.Get() || !my_iter.Get())
      continue;

    int yy = static_cast<int>( std::floor( y_iter.Get() * scale_factor ) );
    xy_values.push_back( ZRT( yy, x_iter.Get() ) );
  }

  std::sort( xy_values.begin(), xy_values.end() );

  double sum_x = 0, sum_x2 = 0, cnt_x = 0;
  double sum_xi = 0, sum_xi2 = 0, cnt_xi = 0;
  int yy = xy_values.front().first;
  double sum_var_xi = 0;

  //Loop over pairs of values
  for( VZRT::iterator ii = xy_values.begin(); ii != xy_values.end(); )
  {
    //Increment values
    if( ii->second > 10 )
    {
      sum_xi += ii->second;
      sum_xi2 += ii->second * ii->second;
      cnt_xi++;
    }

    ++ii; //Next pair

    //If there is a new symbol, calculate variance and restart count
    if( ii == xy_values.end() || yy != ii->first )
    {
      //New symbol - if not at end
      if( ii != xy_values.end() )
        yy = ii->first;

      //Only calculate variance if there is a count
      if( cnt_xi > 2 )
      {
        //Calculate variance
        double mean_xi = sum_xi / cnt_xi;
        double var_xi = sum_xi2 / cnt_xi - mean_xi * mean_xi;

        //Add to variance sum
        sum_var_xi += var_xi;
      }

      //Increment large counters
      sum_x += sum_xi;
      sum_x2 += sum_xi2;
      cnt_x += cnt_xi;

      //Reset counts
      sum_xi = sum_xi2 = cnt_xi = 0;
    }
  }

  //Check that some correlatrions did occur
  if( cnt_x > 2 && sum_var_xi > 0 )
  {
    //Get total variance
    double mean_x = sum_x / cnt_x;
    double var_x = sum_x2 / cnt_x - mean_x * mean_x;

    //Return correlation ratio
    double cr = 1.0 - sum_var_xi / var_x;
    return cr > .4 ? cr : -2;
  }
  else
    return -2;
}

//----------------------------------------------------------------------------------------------------

template<typename ImageType>
double
BlockMatcher<ImageType>
::sumOfSquareDifference(
    itk::ImageRegionConstIterator<ImageType> baseImageIterator,
    itk::ImageRegionConstIterator<ImageType> searchImageIterator,
    MaskConstIterator mx_iter,
    MaskConstIterator my_iter
) const
{
  double sum_square_differences = 0;
  double count = 0;

  for( baseImageIterator.GoToBegin(), searchImageIterator.GoToBegin();
      !baseImageIterator.IsAtEnd(); ++baseImageIterator, ++searchImageIterator,
      ++mx_iter, ++my_iter )
  {
    if(!mx_iter.Get() || !my_iter.Get())
      continue;

    double difference =  baseImageIterator.Get() - searchImageIterator.Get();
    sum_square_differences += difference * difference;
    ++count;
  }
  if( count == 0 )
  {
    return -2;
  }
  else
  {
    return sum_square_differences / count;
  }
}

//----------------------------------------------------------------------------------------------------
#ifdef USE_NPW
template<typename ImageType>
double BlockMatcher<ImageType>::mutualInformation(
    itk::ImageRegionConstIterator<ImageType> baseImageIterator,
    itk::ImageRegionConstIterator<ImageType> searchImageIterator,
    MaskConstIterator mx_iter,
    MaskConstIterator my_iter
) const
  {

  NPWStopwatch sw;

  const int bins = 64;
  float * hist = new float [bins*bins];
  for (int i = 0; i < bins*bins; ++i)
  {
    hist[i] = 0.0f;
  }

  /*
   * get binsizes
   */
  baseImageIterator.GoToBegin();
  searchImageIterator.GoToBegin();
  PixelType min1 = baseImageIterator.Get(); // 1 = base
  PixelType max1 = min1;
  PixelType min2 = searchImageIterator.Get(); // 2 = search
  PixelType max2 = min2;
  unsigned int count = 0;
  while ( !baseImageIterator.IsAtEnd() )
  {
    PixelType val1 = baseImageIterator.Get();
    PixelType val2 = searchImageIterator.Get();
    min1 = std::min(val1, min1);
    max1 = std::max(val1, max1);
    min2 = std::min(val2, min2);
    max2 = std::max(val2, max2);
    ++baseImageIterator;
    ++searchImageIterator;
    ++count;
  }

  if ( min1 == max1 || min2 == max2 )
  {
    // one image only contains one intensity, MI = 0.0, so return -2.0
    if (hist)
    {
      delete[] hist;
      hist = 0;
    }
    return -2.0;
  }

  double binsize1 = static_cast< double > (max1 - min1) / static_cast< double > (bins);
  double binsize2 = static_cast< double > (max2 - min2) / static_cast< double > (bins);

  /*
   * fill histogram
   */
  double weight = 1.0 / static_cast< double > (count);
  for ( baseImageIterator.GoToBegin(), searchImageIterator.GoToBegin();
      !baseImageIterator.IsAtEnd();
      ++baseImageIterator, ++searchImageIterator, ++mx_iter, ++my_iter )
  {
    if(!mx_iter.Get() || !my_iter.Get())
      continue;

    int x = static_cast< int > ( ( baseImageIterator.Get()   - min1 ) / binsize1);
    int y = static_cast< int > ( ( searchImageIterator.Get() - min2 ) / binsize2);
    x = std::max( 0, std::min(x, bins - 1) );
    y = std::max( 0, std::min(y, bins - 1) );

    hist[x + y*bins] += (float)weight;
  }

  /*
   * get Mutual Information from histogram
   */

  double MI = NPW::getMutualInformationFromHistogram(hist,bins);

  if (hist)
  {
    delete[] hist;
    hist = 0;
  }

  //std::cout << "\tmi took " << sw.lap() << " s" << std::endl;

  return ( 0.0 == MI ? -2.0 : MI );
  }
#endif // USE_NPW

#ifdef USE_NPW
#ifdef USE_OPENCL
//---------------------------------------------------------------------------
/*
 * implement NPW mutual information with OpenGL & CUDA
 */
template<typename ImageType>
double BlockMatcher<ImageType>::npwMutualInformation(
    itk::ImageRegionConstIterator<ImageType> baseImageIterator,
    itk::ImageRegionConstIterator<ImageType> searchImageIterator,
    MaskConstIterator /*mx_iter*/,
    MaskConstIterator /*my_iter*/
    ) const
    {
  //std::cout << "BlockMatcher::npwMutualInformation() -- Using the CUDA implementation" << std::endl;
  NPWStopwatch sw;

  /*
   * get ITK iterators as float arrays and get the dimensions
   */
  itk::Size<3u> itkSize = baseImageIterator.GetRegion().GetSize();

  int voxelCount = itkSize[0] * itkSize[1] * itkSize[2];

  PixelType * baseData   = new PixelType [ voxelCount ];
  PixelType * searchData = new PixelType [ voxelCount ];

  if ( !baseData || !searchData )
  {
    std::cout << "BlockMatcher::npwMutualInformation() null pointer buffer!" << std::endl;
    return -2.0;
  }


  baseImageIterator.GoToBegin();
  for (int i = 0; !baseImageIterator.IsAtEnd(); ++i)
  {
    baseData[i] = baseImageIterator.Get();
    ++baseImageIterator;
  }
  searchImageIterator.GoToBegin();
  for (int i = 0; !searchImageIterator.IsAtEnd(); ++i)
  {
    searchData[i] = searchImageIterator.Get();
    ++searchImageIterator;
  }
  /*
   * get the MI
   */
  float MI = npwHistogramCreator->getMutualInformation(
      baseData, searchData,
      itkSize[0],itkSize[1],itkSize[2] );

  //std::cout << "BlockMatcher::npwMutualInformation() got MI: " << MI << std::endl;
  /*
   * clean up
   */
  if (baseData)   { delete [] baseData;   baseData = 0; }
  if (searchData) { delete [] searchData; searchData = 0; }

  //std::cout << "\tnpwmi took " << sw.lap() << " s" << std::endl;

  /*
   * and return
   */
  return ( MI == 0.0f ? -2.0f : MI );

    }

#else // USE_OPENCL (no GPU implementation below)

//---------------------------------------------------------------------------
/*
 * implement NPW mutual information with openGL (but not CUDA)
 */
template<typename ImageType>
double BlockMatcher<ImageType>::npwMutualInformation(
    itk::ImageRegionConstIterator<ImageType> baseImageIterator,
    itk::ImageRegionConstIterator<ImageType> searchImageIterator,
    MaskConstIterator mx_iter,
    MaskConstIterator my_iter
) const
  {
  //std::cout << "BlockMatcher::npwMutualInformation() -- Using the non-CUDA implementation" << std::endl;
#ifdef DB_TIME
  NPWStopwatch imageConvertStopwatch;
#endif
  /*
   * get images from iterators
   */
  itk::Size<3u> itkSize = baseImageIterator.GetRegion().GetSize();

  NPWGLImage baseImage  ( itkSize[0], itkSize[1], itkSize[2] );
  NPWGLImage searchImage( itkSize[0], itkSize[1], itkSize[2] );

  int i = 0;
  for( baseImageIterator.GoToBegin(), searchImageIterator.GoToBegin();
      ! baseImageIterator.IsAtEnd();
      ++baseImageIterator, ++searchImageIterator, ++mx_iter, ++my_iter)
  {
    if(!mx_iter.Get() || !my_iter.Get())
      continue;

    float baseVal   = static_cast<float>( baseImageIterator.Get()   );
    float searchVal = static_cast<float>( searchImageIterator.Get() );

    baseImage.SetValueAtIndex  (i, baseVal );
    searchImage.SetValueAtIndex(i, searchVal );
    ++i;
  }

#ifdef DB_TIME
  imageConvertTime += imageConvertStopwatch.lap();
#endif

  /*
   * get histogram
   */
#ifdef DB_TIME
  NPWStopwatch histogramAcquisitionStopwatch;
#endif
  NPWGLImage hist( this->m_NPWbins, this->m_NPWbins );

  // Order is important here! Base image first, search image second
  int result = npwHistogramCreator->getHistogram( &baseImage, &searchImage, &hist );
#ifdef DB_TIME
  histogramAcquisitionTime += histogramAcquisitionStopwatch.lap();
#endif
  if ( result != 1 )
  {
    std::cout << "*** ERROR *** : could not correctly render joint histogram" << std::endl;
    return -2.0;
  }
  else
  {
    /*
     * get mutual information
     */
#ifdef DB_TIME
    NPWStopwatch MIAcquisitionStopwatch;
#endif
    //TODO
    double MI = NPW::getMutualInformationFromHistogram( hist.GetData(), this->m_NPWbins );
#ifdef DB_TIME
    MIAcquisitionTime += MIAcquisitionStopwatch.lap();
#endif

#ifdef DB_TIME
    ++metricCounter;
    timeCounter += (double)stopwatch.lap();

    int metricLimit = 53166;
    if ( metricCounter > metricLimit )
    {
      std::cout << "result for " << metricLimit
          << " (" <<  searchImage.dim(0) << "x"
          <<  searchImage.dim(1) << "x"
          <<  searchImage.dim(2) << ") blocks" << std::endl;
      npwHistogramCreator->getRenderer()->printCounters();
      std::cout << "total time spent in npw blockmetric : " << timeCounter << " s" << std::endl;

      std::cout << "image converting : \ttotal time = " << imageConvertTime << "\tavg time = "
          << imageConvertTime/(double)metricCounter << std::endl;

      std::cout << "histogram acquisition : \ttotal time = " << histogramAcquisitionTime << "\tavg time = "
          << histogramAcquisitionTime/(double)metricCounter << std::endl;

      std::cout << "MI : \ttotal time = " << MIAcquisitionTime
          << "\tavg time = " << MIAcquisitionTime/(double)metricCounter << std::endl;

      delete npwHistogramCreator;
      throw std::runtime_error("blablabla");
    }
#endif

    return ( MI == 0.0 ? -2.0 : MI );
  }
    }

#endif // USE_OPENCL

#endif // USE_NPW
//---------------------------------------------------------------------------

template < typename ImageType >
void BlockMatcher<ImageType>
::BeforeGetOffsets()
{
  RegionType baseRegion = this->m_BaseImage->GetLargestPossibleRegion();

  if(!this->GetBaseMask())
  {
    MaskPointer mask = MaskType::New();
    mask->SetRegions( baseRegion );
    mask->Allocate();
    mask->FillBuffer(255);
    this->SetBaseMask( mask );
  }
  if(!this->GetSearchMask())
  {
    MaskPointer mask = MaskType::New();
    mask->SetRegions( baseRegion );
    mask->Allocate();
    mask->FillBuffer(255);
    this->SetSearchMask( mask );
  }
}

template < typename ImageType >
typename BlockMatcher<ImageType>::RegionType
BlockMatcher<ImageType>
::ComputeRegionForMatching(RegionType baseRegion) const
{
  int nNhoodsTotal = 1;
  SizeType nNhoods;
  //Exclude fractional region at end
  //Remaining region includes STARTING points of regions only
  {
    IndexType baseIndex = baseRegion.GetIndex();
    SizeType baseRegionSize = baseRegion.GetSize();
    for (int ii = 0; ii < ImageType::ImageDimension; ++ii)
    {
      nNhoods[ii] = (baseRegionSize[ii] - this->m_BlockWidth)
          / this->m_NhoodGap + 1;
      int tRegionSize = this->m_NhoodGap * (nNhoods[ii] - 1)
          + this->m_BlockWidth;
      int tGapAtEnd = baseRegionSize[ii] - tRegionSize;
      if (tGapAtEnd < 0)
      {
        std::stringstream ss;
        ss << "ERROR: Negative Gap at End - this should be impossible."
            << " RegionSize: " << baseRegionSize[ii] << " BlockWidth: "
            << this->m_BlockWidth << " Nhood Gap: " << this->m_NhoodGap
            << " No. Nhoods: " << nNhoods[ii] << " New RegionSize: "
            << tRegionSize << " GapAtEnd: " << tGapAtEnd << std::endl;
        throw std::runtime_error(ss.str());
      }
      //Size of the base image region containing blocks
      baseRegionSize[ii] = tRegionSize;
      //Centre the region in the image to avoid bias away from centre
      baseIndex[ii] += tGapAtEnd / 2;

      //Because Sparse iterator defines top left front corner and does
      //not know about length of the block remove an entire block less 1 pixel
      baseRegionSize[ii] -= this->m_BlockWidth - 1;

      //baseRegionSize[ii] = baseRegionSize[ii] -
      //std::max(this->m_BlockWidth,this->m_BlockGap) - this->m_BlockWidth/2;
      //baseIndex[ii] += this->m_BlockWidth/2; //Ensure blocks are centred

      //Accumulate number of Neighbourhoods
      nNhoodsTotal *= nNhoods[ii];
    }

    baseRegion.SetSize(baseRegionSize);
    baseRegion.SetIndex(baseIndex);
  }

  //Display some info
  if( this->m_Verbosity >= 2 )
  {
    std::cout << "base REGION Size: " << baseRegion.GetSize() << " ";
    std::cout << "base REGION Index: " << baseRegion.GetIndex() << " ";

    std::cout << "base IMAGE Size: "
        << this->m_BaseImage->GetLargestPossibleRegion().GetSize() << " ";
    std::cout << "Search IMAGE Size: "
        << this->m_SearchImage->GetLargestPossibleRegion().GetSize() << std::endl;

    std::cout << "Nhood gap: "<<this->m_NhoodGap << ". No. Nhoods: "
        <<nNhoodsTotal<<" : "<<nNhoods;
    std::cout << " Nhood width: " << this->m_NhoodWidth
        << " Block gap: " << this->m_BlockGap
        << " Block width: " << this->m_BlockWidth << std::endl;
  }

  return baseRegion;
}


template < typename ImageType >
int BlockMatcher<ImageType>
::ComputeNhoodTotal(RegionType matchingRegion) const
{
  int nNhoodsTotal = 1;
  SizeType regionSize = matchingRegion.GetSize();
  for (int ii = 0; ii < ImageType::ImageDimension; ++ii) {
    nNhoodsTotal *= regionSize[ii]
              / this->m_NhoodGap + 1;
  }

  return nNhoodsTotal;
}

template < typename ImageType >
void BlockMatcher<ImageType>
::GetOffsets(
    PointListType &blockPositionsOut,
    PointListType &blockMatchPositionsOut
)
{
  this->BeforeGetOffsets();
  //Iterate over neighbourhoods in Base image (using nhood gap)
  RegionType baseRegion = this->m_BaseImage->GetLargestPossibleRegion();
  baseRegion = this->ComputeRegionForMatching(baseRegion);
  this->GetOffsets(baseRegion, blockPositionsOut, blockMatchPositionsOut);
}

template < typename ImageType >
void BlockMatcher<ImageType>
::GetOffsets( RegionType baseRegion, PointListType &blockPositionsOut,
      PointListType &blockMatchPositionsOut ) const
{
  //Make an iterator to go through all possible nhoods
  SparseConstIterator nhoodIterator( this->m_BaseImage, baseRegion );
  nhoodIterator.SetIteratorGap( this->m_NhoodGap );
  int nNhoodsTotal = this->ComputeNhoodTotal(baseRegion);

  //Ensure vectors are empty and reserve space for all the neighbourhoods
  IndexListType blockPositions, blockMatchPositions;

  blockPositions.clear();
  blockMatchPositions.clear();
  blockPositions.reserve( nNhoodsTotal );
  blockMatchPositions.reserve( nNhoodsTotal );

  //Create a vector to store the scores
#ifdef USE_FLOAT_SCORE_BUFFER
  std::vector<float> blockMatchScores, blockMatchScores2ndBest;
#else
  std::vector<double> blockMatchScores, blockMatchScores2ndBest;
#endif
  blockMatchScores.reserve( nNhoodsTotal );
  blockMatchScores2ndBest.reserve( nNhoodsTotal );

  // Loop over the positions within the base image.
  //These are called neighbourhoods, and the iterator points
  // to the neighbourhood centres
  if( this->m_Verbosity >= 3 )
    std::cout << "Initial position   Final position   Rho" << std::endl;

//------------------------------------------------------------------------
//Implement variance cut-off (not used)
#if 0
  //Get variances of all nhood base blocks (so boring blocks can be culled)
  std::vector<double> base_block_variances;
  base_block_variances.reserve( static_cast<int>(nNhoodsTotal * 1.1) );
  for( nhoodIterator.GoToBegin(); !nhoodIterator.IsAtEnd(); ++nhoodIterator )
  {
    //Create a block for looking at regions surrounding each neighbourhood centre
    // but in the base image (not the search image)
    IndexType baseIndex = nhoodIterator.GetIndex();

    SizeType blockSize;
    blockSize.Fill( this->m_BlockWidth );
    RegionType baseBlock;
    baseBlock.SetIndex( baseIndex );
    baseBlock.SetSize( blockSize );

    //We need an iterator to go through the voxels of each block
    //One iterator for the base image block
    ConstIterator baseImageIterator( this->m_BaseImage, baseBlock );

    double count = 0, sum_int = 0, sum_int_sqr = 0;

    for( baseImageIterator.GoToBegin(); !baseImageIterator.IsAtEnd();
        ++baseImageIterator )
    {
      count++;
      sum_int += baseImageIterator.Get();
      sum_int_sqr += baseImageIterator.Get() * baseImageIterator.Get();
    }

    double mean = sum_int / count;
    double variance = sum_int_sqr / count - mean*mean;

    base_block_variances.push_back( variance );
  }

  double variance_cutoff;
  {
    //Store the sorted variances
    std::vector<double> base_block_variances_sorted(base_block_variances);
    //What would the position of the nth element be
    unsigned int cutoff_index  =
        static_cast<unsigned int>( this->m_PortionMatchesKept * base_block_variances_sorted.size() );
    //Put the nth element in place
    std::nth_element( base_block_variances_sorted.begin(),
        base_block_variances_sorted.begin()+ cutoff_index,
        base_block_variances_sorted.end() );

    variance_cutoff = cutoff_index >= base_block_variances_sorted.size() ?
        base_block_variances_sorted.back() :
        base_block_variances_sorted[cutoff_index];
  }
#endif
//------------------------------------------------------------------------

  unsigned int nhood_index = 0;
  for( nhoodIterator.GoToBegin(); !nhoodIterator.IsAtEnd(); ++nhoodIterator, ++nhood_index )
  {
    //If the variance is too low, skip this block
    //if( base_block_variances[nhood_index] <= variance_cutoff )
    //  continue;

    //Get starting index
    IndexType nhoodStart, nhoodEnd;
    IndexType searchStart;
    searchStart.Fill(0); //to be like a2k
    IndexType searchEnd;

    for( int ii=0; ii<ImageType::ImageDimension; ++ii )
      searchEnd[ii] =
          this->m_BaseImage->GetLargestPossibleRegion().GetSize()[ii] -
          this->m_BlockWidth + 1;

    SizeType currentNhoodSize;

    //Get the start and end of the base iterator to ensure that
    //the nhoods are centred on the central nhood -
    //ITK makes the code look ugly. Sorry
    for( int ii=0; ii<ImageType::ImageDimension; ++ii )
    {
      nhoodStart[ii] = nhoodIterator.GetIndex()[ii] - ((this->m_NhoodWidth)/2) * this->m_BlockGap;
      nhoodEnd[ii] = nhoodStart[ii] + this->m_NhoodWidth * this->m_BlockGap;

      //Make sure we start within image boundaries
      while( nhoodStart[ii] < searchStart[ii] )
        nhoodStart[ii] += this->m_BlockGap;
      while( nhoodEnd[ii] > searchEnd[ii] ) //ND >= Can include end
        nhoodEnd[ii] -= this->m_BlockGap;

      currentNhoodSize[ii] = nhoodEnd[ii] - nhoodStart[ii];
    }

    //Create a block within the search image
    RegionType searchRegion;
    searchRegion.SetSize( currentNhoodSize );
    searchRegion.SetIndex( nhoodStart );

    //Create a sparse iterator within the search images for the different
    //locations the block will be within
    SparseConstIterator blockIterator( this->m_SearchImage, searchRegion );
    blockIterator.SetIteratorGap( this->m_BlockGap );
    blockIterator.GoToBegin();

    //Store the first block positions
    blockPositions.push_back( nhoodIterator.GetIndex() );
    blockMatchPositions.push_back( blockIterator.GetIndex() );
    blockMatchScores.push_back(-3);

    blockMatchScores2ndBest.push_back(-3);

    //Tell the user what going on
    if( this->m_Verbosity >= 3 )
    {
      std::cout
      << nhoodIterator.GetIndex()[0] << " "
      << nhoodIterator.GetIndex()[1] << " "
      << nhoodIterator.GetIndex()[2] << " ";
    }

    //Loop over positions WITHIN the neighbourhood
    //The base image block will move around the centre of the neighbourhood,
    //Hence this inner loop
    for( blockIterator.GoToBegin(); !blockIterator.IsAtEnd(); ++blockIterator )
    {
      //Create a block for looking at regions surrounding each neighbourhood centre
      // but in the base image (not the search image)
      IndexType baseIndex = nhoodIterator.GetIndex();
      SizeType blockSize;
      blockSize.Fill( this->m_BlockWidth );
      RegionType baseBlock;
      baseBlock.SetIndex( baseIndex );
      baseBlock.SetSize( blockSize );

      //We need an iterator to go through the voxels of each block
      //One iterator for the base image block
      ConstIterator baseImageIterator( this->m_BaseImage, baseBlock );
      MaskConstIterator baseMaskIterator( this->GetBaseMask(), baseBlock );

      //And one iterator for the search image block.
      //The search image block's position is set by the iterator over
      //neighbourhoods in the search image
      IndexType searchIndex = blockIterator.GetIndex();
      RegionType searchBlock;
      searchBlock.SetIndex( searchIndex );
      searchBlock.SetSize( blockSize );

      ConstIterator searchImageIterator( this->m_SearchImage, searchBlock );
      MaskConstIterator searchMaskIterator( this->GetSearchMask(), searchBlock );

      if( !metricFunc ) {
        std::ostringstream ss;
        ss << "ERROR: blockmetric function null pointer exception. Exiting!" << std::endl; \
        throw std::runtime_error(ss.str());
      }
      double metricValue = (this->*metricFunc)( baseImageIterator, searchImageIterator,
          baseMaskIterator, searchMaskIterator );

      //Store the best match so far
      if( blockMatchScores.back() < metricValue )
      {
        blockMatchScores2ndBest.back() = blockMatchScores.back();
        blockMatchScores.back() = metricValue;
        blockMatchPositions.back() = blockIterator.GetIndex();
      }
      else if ( blockMatchScores2ndBest.back() < metricValue )
      {
        blockMatchScores2ndBest.back() = metricValue;
      }
    }

    if( this->m_Verbosity >= 3 )
      std::cout
      << blockMatchPositions.back()[0] << " "
      << blockMatchPositions.back()[1] << " "
      << blockMatchPositions.back()[2] << " "
      << blockMatchScores.back()
      << std::endl;
  }

  this->PostProcessOffsets(blockPositions, blockMatchPositions,
                            blockMatchScores, blockMatchScores2ndBest,
                            blockPositionsOut, blockMatchPositionsOut );
}

#ifdef USE_FLOAT_SCORE_BUFFER
template < typename ImageType >
void BlockMatcher<ImageType>
::PostProcessOffsets(
    const IndexListType &blockPositions, const IndexListType &blockMatchPositions,
    const std::vector<float> &blockMatchScores, const std::vector<float> &blockMatchScores2ndBest,
    PointListType &blockPositionsOut, PointListType &blockMatchPositionsOut ) const
#else
template < typename ImageType >
void BlockMatcher<ImageType>
::PostProcessOffsets(
    const IndexListType &blockPositions, const IndexListType &blockMatchPositions,
    const std::vector<double> &blockMatchScores, const std::vector<double> &blockMatchScores2ndBest,
    PointListType &blockPositionsOut, PointListType &blockMatchPositionsOut ) const
#endif
{
  //Copy all elements up to cutoff_similarity out - because we're selecting particular elements
  //We have to do this manually
  blockMatchPositionsOut.clear();
  blockMatchPositionsOut.reserve( blockMatchScores.size() );
  blockPositionsOut.clear();
  blockPositionsOut.reserve( blockMatchScores.size() );

  int nLowScores = 0, nMultimatch = 0;
  double score_avg=0.f, score_min = 1e200, score_max = -1;

#ifdef DEBUG_GPU_CPU_CONSISTENCY
  std::cout << "USING itkBlockMatcher.txx" << std::endl;
  this->scoresOut.clear();
  this->scoresOut.reserve(blockMatchScores.size());
  this->scoresOut2nd.clear();
  this->scoresOut2nd.reserve(blockMatchScores.size());
#endif

  for( unsigned int currentIndex = 0; currentIndex<blockMatchScores.size(); ++currentIndex )
  {
    double currentBestScore = blockMatchScores[currentIndex];
    double current2ndBestScore = blockMatchScores2ndBest[currentIndex];

    //If we don't have multiple good matches, keep the block
    if( currentBestScore > 1e-7 )
    //if (true)
    {

      if( fabs( currentBestScore - current2ndBestScore ) > 1e-4 )
      //if (true)
      {
        //Note these are absolute positions
        blockPositionsOut.push_back( this->blockPositionToPoint( blockPositions[currentIndex] ) );

        //Note these are absolute positions
        blockMatchPositionsOut.push_back( this->blockPositionToPoint( blockMatchPositions[currentIndex] ) );
        score_avg += currentBestScore;
        score_min = std::min(currentBestScore,score_min);
        score_max = std::max(currentBestScore,score_max);

#ifdef DEBUG_GPU_CPU_CONSISTENCY
        this->scoresOut.push_back(currentBestScore); //DEBUG ONLY, SINGLE THREAD, RACE CONDITION IN MULTITHREADING
        this->scoresOut2nd.push_back(current2ndBestScore);
#endif

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
    std::ostringstream ss;
    ss<< "WARNING: No block matches occurred. Images probably have"
        <<" zero overlap. Try supplying a better initial transform."
        <<std::endl;
    throw std::runtime_error(ss.str());
  }
}


} // end namespace ITK

#endif
