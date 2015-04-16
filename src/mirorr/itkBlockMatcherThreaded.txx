/*=========================================================================
Program: mirorr
Module: itkBlockMatcherThreaded.txx
Author: David Rivest-Henault
Created: 02 Nov 2012

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include "itkBlockMatcherThreaded.h"


namespace itk
{
template < typename ImageType >
void
BlockMatcherThreaded<ImageType>::GetOffsets( SuperPointListType &blockPositions,
    SuperPointListType &blockMatchPositions )
{
  if (this->m_Verbosity >=2) {
    std::cout
        << "Executing BlockMatcherThreaded::GetOffsets(...) with "
        << this->GetNumberOfThreads() << " thread(s)" << std::endl;
  }
  // Set up the multithreaded processing
  ThreadStruct str;
  str.Matcher = this;

  this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(this->ThreaderCallback, &str);

  this->BeforeGetOffsets();
  this->InitializeLists();

  // multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();

  //Reduction
  this->ReduceLists(blockPositions, blockMatchPositions);
}

template<class ImageType>
void
BlockMatcherThreaded<ImageType>
::GetOffsets( SuperRegionType baseRegion, SuperPointListType &blockPositionsOut,
    SuperPointListType &blockMatchPositionsOut ) const
{
  Superclass::GetOffsets(baseRegion, blockPositionsOut, blockMatchPositionsOut);
}


template<class ImageType>
void
BlockMatcherThreaded<ImageType>::InitializeLists()
{
  this->blockMatchPositionsFractions.resize(this->GetMultiThreader()->GetNumberOfThreads());
  for (size_t i = 0; i<this->blockMatchPositionsFractions.size();i++) {
    this->blockMatchPositionsFractions[i].resize(0);
  }
  this->blockPositionsFractions.resize(this->GetMultiThreader()->GetNumberOfThreads());
  for (size_t i = 0; i<this->blockPositionsFractions.size();i++) {
    this->blockPositionsFractions[i].resize(0);
  }
}

template<class ImageType>
void
BlockMatcherThreaded<ImageType>
::ReduceLists(PointListType &blockPositions, PointListType &blockMatchPositions)
{
  blockPositions.clear();
  for (size_t i=0; i<this->blockPositionsFractions.size(); i++) {
    PointListType &ref = this->blockPositionsFractions[i];
    blockPositions.insert(blockPositions.end(), ref.begin(), ref.end());
  }

  blockMatchPositions.clear();
  for (size_t i=0; i<this->blockMatchPositionsFractions.size(); i++) {
    PointListType &ref = this->blockMatchPositionsFractions[i];
    blockMatchPositions.insert(blockMatchPositions.end(), ref.begin(), ref.end());
  }
}


// DRH: Split the image region into a maximum of 'threadCount' sub regions. Region
// splitting is conducted along the greatest dimensions  of the
// 'outputPtr->GetLargestPossibleRegion()'.
// The region associated with 'threadId' is stored in 'splitRegion', and the method
// return the maximum number of thread that can be used.
template< class ImageType >
unsigned int
BlockMatcherThreaded< ImageType >
::SplitRequestedRegion(unsigned int threadId, unsigned int threadCount, RegionType & splitRegion)
{
  // Get the output pointer
  typename ImageType::ConstPointer outputPtr = this->m_BaseImage;

  // Initialize the splitRegion to the output requested region
  splitRegion = this->ComputeRegionForMatching(outputPtr->GetLargestPossibleRegion());
  //splitRegion = this->ComputeRegionForMatching(outputPtr->GetRequestedRegion()); //DRH: Better?
  typename ImageType::IndexType  splitIndex = splitRegion.GetIndex();
  typename ImageType::SizeType splitSize = splitRegion.GetSize();

  const typename ImageType::SizeType requestedRegionSize = splitSize;

  // split on the largest dimension available
  int splitAxis = (splitSize[0] > splitSize[1]) \
      ? ( (splitSize[0] > splitSize[2]) ? 0 : 2 ) \
          : ( (splitSize[1] > splitSize[2]) ? 1 : 2 );

  if (splitSize[splitAxis] - this->m_BlockWidth < this->m_NhoodGap) {
    itkDebugMacro("  Cannot Split");
    return 1;
  }

  // compute the number of neighbourhoods along that axis
  unsigned int nNhoods = (requestedRegionSize[splitAxis] - 1) / this->m_NhoodGap + 1;

  // determine the actual number of neighbourhoods that will be processed per thread
  unsigned int nNhoodsPerThread = nNhoods / threadCount;
  unsigned int nNhoodsResidual = nNhoods % threadCount;

  unsigned int maxThreadIdUsed = (nNhoodsPerThread > 0) ? threadCount : (nNhoodsResidual-1);

  // determine the actual number of neighbourhoods that will be processed for THIS thread
  unsigned int nNhoodsThisThread = nNhoodsPerThread;
  if (threadId < nNhoodsResidual) {
    nNhoodsThisThread++;
  }

  // determine the index of the first neighbourhoods that will be processed for THIS thread
  // Note: If there is a nNhoodsResidual, the first 'nNhoodsResidual' threads process an
  //       extra neighbourhood.
  unsigned int nNhoodIdx = threadId * nNhoodsPerThread;
  nNhoodIdx += (threadId < nNhoodsResidual) ? threadId : nNhoodsResidual;

  // Split the region
  splitIndex[splitAxis] += this->m_NhoodGap * nNhoodIdx;
  splitSize[splitAxis] = 1 + this->m_NhoodGap * (nNhoodsThisThread-1);

  // set the split region
  splitRegion.SetIndex(splitIndex);
  splitRegion.SetSize(splitSize);

  itkDebugMacro("  Split Piece: " << splitRegion);

  return maxThreadIdUsed + 1;
}


// Callback routine used by the threading library. This routine just calls
// the ThreadedGenerateData method after setting the correct region for this
// thread.
template<class ImageType>
  ITK_THREAD_RETURN_TYPE BlockMatcherThreaded<ImageType>
  ::ThreaderCallback(void *arg)
{
  ThreadStruct *str;
  ThreadIdType total, threadId, threadCount;

  threadId = ( (MultiThreader::ThreadInfoStruct *)( arg ) )->ThreadID;
  threadCount = ( (MultiThreader::ThreadInfoStruct *)( arg ) )->NumberOfThreads;

  str = (ThreadStruct *)( ( (MultiThreader::ThreadInfoStruct *)( arg ) )->UserData );

  // execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  typename ImageType::RegionType splitRegion;
  total = str->Matcher->SplitRequestedRegion(threadId, threadCount,
      splitRegion);

  if ( threadId < total )
  {
    if (threadId > static_cast<ThreadIdType>(str->Matcher->blockPositionsFractions.size()) ||
        threadId > static_cast<ThreadIdType>(str->Matcher->blockMatchPositionsFractions.size()) )
    {
      std::stringstream ss;
      ss
       << "ERROR: threadId bigger than expected. "
       << " threadId: " << threadId
       << " number of slots available: " << str->Matcher->blockPositionsFractions.size() << std::endl;
      throw std::runtime_error(ss.str());
    }

    PointListType &blockPositionsOut = str->Matcher->blockPositionsFractions[threadId];
    PointListType &blockMatchPositionsOut = str->Matcher->blockMatchPositionsFractions[threadId];

    str->Matcher->GetOffsets(splitRegion, blockPositionsOut, blockMatchPositionsOut);
  }
  // else
  //   {
  //   otherwise don't use this thread. Sometimes the threads dont
  //   break up very well and it is just as efficient to leave a
  //   few threads idle.
  //   }

  return ITK_THREAD_RETURN_VALUE;
}

} //namespace itk

//
