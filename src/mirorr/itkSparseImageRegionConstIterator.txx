/*=========================================================================
Program: mirorr
Module: itkSparseImageRegionConstIterator.txx
Author: Nicholas Dowson
Created: Mon 11 Feb 2009

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/** An iterator for sampling images at fixed intervals
Todo:
 * Allow different spacings for different dimensions
 * Allow different initial offsets
 */

#ifndef __itkSparseImageRegionConstIterator_txx
#define __itkSparseImageRegionConstIterator_txx

#include "itkSparseImageRegionConstIterator.h"
#include <iostream>

#ifdef __INTEL_COMPILER
/* pointless comparision of unsigned integer with zero */
#pragma warning(disable: 186)
#endif

namespace itk
{

//----------------------------------------------------------------------------
// Increment when the fastest moving direction has reached its bound.
// This method should *ONLY* be invoked from the operator++() method.
template<class TImage>
void
SparseImageRegionConstIterator<TImage>
::Increment()
{
  // We have reached the end of the span (row), need to wrap around.

  // First back up N pixels, because we are going to use a different
  // algorithm to compute the next pixel
  this->m_Offset -= m_IteratorGap;

  // Get the index of the last pixel on the span (row)
  typename ImageIterator<TImage>::IndexType
  ind = this->m_Image->ComputeIndex( static_cast<OffsetValueType>(this->m_Offset) );

  const typename ImageIterator<TImage>::IndexType&
  startIndex = this->m_Region.GetIndex();
  const typename ImageIterator<TImage>::SizeType&
  size = this->m_Region.GetSize();

  // Increment along a row, then wrap at the end of the region row.
  unsigned int dim;

  // Check to see if we are past the last pixel in the region
  // Note that ++ind[0] moves to the next pixel along the row.
  ind[0] += m_IteratorGap;
  bool done = (ind[0] >= startIndex[0] +
      static_cast<IndexValueType>(size[0]) - 1);
  for (unsigned int i=1; done && i < ImageIteratorDimension; i++)
    {
    done = ( static_cast<int>(ind[i]+m_IteratorGap) > startIndex[i] +
        static_cast<IndexValueType>(size[i]) - 1);
    }//!Check this works

  // if the iterator is outside the region (but not past region end) then
  // we need to wrap around the region
  dim = 0;
  if (!done)
    {
    while ( ( dim+1 < ImageIteratorDimension )
        && (ind[dim] > startIndex[dim] +  static_cast<IndexValueType>(size[dim]) - 1) )
      {
      ind[dim] = startIndex[dim];
      ind[++dim] += m_IteratorGap;
      }
    }
  this->m_Offset = this->m_Image->ComputeOffset( ind );
  this->m_SpanEndOffset = this->m_Offset + static_cast<long>(size[0]);
  this->m_SpanBeginOffset = this->m_Offset;
}


//----------------------------------------------------------------------------
// Decrement when the fastest moving direction has reached its bound.
// This method should *ONLY* be invoked from the operator--() method.
template<class TImage>
void
SparseImageRegionConstIterator<TImage>
::Decrement()
{
  // We have pasted the beginning of the span (row), need to wrap around.

  // First move forward one pixel, because we are going to use a different
  // algorithm to compute the next pixel
  this->m_Offset += m_IteratorGap;

  // Get the index of the first pixel on the span (row)
  typename ImageIterator<TImage>::IndexType
  ind = this->m_Image->ComputeIndex( static_cast<IndexValueType>(this->m_Offset) );

  const typename ImageIterator<TImage>::IndexType&
  startIndex = this->m_Region.GetIndex();
  const typename ImageIterator<TImage>::SizeType&
  size = this->m_Region.GetSize();

  // Deccrement along a row, then wrap at the beginning of the region row.
  bool done;
  unsigned int dim;

  // Check to see if we are past the first pixel in the region
  // Note that --ind[0] moves to the previous pixel along the row.
  done = ( (ind[0]-=m_IteratorGap) = startIndex[0] - 1);
  for (unsigned int i=1; done && i < ImageIteratorDimension; i++)
    {
    done = (ind[i] == startIndex[i]);
    }

  // if the iterator is outside the region (but not past region begin) then
  // we need to wrap around the region
  dim = 0;
  if (!done)
    {
    while ( (dim < ImageIteratorDimension - 1)
        && (ind[dim] < startIndex[dim]) )
      {
      ind[dim] = startIndex[dim] +
        ((static_cast<IndexValueType>(size[dim])-1) / m_IteratorGap) * m_IteratorGap;
      ind[++dim] -= m_IteratorGap;
      }
    }
  this->m_Offset = this->m_Image->ComputeOffset( ind );
  this->m_SpanEndOffset = this->m_Offset + 1;
  this->m_SpanBeginOffset = this->m_SpanEndOffset - (static_cast<long>(size[0]-1)/m_IteratorGap)*m_IteratorGap - 1;
}

} // end namespace itk

#endif
