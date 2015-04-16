/*=========================================================================
Program: mirorr
Module: itkSparseImageRegionConstIterator.h
Author: Nicholas Dowson
Created: Mon 11 Feb 2009

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
#ifndef __itkSparseImageRegionConstIterator_h
#define __itkSparseImageRegionConstIterator_h

#include <itkImageRegionConstIterator.h>

#if defined(_MSC_VER) && !defined(ITK_LEAN_AND_MEAN)
#define ITK_LEAN_AND_MEAN
#endif

namespace itk
{

/*! \class SparseImageRegionConstIterator
 * \brief An iterator for sampling images with a gap (in pixels) between
 * each sample. The gap taken is the same in all dimensions. In other words,
 * in a 40x40x40 image where a normal iterator would touch 40^3=64000 voxels
 * a sparse iterator with a gap of 2 would touch only 20^3=8000 voxels.
 *
 * A default gap of 1 is always used. The points touched by the iterator
 * are aligned, i.e. the points form an undeformed rectangular grid in N-D
 * space.
 * The image iterator will always touch the first voxel in the image buffer.
 * This is somewhat restrictive and may be fixed in the future.
 * The iterator is templated on the image type, as the pixel type and
 * image dimension need to to be known beforehand.
 *
 * \par Note on the image size and the gap
 * The gap does NOT need need to be a divisor of the image size.
 * In a 5x5 image with voxels A1 B1 C1 D1 E1, A2..., A3..., ..., A5 ...
 * and a gap of 2, the iterator with touch  voxels A1 C1 E1, A3 C3 E3, A5 C5 E5.
 * If the gap is 3, the iterator will touch voxels A1 D1, A4 D4.
 *
 * \par Inputs
 * Two inputs are required:
 * a region within an image, and
 * a gap (in voxels) to take between successive positions that
 * are touched by the iterator.
 *
 * \par To do
 * 1. Allow different spacings for different dimensions.
 * 2. Allow different initial offsets.
 *
 * \ingroup ITKModulesCommon
 */

template<typename TImage>
class ITK_EXPORT SparseImageRegionConstIterator : public ImageRegionConstIterator<TImage>
{
public:
  /** Standard class type definitions. */
  typedef SparseImageRegionConstIterator    Self;
  typedef ImageRegionConstIterator<TImage>  Superclass;

  /** Dimension of the image the iterator walks.  This constant is needed so
   * functions that are templated over image iterator type (as opposed to
   * being templated over pixel type and dimension) can have compile time
   * access to the dimension of the image that the iterator walks. */
  itkStaticConstMacro(ImageIteratorDimension, unsigned int,
      Superclass::ImageIteratorDimension);

  /** Index typedef support. While this was already typdef'ed in the superclass
   * it needs to be redone here for this subclass to compile properly with gcc. */
  typedef typename Superclass::IndexType    IndexType;

  /** Size typedef support. While this was already typdef'ed in the superclass
   * it needs to be redone here for this subclass to compile properly with gcc. */
  typedef typename Superclass::SizeType SizeType;

  /** Region typedef support. */
  typedef typename Superclass::RegionType   RegionType;

  /** Image typedef support. While this was already typdef'ed in the superclass
   * it needs to be redone here for this subclass to compile properly with gcc. */
  typedef typename Superclass::ImageType    ImageType;

  /** PixelContainer typedef support. Used to refer to the container for
   * the pixel data. While this was already typdef'ed in the superclass
   * it needs to be redone here for this subclass to compile properly with gcc. */
  typedef typename Superclass::PixelContainer         PixelContainer;
  typedef typename Superclass::PixelContainerPointer  PixelContainerPointer;

  /** Internal Pixel Type */
  typedef typename Superclass::InternalPixelType   InternalPixelType;

  /** External Pixel Type */
  typedef typename Superclass::PixelType    PixelType;

  /**  Accessor type that convert data between internal and external
   *  representations. */
  typedef typename Superclass::AccessorType AccessorType;

#if ITK_VERSION_MAJOR < 4
  typedef typename Superclass::IndexValueType IndexValueType;
  typedef typename Superclass::IndexValueType OffsetValueType;
#else
  typedef typename ::itk::IndexValueType IndexValueType;
  typedef typename ::itk::IndexValueType OffsetValueType;
#endif
  

  
  /** Run-time type information (and related methods). */
  itkTypeMacro(SparseImageRegionConstIterator, ImageIterator);

  /** Default constructor. This is needed because we provide a cast
   * constructor. By default the iterator will behave like a normal
   * iterator with a step size of one. */
  SparseImageRegionConstIterator() : ImageRegionConstIterator<TImage>()
    {
    SetIteratorGap( 1 );
    }

  /** This constructor establishes an iterator to walk through a
   * particular region within an image. By default the iterator
   * will behave like a normal
   * iterator with a step size of one. */
  SparseImageRegionConstIterator( const ImageType *ptr,
      const RegionType &region)
  : ImageRegionConstIterator<TImage>(ptr, region)
    {
    SetIteratorGap( 1 ); //By default behave like a normal iterator
    }

  /** Constructor that can be used to cast from an ImageIterator to an
   * ImageRegionConstIterator. Many routines return an ImageIterator but for a
   * particular task, you may want an ImageRegionConstIterator.  Rather than
   * provide overloaded APIs that return different types of Iterators, itk
   * returns ImageIterators and uses constructors to cast from an
   * ImageIterator to a ImageRegionConstIterator. */
  SparseImageRegionConstIterator( const ImageIterator<TImage> &it)
    {
    this->ImageConstIterator<TImage>::operator=(it);

    IndexType ind = this->GetIndex(); //Get current position
    //Make sure index lines up with iterator gap
    SetIndex( ind );
    //    for( unsigned int ii=0; ii<IndexType::GetIndexDimension(); ++ii )
    //      ind[ii] = floor( ind[ii] / m_IteratorGap ) * m_IteratorGap;
    //
    //    this->m_SpanEndOffset = this->m_Offset + static_cast<long>(this->m_Region.GetSize()[0])
    //      - (ind[0] - this->m_Region.GetIndex()[0]);
    //    this->m_SpanBeginOffset = this->m_SpanEndOffset
    //      - static_cast<long>(this->m_Region.GetSize()[0]);
    }

  /** Constructor that can be used to cast from an ImageConstIterator to an
   * ImageRegionConstIterator. Many routines return an ImageIterator but for a
   * particular task, you may want an ImageRegionConstIterator.  Rather than
   * provide overloaded APIs that return different types of Iterators, itk
   * returns ImageIterators and uses constructors to cast from an
   * ImageIterator to a ImageRegionConstIterator. */
  SparseImageRegionConstIterator( const ImageConstIterator<TImage> &it)
    {
    this->ImageConstIterator<TImage>::operator=(it);

    IndexType ind = this->GetIndex(); //Get current position coordinates
    //Make sure index lines up with iterator gap
    SetIndex( ind );

    }

  /** Set the index. No bounds checking is performed. This is overridden
   * from the parent because we have an extra ivar.
   * \sa GetIndex */
  void SetIndex( const IndexType & _ind )
    {
    IndexType ind;
    //Make sure index lines up with iterator gap
    for( unsigned int ii=0; ii<IndexType::GetIndexDimension(); ++ii )
      ind[ii] = ( _ind[ii] / m_IteratorGap ) * m_IteratorGap;
    Superclass::SetIndex(ind);
    }

  /** Increment (prefix) the fastest moving dimension of the iterator's index.
   * This operator will constrain the iterator within the region (i.e. the
   * iterator will automatically wrap from the end of the row of the region
   * to the beginning of the next row of the region) up until the iterator
   * tries to moves past the last pixel of the region.  Here, the iterator
   * will be set to be one pixel past the end of the region.
   * \sa operator++(int) */
  Self &
  operator++()
    {
    if( (this->m_Offset+=m_IteratorGap) >= this->m_SpanEndOffset)
      {
      this->Increment();
      }
    return *this;
    }

  /** Decrement (prefix) the fastest moving dimension of the iterator's index.
   * This operator will constrain the iterator within the region (i.e. the
   * iterator will automatically wrap from the beginning of the row of the region
   * to the end of the next row of the region) up until the iterator
   * tries to moves past the first pixel of the region.  Here, the iterator
   * will be set to be one pixel past the beginning of the region.
   * \sa operator--(int) */
  Self & operator--()
    {
    if( (this->m_Offset-=m_IteratorGap) < this->m_SpanBeginOffset)
      {
      this->Decrement();
      }
    return *this;
    }

  /** Set the gap (in voxels) between each pair of successive samples
   * as the iterator traverses over the image.
   *
   * @param[in] gap The required gap.
   *
   * \par To do:
   * - Ensure gap is short enough to fit within image
   */
  void SetIteratorGap( unsigned long gap )
    {
    m_IteratorGap = gap;

    //Change end offset
    IndexType ind(this->m_Region.GetIndex());
    SizeType size(this->m_Region.GetSize());
    for (unsigned int i=0; i < ImageIteratorDimension; ++i)
      {
      ind[i] += ((static_cast<int>(size[i]) - 1) / m_IteratorGap) * m_IteratorGap;
      }
    this->m_EndOffset = this->m_Image->ComputeOffset( ind );
    this->m_EndOffset += m_IteratorGap;
    }

  /** Return the gap (in voxels) between each pair of successive samples
   * as the iterator traverses over the region.
   */
  unsigned long GetIteratorGap() { return m_IteratorGap; }

protected:
  unsigned long m_IteratorGap; //Stores the gap (in voxels) between two samples

private:
  void Increment(); /// Advance in a direction other than the fastest moving dimension
  void Decrement(); /// Go back in a direction other than the fastest moving dimension

};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_SparseImageRegionConstIterator(_, EXPORT, x, y) namespace itk { \
    _(1(class EXPORT SparseImageRegionConstIterator< ITK_TEMPLATE_1 x >)) \
    namespace Templates { typedef SparseImageRegionConstIterator< ITK_TEMPLATE_1 x > SparseImageRegionConstIterator##y; } \
}


#ifndef ITK_MANUAL_INSTANTIATION
# include "itkSparseImageRegionConstIterator.txx"
#endif

#if defined(_MSC_VER)
#undef ITK_LEAN_AND_MEAN
#endif

#endif
