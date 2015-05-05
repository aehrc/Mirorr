/*=========================================================================
Program: MILX MixView
Module: itkCimgUtil.h
Author: Daan Broekhuizen
Modified by:
Language: C++
Created: Tue 25 May 2010 12:23:13 EST

Copyright: (c) 2009 CSIRO, Australia.

This software is protected by international copyright laws.
Any unauthorised copying, distribution or reverse engineering is prohibited.

Licence:
All rights in this Software are reserved to CSIRO. You are only permitted
to have this Software in your possession and to make use of it if you have
agreed to a Software License with CSIRO.

BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
=========================================================================*/
#ifndef __ITKCIMGUTIL_H__
#define __ITKCIMGUTIL_H__

#include <boost/shared_ptr.hpp>

//Useful routines for converting between itk and cimg

namespace itk
{
  namespace cimg = cimg_library;

  template <typename TVoxel>
    boost::shared_ptr< cimg::CImg<TVoxel> >
    get_itk2cimg( typename itk::OrientedImage<TVoxel,3>::ConstPointer input )
    {
      if( input.IsNull() )
        throw std::invalid_argument("itk2cimg Error: Input itkImagePointer is NULL.\n");

      typename itk::OrientedImage<TVoxel,3>::RegionType region = input->GetLargestPossibleRegion();
      typename itk::OrientedImage<TVoxel,3>::SizeType size = region.GetSize();
      const TVoxel * buffer = input->GetBufferPointer();

      boost::shared_ptr< cimg::CImg<TVoxel> >
        out( new typename cimg::CImg<TVoxel>::CImg( buffer, size[0], size[1], size[2], 1, false ) );

      return out;
    }

  template <typename TVoxel>
    boost::shared_ptr< cimg::CImg<TVoxel> >
    get_itk2cimg( typename itk::OrientedImage<TVoxel,3>::Pointer input )
    {
      if( input.IsNull() )
        throw std::invalid_argument("itk2cimg Error: Input itkImagePointer is NULL.\n");

      typename itk::OrientedImage<TVoxel,3>::RegionType region = input->GetLargestPossibleRegion();
      typename itk::OrientedImage<TVoxel,3>::SizeType size = region.GetSize();
      const TVoxel * buffer = input->GetBufferPointer();

      boost::shared_ptr< cimg::CImg<TVoxel> >
        out( new typename cimg::CImg<TVoxel>::CImg( buffer, size[0], size[1], size[2], 1, false ) );

      return out;
    }

  template <typename TVoxel>
    boost::shared_ptr< cimg::CImg<TVoxel> >
    itk2cimg( typename itk::OrientedImage<TVoxel,3>::Pointer input )
    {
      if( input.IsNull() )
        throw std::invalid_argument("itk2cimg Error: Input itkImagePointer is NULL.\n");

      typename itk::OrientedImage<TVoxel,3>::RegionType region = input->GetLargestPossibleRegion();
      typename itk::OrientedImage<TVoxel,3>::SizeType size = region.GetSize();
      const TVoxel * buffer = input->GetBufferPointer();

      boost::shared_ptr< cimg::CImg<TVoxel> >
        out( new typename cimg::CImg<TVoxel>::CImg( buffer, size[0], size[1], size[2], 1, true ) );

      return out;
    }

  template <typename TVoxel>
  boost::shared_ptr< const cimg::CImg<TVoxel> >
  itk2cimg( const typename itk::OrientedImage<TVoxel,3u>::ConstPointer &input )
  {
    if( input.IsNull() )
      throw std::invalid_argument("itk2cimg Error: Input itkImagePointer is NULL.\n");
    const typename itk::OrientedImage<TVoxel,3u>::RegionType region = input->GetLargestPossibleRegion();
    const typename itk::OrientedImage<TVoxel,3u>::SizeType size = region.GetSize();
    const TVoxel * buffer = input->GetBufferPointer();

    boost::shared_ptr< const cimg::CImg<TVoxel> >
    out( new const typename cimg::CImg<TVoxel>::CImg( buffer, size[0], size[1], size[2], 1, true ) );

    return out;
  }


  template <typename TVoxel>
  typename itk::OrientedImage<TVoxel,3>::Pointer
  cimg2itk( boost::shared_ptr< const cimg::CImg<TVoxel> > input )
  {
      //Create a new ITK image
      typename itk::OrientedImage<TVoxel,3>::Pointer out = itk::OrientedImage<TVoxel,3>::New();

      typedef typename itk::OrientedImage<TVoxel,3>::RegionType RegionType;
      typedef typename itk::OrientedImage<TVoxel,3>::SizeType SizeType;
      SizeType size;
      size[0] = input->dimx();
      size[1] = input->dimy();
      size[2] = input->dimz();

      RegionType region( size );

      out->SetRegions( region );

      //Allocate the memory
      out->Allocate();

      //Copy the voxels
      itk::ImageRegionIterator< itk::OrientedImage<TVoxel,3> > outIter( out, region );
      const TVoxel* inIter = input->begin();

      for( outIter.GoToBegin(); !outIter.IsAtEnd(); ++outIter )
    outIter.Set( *inIter++ );

      return out;
  }
}

#endif //__ITKCIMGUTIL_H__
