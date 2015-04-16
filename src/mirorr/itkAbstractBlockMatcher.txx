/*=========================================================================
Program: mirorr
Module: itkAbstractBlockMatcher.txx
Author: David Rivest-Henault
Created: 02 Nov 2012

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
#ifndef __itkAbstractBlockMatcher_txx
#define __itkAbstractBlockMatcher_txx

#include "itkAbstractBlockMatcher.h"

namespace itk
{

  template < typename ImageType >
  AbstractBlockMatcher<ImageType>::AbstractBlockMatcher():
  m_BlockGap( 1 ),
  m_BlockWidth( 4 ),
  m_NhoodGap( 2 ),
  m_NhoodWidth( 5 ),
  m_Verbosity( 0 ),
  m_PortionMatchesKept( 0.5 ),
  m_Padding( 0 ),
#ifdef USE_NPW
    m_NPWbins( 32 ),
#ifndef USE_OPENCL
    m_NPWscaleFactor( 2 ),
    m_NPWshapeExpansion( false ),
#endif
#endif
  m_BaseImage(0),
  m_SearchImage(0)
  {
  }

  template < typename ImageType >
  void AbstractBlockMatcher<ImageType>::SetBaseImage( const ImageType * movingImage )
  {
    if( this->m_Verbosity >= 2 )
      std::cout << "setting Moving Image to " << movingImage << ". " << movingImage->GetLargestPossibleRegion().GetSize() << std::endl;

    if (this->m_BaseImage.GetPointer() != movingImage )
      {
        this->m_BaseImage = movingImage;
        // Process object is not const-correct so the const_cast is required here
        this->ProcessObject::SetNthInput(0, const_cast< ImageType *>( movingImage ) );
        this->Modified();
      }
  }

  template < typename ImageType >
  void AbstractBlockMatcher<ImageType>::SetSearchImage( const ImageType * fixedImage )
  {
    if( this->m_Verbosity >= 2 )
      std::cout << "setting Fixed Image to " << fixedImage << ". " << fixedImage->GetLargestPossibleRegion().GetSize() << std::endl;

    if (this->m_SearchImage.GetPointer() != fixedImage )
      {
        this->m_SearchImage = fixedImage;
        // Process object is not const-correct so the const_cast is required here
        this->ProcessObject::SetNthInput(1, const_cast< ImageType *>( fixedImage ) );
        this->Modified();
      }
  }

  template < typename ImageType >
  void AbstractBlockMatcher<ImageType>::SetPadding( unsigned int in )
  {
    m_Padding = (in<50)?in:50;
    m_Padding = (m_NhoodGap*m_NhoodWidth/2+m_BlockWidth/2+2<m_Padding)?m_Padding:m_NhoodGap*m_NhoodWidth/2+m_BlockWidth/2+2;
  }


  template < typename ImageType >
  void AbstractBlockMatcher<ImageType>::SetPortionMatchesKept( double in )
  {
    if( this->m_Verbosity >= 3 )
      std::cout << "..:: SetPortionMatchesKept ::.." << std::endl;
    if( in == 0.0 )
      in = 0.5;

    if( in < 0.0 || in > 1.0 )
      {
        std::stringstream ss;
        ss << "itkBlockMatcher::SetPortionMatchesKept input is "
            << in << ". It should be between 0 and 1.0.";
        throw std::out_of_range( ss.str() );
      }
    m_PortionMatchesKept = in;
  }

  template < typename ImageType >

  void AbstractBlockMatcher<ImageType>::SwapImageOrder()
  {
    typename ImageType::ConstPointer tmpImg = this->m_BaseImage;
    this->m_BaseImage = this->m_SearchImage;
    this->m_SearchImage = tmpImg;

    typename MaskType::Pointer tmpMask = this->m_BaseMask;
    this->m_BaseMask = this->m_SearchMask;
    this->m_SearchMask = tmpMask;
  }

  template < typename ImageType >

  void AbstractBlockMatcher<ImageType>::displayInfo()
  {
    if( this->m_Verbosity >= 2 )
      std::cout << "..:: displayInfo ::.." << std::endl;

    std::cout << "Block Gap: " << m_BlockGap << "px  "
        << "Width: " << m_BlockWidth << "px  "
        << "N'hood Gap: " << m_NhoodGap << "px  "
        << "Width: " << m_NhoodWidth << "blk. "
        << "Portion Kept: " <<  m_PortionMatchesKept << " "
        << std::endl;
  }

  template < typename ImageType >
  void AbstractBlockMatcher<ImageType>::PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf(os,indent);

    os << indent << "Block Gap: " << m_BlockGap << "px  "
        << "Width: " << m_BlockWidth << "px  "
        << "N'hood Gap: " << m_NhoodGap << "px  "
        << "Width: " << m_NhoodWidth << "blk. "
        << "Portion Kept: " <<  m_PortionMatchesKept << " "
        << std::endl;
  }

}

#endif
