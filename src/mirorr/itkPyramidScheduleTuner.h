/*=========================================================================
Program: mirorr
Module: milxScheduleTuner.h
Author: Nicholas Dowson
Created: 19 Oct 2010

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#ifndef __ITKSCHEDULETUNER_H_
#define __ITKSCHEDULETUNER_H_

#include <itkRecursiveMultiResolutionPyramidImageFilter.h>
#include <iostream>
#include <iomanip>
#include <itkObject.h>

namespace itk
{

/** Based on the sizes of two images, create a pyramid schedule
 *  that optimally matches the schedules based on the spacings
 *  within the two images and their sizes
 */
template <unsigned int DIMENSION = 3 >
class ITK_EXPORT PyramidScheduleTuner :
public Object
{
public:
  /** Standard class typedefs. */
  typedef PyramidScheduleTuner Self;
  typedef LightObject SuperClass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PyramidScheduleTuner, LightObject);

  //!Useful type definitions
  //static const unsigned int DIMENSION = 3;
  static const unsigned int Dimension = DIMENSION;
  typedef itk::Size<DIMENSION> SizeType;
  typedef itk::Vector<double, DIMENSION> SpacingType;
  typedef itk::Array2D<unsigned int> ScheduleType;

  //! Create the pyramid schedules
  void Update();

  //! Set and get images spacings and sizes
  itkSetMacro( FixedSpacing, SpacingType )
  itkGetConstReferenceMacro( FixedSpacing, SpacingType )
  itkSetMacro( MovingSpacing, SpacingType )
  itkGetConstReferenceMacro( MovingSpacing, SpacingType )

  itkSetMacro( FixedSize, SizeType )
  itkGetConstReferenceMacro( FixedSize, SizeType )
  itkSetMacro( MovingSize, SizeType )
  itkGetConstReferenceMacro( MovingSize, SizeType )

  void SetBothSpacings( const SpacingType & in )
  {
    SetFixedSpacing( in );
    SetMovingSpacing( in );
  }
  void SetBothSizes( const SizeType & in )
  {
    SetFixedSize( in );
    SetMovingSize( in );
  }

  //! Crop a level
  int GetCroppedLevel( int level, int default_level = 1 );

  //! Set / Get the minimum allowed length along any dimension
  //! when subsampling an image
  itkSetMacro( MinLength, unsigned int )
  itkGetMacro( MinLength, unsigned int )

  //! Get the full schedules that have been generated
  ScheduleType GetFullFixedSchedule() const { return m_FixedSchedule; };
  ScheduleType GetFullMovingSchedule() const { return m_MovingSchedule; };

  //! Get the schedules after cropping with min and max level
  ScheduleType GetFixedSchedule() { return CropSchedule( m_FixedSchedule ); };
  ScheduleType GetMovingSchedule() { return CropSchedule( m_MovingSchedule ); };

  //! Set / Get the minimum and maximum levels of the pyramid to use
  //! 0 means use default, negative means count backwards from the top
  //! i.e. 1 means top level with most downsampling and
  //! -1 means bottom level with least downsampling
  itkSetMacro( LevelMin, int )
  itkGetMacro( LevelMin, int )
  itkSetMacro( LevelMax, int )
  itkGetMacro( LevelMax, int )
  itkSetMacro( MaxLevelNumber, unsigned int )
  itkGetMacro( MaxLevelNumber, unsigned int )

  void PrintPyramidInfo(std::ostream& os, int verbosity=1) const;

protected:
  PyramidScheduleTuner() : m_MinLength(32), m_MaxLevelNumber(1024), m_LevelMin(0), m_LevelMax(0) {};
  virtual ~PyramidScheduleTuner() {};
  void PrintSelf(std::ostream&os, Indent indent) const;

  void UpdateVintage();   // Implement the original optimised, but asymmetric, pyramid schedule
  void UpdateSymmetric(); // Implement an optimised symmetric pyramid schedule
  void UpdateSimple();    // Implement the simplest pyramid schedule: just dived the image size by 2 at each level

  //! Sizes of the fixed and moving images
  SizeType m_FixedSize;
  SizeType m_MovingSize;
  //! Spacings within the fixed and moving images
  SpacingType m_FixedSpacing;
  SpacingType m_MovingSpacing;

  //!Full schedules for the fixed and moving images
  ScheduleType m_FixedSchedule;
  ScheduleType m_MovingSchedule;

  //!The minimum number of voxels along any one dimension of the image
  unsigned int m_MinLength;
  unsigned int m_MaxLevelNumber;

  //!The minimum and maximum pyramid levels to use
  int m_LevelMin;
  int m_LevelMax;

  /** The pyramid levels after being corrected into positive non-zero values */
//  unsigned int m_OutputMinLevel, m_OutputMaxLevel;

private:
  PyramidScheduleTuner(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  void GetCroppedMinMax( unsigned int &min, unsigned int &max ) const;

  //! Crop the full schedule based on the minimum and maximum
  ScheduleType CropSchedule( const ScheduleType& in_schedule ) const
  {
    unsigned int minLevel = 0, maxLevel = 0;
    GetCroppedMinMax( minLevel, maxLevel );

    //Copy at schedule based on min and max
    ScheduleType schedule =  ScheduleType( maxLevel - minLevel + 1, DIMENSION );
    unsigned int outlevel = 0;
    for( unsigned int level = minLevel; level<=maxLevel; ++level, ++outlevel )
      for( unsigned int dim = 0; dim<DIMENSION; ++dim )
        schedule[outlevel][dim] = in_schedule[level-1][dim];

    return schedule;
  }
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPyramidScheduleTuner.txx"
#endif

#endif /* MILXSCHEDULETUNER_H_ */
