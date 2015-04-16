/*=========================================================================
Program: mirorr
Module: milxScheduleTuner.txx
Author: Nicholas Dowson
Created: 24 Nov 2010

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#ifndef ITKPYRAMIDSCHEDULETUNER_TXX_
#define ITKPYRAMIDSCHEDULETUNER_TXX_

#include "itkPyramidScheduleTuner.h"
#include "itkIOUtils.h"

template <unsigned int DIMENSION >
int
itk::PyramidScheduleTuner<DIMENSION>::
GetCroppedLevel( int level, int default_level )
{
  int totalLevels = m_FixedSchedule.rows();

  //Make Local copies of level min and max
  int f_level = level;

  //Catch errors
  if( f_level < (-totalLevels) )
  {
    std::cerr << "WARNING: f_level (" << f_level << ") below min allowed ("
        << -totalLevels << "). Resetting to 0." << std::endl;
    f_level = 0;
  }
  if( f_level > totalLevels )
  {
    std::cerr << "WARNING: f_level (" << f_level << ") above max allowed ("
        << totalLevels << "). Setting to 0." << std::endl;
    f_level = 0;
  }

  //Set minimum level counting from top and bottom
  if( f_level == 0 )
    level = default_level;
  else if( f_level <  0 )
    level = totalLevels + 1 + f_level;
  else
    level = f_level;

  return level;
}

template <unsigned int DIMENSION >
void
itk::PyramidScheduleTuner<DIMENSION>::
GetCroppedMinMax( unsigned int &minLevel, unsigned int &maxLevel ) const
{
  int totalLevels = m_FixedSchedule.rows();

  //Make Local copies of level min and max
  int f_LevelMin = m_LevelMin;
  int f_LevelMax = m_LevelMax;

  //MINIMUM level
  minLevel = 0;
  //Catch errors
  if( f_LevelMin < (-totalLevels) )
  {
    std::cerr << "WARNING: f_LevelMin (" << f_LevelMin << ") below min allowed ("
        << -totalLevels << "). Resetting to 0." << std::endl;
    f_LevelMin = 0;
  }
  if( f_LevelMin > totalLevels )
  {
    std::cerr << "WARNING: f_LevelMin (" << f_LevelMin << ") above max allowed ("
        << totalLevels << "). Setting to 0." << std::endl;
    f_LevelMin = 0;
  }


  //Set minimum level counting from top and bottom
  if( f_LevelMin == 0 )
    minLevel = 1;
  else if( f_LevelMin <  0 )
    minLevel = totalLevels + 1 + f_LevelMin;
  else
    minLevel = f_LevelMin;

  //MAXIMUM level
  maxLevel = 0;
  //Catch errors
  if( f_LevelMax < (-totalLevels) )
  {
    std::cerr << "WARNING: f_LevelMax (" << f_LevelMax << ") below min allowed ("
        << -totalLevels << "). Resetting to 0." << std::endl;
    f_LevelMax = 0;
  }
  if( f_LevelMax > totalLevels )
  {
    std::cerr << "WARNING: f_LevelMax (" << f_LevelMax << ") above max allowed ("
        << totalLevels << "). Setting to 0." << std::endl;
    f_LevelMax = 0;
  }

  //Set maximum level counting from top and bottom
  if( f_LevelMax == 0 )
    maxLevel = totalLevels;
  else if( f_LevelMax <  0 )
    maxLevel = totalLevels + 1 + f_LevelMax;
  else
    maxLevel = f_LevelMax;

  //Check min is not above max
  if( minLevel > maxLevel )
  {
    std::cerr << "WARNING: maxLevel (" << maxLevel << ") below maxLevel ("
        << minLevel << "). Setting to minLevel." << std::endl;
    maxLevel = minLevel;
  }

}

template <unsigned int DIMENSION >
void
itk::PyramidScheduleTuner<DIMENSION>::Update()
{
#ifdef DRH_USE_OLD_STUFF
  this->UpdateVintage();
#else
  this->UpdateSymmetric();
  //this->UpdateSimple();
#endif
}

template <unsigned int DIMENSION >
void
itk::PyramidScheduleTuner<DIMENSION>::UpdateVintage() //DRH - To delete once we're happy with UpdateSymmetric()
{
  //Set up initial shrink factors, so that shrinkage gives minSize^3 volume at
  //highest pyramid level

  //First get min_images of max_dim dimension in images
  unsigned int longestMovingLength =
      *std::max_element( m_MovingSize.GetSize(), m_MovingSize.GetSize() + DIMENSION );
  unsigned int longestFixedLength =
      *std::max_element( m_FixedSize.GetSize(), m_FixedSize.GetSize() + DIMENSION );
  unsigned int longestLength = std::min(longestMovingLength,longestFixedLength);

  //The highest possible octave n in the image pyramid,
  //such that longestlength * 0.5^n is greater than m_MinLength
  const double maxPossibleLevelD =
      std::max(1.0, //In case image with less than m_MinLength occurs
          std::floor( std::log(double(longestLength/m_MinLength))/std::log(2.0) )+1.0 );
  const unsigned int maxPossibleLevel =
      static_cast<unsigned int>( maxPossibleLevelD );

  //Extract the schedule so we can modify it appropriately.
  //This is a 2D matrix which is levels by dimensions in size
  //Access using schedule[level][dim]
  m_MovingSchedule.SetSize( maxPossibleLevel, DIMENSION );
  m_MovingSchedule.Fill(1);

  //Check of moving image has hhigher resolution than fixed image.
  //If so, lower resolution appropriately
  for( unsigned int dim = 0; dim != DIMENSION; ++dim )
  {
    int levels = std::max( 0, static_cast<int>(
        std::log(m_FixedSpacing[dim]/m_MovingSpacing[dim])/std::log(2.0) ) );
    m_MovingSchedule[maxPossibleLevel-1][dim] <<= levels;
  }

  //Initial schedule
  for( unsigned int dim = 0; dim != DIMENSION; ++dim )
    for( int level=maxPossibleLevel-2; level>=0; --level )
      m_MovingSchedule[level][dim] = m_MovingSchedule[level+1][dim]*2;

  //For each dimension specify the sampling rate at each level
  for( unsigned int dim = 0; dim != DIMENSION; ++dim )
  {
    //Get the length of current dimension at base of pyramid
    //const unsigned int length = m_MovingSize[dim]; //DRH
    const unsigned int length = std::min(m_MovingSize[dim], m_FixedSize[dim]); //DRH

    //Get the max level
    const double maxLevelD =
        std::max(1.0, std::floor( std::log(double(length/m_MinLength))/std::log(2.0) )+1.0 );
    const unsigned int maxLevel = static_cast<unsigned int>( maxLevelD );

    //Clip schedule to ensure that max length is never exceeded Note
    //ordering is in descending order of level, i.e. coarsest sampling
    //first
    for( unsigned int level=maxLevel; level<maxPossibleLevel; ++level )
      m_MovingSchedule[maxPossibleLevel-level-1][dim] =
          m_MovingSchedule[maxPossibleLevel-maxLevel][dim];
  }

  //Get appropriate levels in fixed pyramid:
  //i.e. the spacing should be larger at level but no more than twice the
  //moving spacing
  m_FixedSchedule = m_MovingSchedule;

  //For each dimension specify the sampling rate at each level
  for( unsigned int dim = 0; dim != DIMENSION; ++dim )
  {
    //Get the length of current dimension at base of pyramid
    const double movingBaseSpacing =  m_MovingSpacing[dim];
    const double fixedBaseSpacing =  m_FixedSpacing[dim];

    //Get max sample rate for current dimension
    //const unsigned int fixedLength = m_FixedSize[dim];
    //const double maxLevelD = std::max(0.0, //In case image with less than m_MinLength occurs
    //                                  std::floor( std::log(fixedLength/m_MinLength)/std::log(2) )+0.0 );
    //const double maxFixedSampleRateD = std::pow( 2.0, maxLevelD );
    //const unsigned int maxFixedSampleRate = static_cast<unsigned int>(maxFixedSampleRateD);

    //Clip schedule to ensure this length is never exceeded
    for( unsigned int level=0; level!=maxPossibleLevel; ++level )
    {
      //Calculate the current moving Spacing
      const double movingSampleRate = m_MovingSchedule[level][dim];
      const double movingSpacing = movingBaseSpacing * movingSampleRate;

      //Calculate how many doublings of fixed spacing are needed to get equivalent
      //ceiling is used because we want a LOWER sampling rate
      double log2diff = std::log( movingSpacing / fixedBaseSpacing ) / std::log(2.0);
      const double fixedSampleRateD = std::pow( 2.0, std::floor( 0.5 + log2diff ) );
      //Ensure we never get values less than 1 or enough to decrease image size by less than 16
      unsigned int fixedSampleRate = static_cast<unsigned int>( std::max( 1.0, fixedSampleRateD ) );

      //Rely on moving sample rate to limit amount of downsampling
      m_FixedSchedule[level][dim] = fixedSampleRate; //std::min( maxFixedSampleRate, fixedSampleRate);
    }

    if( 0 == m_LevelMin ) m_LevelMin = 1;
    if( 0 == m_LevelMax ) m_LevelMax = m_FixedSchedule.rows();
  }
}

// Return the number of time 'n' that 'base' can be doubled while still being smaller than dimension
// Formally: return the max 'n' subject to ('base' * 2^'n' <= 'dimension')
inline double __computeNumberOfOctaves(double dimension, double base)
{
  return std::floor( std::log(double(dimension/base))/std::log(2.0) );
}

// Return 1 + max(0, the number of octave (see above))
inline unsigned int __computeNumberOfLevels(double dimension, double base)
{
  return 1u + static_cast<unsigned int>(std::max(0., __computeNumberOfOctaves(dimension, base)));
}

// Return the best sampling rate to align the smallSpacing data with the bigSpacing data
// Note assume that bigSpacing > smallSpacing
inline unsigned int __computeBestSamplingRate(double bigSpacing, double smallSpacing)
{
  const double MAGIC = 4./3.; // Chosen such that ||bigSpacing - smallSpacing * 2^exponent|| is minimal
  unsigned int exponent = static_cast<unsigned int>(std::max(0., __computeNumberOfOctaves(bigSpacing*MAGIC, smallSpacing)));
  return std::pow<unsigned int>(2, exponent);
}

typedef itk::Array2D<unsigned int> __ScheduleType;
inline void __debugPrintSchedules(__ScheduleType scheduleMoving, __ScheduleType scheduleFixed) {
  const unsigned int NLEVELS   = scheduleMoving.rows();
  const unsigned int DIMENSION = scheduleMoving.cols();

  std::cout << "# DRH - Level" << std::setw(18) << "ScheduleMoving" << std::setw(18) << "Schedule Fixed" << std::endl;
  for (unsigned int level=0; level<NLEVELS; level++) {
    std::cout << "# " << std::setw(9) << level << "        ";
    for( unsigned int dim = 0; dim != DIMENSION; ++dim ) {
      std::cout << std::setw(3) << scheduleMoving[level][dim];
    }
    std::cout << "        ";
    for( unsigned int dim = 0; dim != DIMENSION; ++dim ) {
      std::cout << std::setw(3) << scheduleFixed [level][dim];
    }
    std::cout << "\n";
  }
}

// DRH, 2013-08-29
// Member function itk::PyramidScheduleTuner<DIMENSION>::UpdateSymmetric()
//
// Set up shrink factors, so that shrinkage gives approx. m_MinLength^3 volume at
// highest pyramid level
//
// input:  m_MovingSize, m_FixedSize, m_FixedSpacing, m_MovingSpacing, m_MinLength
// output: m_MovingSchedule, m_FixedSchedule, m_LevelMin, m_LevelMax
template <unsigned int DIMENSION >
void
itk::PyramidScheduleTuner<DIMENSION>::UpdateSymmetric()
{
  unsigned int longestLengthMoving = *std::max_element( m_MovingSize.GetSize(), m_MovingSize.GetSize() + DIMENSION );
  unsigned int longestLengthFixed  = *std::max_element( m_FixedSize.GetSize(),  m_FixedSize.GetSize()  + DIMENSION );
  const unsigned int maxPossibleLevelMoving = __computeNumberOfLevels(longestLengthMoving, m_MinLength);
  const unsigned int maxPossibleLevelFixed  = __computeNumberOfLevels(longestLengthFixed,  m_MinLength);
  const unsigned int maxPossibleLevel = 1 + std::max(maxPossibleLevelMoving, maxPossibleLevelFixed);

  ScheduleType scheduleMoving, scheduleFixed; // Schedules specify sampling rates e.g. 1, 2, 4, 8...
  scheduleMoving.SetSize( maxPossibleLevel, DIMENSION );
  scheduleFixed.SetSize( maxPossibleLevel, DIMENSION );
  scheduleMoving.Fill(-1);
  scheduleFixed.Fill(-1);

  for( unsigned int dim = 0; dim != DIMENSION; ++dim )
  {
    unsigned int rIdx = 1;
    // Highest available resolution
    scheduleMoving[maxPossibleLevel-rIdx][dim] = 1u;
    scheduleFixed [maxPossibleLevel-rIdx][dim] = 1u;

    // Second highest resolution -- we perform a crude resolution alignment here
    const double movingBaseSpacing = m_MovingSpacing[dim];
    const double fixedBaseSpacing  = m_FixedSpacing[dim];

    unsigned int maxExponentMoving = std::pow<unsigned int>(2u, __computeNumberOfLevels(m_MovingSize[dim], m_MinLength) - 1u);
    unsigned int maxExponentFixed  = std::pow<unsigned int>(2u, __computeNumberOfLevels(m_FixedSize [dim], m_MinLength) - 1u);

    rIdx = 2;
    if (movingBaseSpacing < fixedBaseSpacing) {
      const unsigned int samplingRate = std::min(__computeBestSamplingRate(fixedBaseSpacing, movingBaseSpacing), maxExponentMoving);
      unsigned int s = (samplingRate>=2u) ? 2u : 1u;
      for (;s<=samplingRate; s*=2) {
        scheduleMoving[maxPossibleLevel-rIdx][dim] = s;
        scheduleFixed [maxPossibleLevel-rIdx][dim] = 1u;
        rIdx++;
      }
    } else if (movingBaseSpacing > fixedBaseSpacing) {
      unsigned int samplingRate = std::min(__computeBestSamplingRate(movingBaseSpacing, fixedBaseSpacing), maxExponentFixed);
      unsigned int s = (samplingRate>=2u) ? 2u : 1u;
      for (;s<=samplingRate; s*=2) {
        scheduleMoving[maxPossibleLevel-rIdx][dim] = 1u;
        scheduleFixed [maxPossibleLevel-rIdx][dim] = s;
        rIdx++;
      }
    } else { // Necessary for strict symmetry, will be removed later...
      scheduleMoving[maxPossibleLevel-rIdx][dim] = 1u;
      scheduleFixed [maxPossibleLevel-rIdx][dim] = 1u;
      rIdx++;
    }

    // Other levels: we simply double the previous one until we reach the {moving, fixed} max possible levels
    for( int level=maxPossibleLevel-rIdx; level>=0; --level ) {
      // Try to increase the sampling rate
      if (2 * scheduleMoving[level+1][dim] > maxExponentMoving || 2 * scheduleFixed[level+1][dim] > maxExponentFixed) {
        scheduleMoving[level][dim] = scheduleMoving[level+1][dim];
        scheduleFixed [level][dim] = scheduleFixed[level+1][dim];
      } else {
        scheduleMoving[level][dim] = 2 * scheduleMoving[level+1][dim];
        scheduleFixed [level][dim] = 2 * scheduleFixed[level+1][dim];
      }
    }
  }

  // Find the first level (level 1 and 2 might be duplicate of each other at this point)
  unsigned int firstLevel = 0;
  {
    int ndiff = 0;
    for( unsigned int dim = 0; dim != DIMENSION; ++dim ) {
      if (scheduleMoving[maxPossibleLevel-1][dim] != scheduleMoving[maxPossibleLevel-2][dim]) {ndiff++;}
      if (scheduleFixed [maxPossibleLevel-1][dim] != scheduleFixed [maxPossibleLevel-2][dim]) {ndiff++;}
    }

    if (ndiff == 0) {
      firstLevel = 1;
    }
  }

  // Find the last _*_significant_*_ levels : [|M(i) - M(i-1)| + |F(i) - F(i-1)|] > 3
  unsigned int lastLevel=maxPossibleLevel;
  for (unsigned int level=2; level<maxPossibleLevel; level++) {
    int ndiff = 0;
    for( unsigned int dim = 0; dim != DIMENSION; ++dim ) {
      if (scheduleMoving[maxPossibleLevel -level -1][dim] != scheduleMoving[maxPossibleLevel -level][dim]) {ndiff++;}
      if (scheduleFixed [maxPossibleLevel -level -1][dim] != scheduleFixed [maxPossibleLevel -level][dim]) {ndiff++;}
    }

    if (ndiff < 3) {
      lastLevel = level;
      break;
    }
  }

  // Apply user-selectable bound
  unsigned int trueNumberOfLevel = lastLevel - firstLevel;
  if (trueNumberOfLevel > m_MaxLevelNumber) {
    trueNumberOfLevel = m_MaxLevelNumber;
  }

  // Store the final schedules (non-redundant levels, starting from the highest resolution)
  m_MovingSchedule.SetSize( trueNumberOfLevel, DIMENSION );
  m_FixedSchedule.SetSize( trueNumberOfLevel, DIMENSION );

  for (unsigned int level=0; level<trueNumberOfLevel; level++) {
    for( unsigned int dim = 0; dim != DIMENSION; ++dim ) {
      m_MovingSchedule[trueNumberOfLevel -level -1][dim] = scheduleMoving[maxPossibleLevel -firstLevel -level -1][dim];
      m_FixedSchedule [trueNumberOfLevel -level -1][dim] = scheduleFixed [maxPossibleLevel -firstLevel -level -1][dim];
    }
  }
  //__debugPrintSchedules(m_MovingSchedule, m_FixedSchedule);

  if( 0 == m_LevelMin ) m_LevelMin = 1;
  if( 0 == m_LevelMax ) m_LevelMax = m_FixedSchedule.rows();
}


// DRH, 2013-08-29
// Member function itk::PyramidScheduleTuner<DIMENSION>::UpdateSimple()
//
// Implement the simplest pyramid schedule: just dived the image size by 2 at each level
// The maximal number of level == level at wich the longuest image dim is just bigger than m_MinLength (default: 32)
//
// input:  m_MovingSize, m_FixedSize, m_FixedSpacing, m_MovingSpacing, m_MinLength
// output: m_MovingSchedule, m_FixedSchedule, m_LevelMin, m_LevelMax
template <unsigned int DIMENSION >
void
itk::PyramidScheduleTuner<DIMENSION>::UpdateSimple()
{
  unsigned int longestLengthMoving = *std::max_element( m_MovingSize.GetSize(), m_MovingSize.GetSize() + DIMENSION );
  unsigned int longestLengthFixed  = *std::max_element( m_FixedSize.GetSize(),  m_FixedSize.GetSize()  + DIMENSION );
  const unsigned int maxPossibleLevelMoving = __computeNumberOfLevels(longestLengthMoving, m_MinLength);
  const unsigned int maxPossibleLevelFixed  = __computeNumberOfLevels(longestLengthFixed,  m_MinLength);
  unsigned int maxPossibleLevel = std::max(maxPossibleLevelMoving, maxPossibleLevelFixed);

  // Apply user-selectable bound
  if (maxPossibleLevel > m_MaxLevelNumber) {
    maxPossibleLevel = m_MaxLevelNumber;
  }

  m_MovingSchedule.SetSize( maxPossibleLevel, DIMENSION );
  m_FixedSchedule.SetSize( maxPossibleLevel, DIMENSION );

  for( unsigned int dim = 0; dim != DIMENSION; ++dim )
  {
    // Highest available resolution
    m_MovingSchedule[maxPossibleLevel-1][dim] = 1u;
    m_FixedSchedule [maxPossibleLevel-1][dim] = 1u;

    unsigned int maxExponentMoving = std::pow<unsigned int>(2u, __computeNumberOfLevels(m_MovingSize[dim], m_MinLength) - 1u);
    unsigned int maxExponentFixed  = std::pow<unsigned int>(2u, __computeNumberOfLevels(m_FixedSize [dim], m_MinLength) - 1u);

    // Other levels: we simply double the previous one until we reach the {moving, fixed} max possible levels
    for( int level=maxPossibleLevel-2; level>=0; --level ) {
      m_MovingSchedule[level][dim] = std::min(2 * m_MovingSchedule[level+1][dim], maxExponentMoving);
      m_FixedSchedule [level][dim] = std::min(2 * m_FixedSchedule[level+1][dim], maxExponentFixed);
    }
  }
  //__debugPrintSchedules(scheduleMoving, scheduleFixed);

  if( 0 == m_LevelMin ) m_LevelMin = 1;
  if( 0 == m_LevelMax ) m_LevelMax = m_FixedSchedule.rows();
}


namespace
{

template <typename ATYPE, typename BTYPE>
vnl_vector<double>
vnl_multiply_elements_(
    const vnl_vector<ATYPE> & aa,
    const vnl_vector<BTYPE> & bb
)
{
  if( aa.size() != bb.size() )
    throw( std::length_error("Input vectors are different sizes.") );

  vnl_vector<double> cc(aa);
  for( unsigned int ii=0; ii<cc.size(); ++ii )
    cc[ii] *= bb[ii];

  return cc;
}

template <unsigned int ASIZE, typename BTYPE>
vnl_vector<double>
vnl_divide_elements_(
    const itk::Size<ASIZE> & aa,
    const vnl_vector<BTYPE> & bb
)
{
  if( ASIZE != bb.size() )
    throw( std::length_error("Input vectors are different sizes.") );

  vnl_vector<double> cc(ASIZE);
  for( unsigned int ii=0; ii<ASIZE; ++ii )
    cc[ii] = aa[ii] / bb[ii];

  return cc;
}
}

template <unsigned int DIMENSION >
void
itk::PyramidScheduleTuner<DIMENSION>::PrintPyramidInfo(std::ostream& os, int verbosity) const
{
  unsigned int min_level = 0, max_level = 0;
  GetCroppedMinMax( min_level, max_level );

  if (verbosity >= 2) {
  os << "Pyramid levels output (specified):"
      <<"\n  MIN "<<min_level<<"("<<m_LevelMin<<"). Note: 0=>min -1=>top 1=>bottom."
      <<"\n  MAX "<<max_level<<"("<<m_LevelMax<<"). Note: 0=>max -1=>top 1=>bottom."
      << std::endl;
  }

  unsigned int levels = m_MovingSchedule.rows();
  os << "Limits  Level   Spacing                                              Sampling rate                Image size"
     << std::endl;

  for( unsigned int level = 0; level < levels; ++level )
  {
    if(min_level==max_level && min_level==level+1)
      os << "1 Lvl->" ;//<< std::endl;
    else
    {
      if(min_level==level+1) {
        os << "1st  ->" ;
      } else if(max_level==level+1) {
        os << "Last ->" ;
      } else {
        os << "       " ;
      }
    }
#if ITK_VERSION_MAJOR < 4
    os
    << "\t" << level+1 << "   "
    << std::setw(5) << std::setprecision(3)
    << "m[" << ::vnl_multiply_elements_<double,unsigned int>( //Moving Spacing
        m_MovingSpacing.Get_vnl_vector(),
        m_MovingSchedule.get_row(level) ) << "] "
    << "f["<< ::vnl_multiply_elements_<double,unsigned int>( //Fixed spacing
        m_FixedSpacing.Get_vnl_vector(),
        m_FixedSchedule.get_row(level) ) << "]  "
    << "m[" << m_MovingSchedule.get_row(level) << "] " //Moving sample rate
    << "f["<< m_FixedSchedule.get_row(level) << "]  " //Fixed sample rate
    << "m[" << ::vnl_divide_elements_<3,unsigned int>( //Moving Size
        m_MovingSize,
        m_MovingSchedule.get_row(level) ) << "] "
    << "f[" << ::vnl_divide_elements_<3,unsigned int>( //Fixed size
        m_FixedSize,
        m_FixedSchedule.get_row(level) ) << "]"
    << std::endl;
#else
    std::streamsize s = os.precision();
    os
    << " " << level+1 << "   "
    << std::setw(5) << std::setprecision(3)
    << "m" << PrettyPrint(::vnl_multiply_elements_<double,unsigned int>( //Moving Spacing
        m_MovingSpacing.GetVnlVector(),
        m_MovingSchedule.get_row(level) ), 6, " ") << " "
    << "f"<< PrettyPrint(::vnl_multiply_elements_<double,unsigned int>( //Fixed spacing
        m_FixedSpacing.GetVnlVector(),
        m_FixedSchedule.get_row(level) ), 6, " ") << "  "
    << "m" << PrettyPrint(m_MovingSchedule.get_row(level), 2, " ") << " " //Moving sample rate
    << "f"<< PrettyPrint(m_FixedSchedule.get_row(level), 2, " ") << "  " //Fixed sample rate
    << std::setprecision(0)
    << "m" << PrettyPrint(::vnl_divide_elements_<3,unsigned int>( //Moving Size
        m_MovingSize,
        m_MovingSchedule.get_row(level) ), 3, " ") << " "
    << "f" << PrettyPrint(::vnl_divide_elements_<3,unsigned int>( //Fixed size
        m_FixedSize,
        m_FixedSchedule.get_row(level) ), 3, " ") << ""
    << std::setprecision(s)
    << std::endl;
#endif
  }
}

template <unsigned int DIMENSION >
void
itk::PyramidScheduleTuner<DIMENSION>::PrintSelf
(std::ostream& os, Indent /*indent*/) const
{
  this->PrintPyramidInfo(os,1);
}

#endif /* ITKPYRAMIDSCHEDULETUNER_TXX_ */
