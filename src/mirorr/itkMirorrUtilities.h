/*=========================================================================
Program: mirorr
Module: itkMirorrUtilities.h
Author: Nicholas Dowson
Created: 11 Jul 2017

Copyright (c) 2009-17 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#ifndef MIRORRPROJECT_ITKMIRORRUTIL_H
#define MIRORRPROJECT_ITKMIRORRUTIL_H

#include <exception>
#include <vector>
#include <string>

#include <itkImage.h>
#include <itkMatrixOffsetTransformBase.h>
#include <itkAffineTransform.h>

#include <boost/algorithm/string.hpp>

namespace itk
{
namespace util
{
/**
 * Get a transform that will rigidly transform an image to the identity space
 * such that its origin is zero and its rotation matrix is the identity
 * (although allowing permutations e.g. [0 1 0; 1 0 0; 0 0 1]
 * MatrixOffsetTransformBase<double, 3,3> tfm =
 *   GetTransformFromImageToIdentitySpace(
 *     image->GetOrigin(), image->GetDirection() );
 *
 * This code is an adapted copy-paste from imgGetIdSpaceTfm
 *
 * @param origin Image Origin
 * @param direction Image Direction Cosines
 * @return Transform
 */
typedef Image<char,3> TImageType;
AffineTransform<double, 3>::Pointer
GetTransformFromImageToIdentitySpace(
        Point<double,3> origin,
        Matrix<double,3,3> direction,
        Image<char,3>::SpacingType spacing,
        Image<char,3>::SizeType size,
        bool do_centre = true
);

/**
 * Parse an input string of a list of numbers to extract and convert individual elements.
 * Several delimiters are supported including "x :,()="
 * @tparam T The element type each element in the list should be converted to
 * @param paramString The input string in format "el1,el2,el3".
 * @param max_dimension The maximum length of the input values. Default is zero.
 * @param out
 */
template <typename T>
void
parseInputToList( std::string paramString, unsigned int max_dimension, std::vector<T> & out)
{
  try {
    std::vector<std::string> strs;
    boost::split(strs, paramString, boost::is_any_of("x :,()=") );

    std::vector<double> tout(std::min<unsigned int>(strs.size(), max_dimension), 0.0);
    for( unsigned int ii=0; ii<max_dimension && ii<strs.size(); ++ii )
    {
      std::stringstream ss(strs[ii]);
      ss >> tout[ii];
    }
    out.resize(tout.size());
    if(out.size() < max_dimension)
      out.resize(max_dimension, 0);
    std::copy( tout.begin(), tout.end(), out.begin() );
  }
  catch( ... )
  {
    std::stringstream ss;
    ss<< "Failed to parse input string: \""<< paramString
      <<"\"\nExiting.\n";
    throw std::runtime_error(ss.str());
  }
}

}

}

#endif //MIRORRPROJECT_ITKMIRORRUTIL_H
