/*=========================================================================
Program: mirorr
Module: itkMirorrUtilities.cxx
Author: Nicholas Dowson
Created: 11 Jul 2017

Copyright (c) 2009-17 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include <itkImage.h>
#include <itkMatrixOffsetTransformBase.h>
#include "itkMirorrUtilities.h"
#include <vnl/vnl_det.h>
#include <itkLandmarkBasedTransformInitializer.h>
#include <itkAffineTransform.h>

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
        bool do_centre
)
{
  //Read the input file and display its orientation
  itk::Matrix<double,3,3> t_dir = direction;
  t_dir.Fill(0.0);
  for(int ii=0; ii<3; ++ii)
  {
    //Find maximum
    std::vector<double> vals(3,0);
    std::vector<bool> vsign(3,0);
    for(int jj=0; jj<3; ++jj)
    {
      vals[jj] = fabs(direction[ii][jj]);
      vsign[jj] = direction[ii][jj]>0;
    }
    int pos = std::max_element(vals.begin(),vals.end())-vals.begin();
    t_dir[ii][pos] = vsign[pos] ? 1.0 : -1.0;
  }

  //BEWARE of negative determinants - they will flip your axes!
  double det = fabs(vnl_det( t_dir.GetVnlMatrix() ));
  for( int ii=0; ii<3; ++ii )
    for( int jj=0; jj<3; ++jj )
      t_dir[ii][jj] /= det;

  itk::Matrix<double,3,3> out_dir = t_dir;

  //Create the transform
  typedef itk::Image<double,3> T_Image;
  typedef T_Image::PointType T_Point;

  //typedef itk::MatrixOffsetTransformBase<double, 3,3> T_Transform;
  typedef itk::AffineTransform<double, 3> T_Transform;
  typedef itk::LandmarkBasedTransformInitializer< T_Transform, T_Image, T_Image > T_LandmarkBasedTransformInitializer;

//Create the transform
  typedef itk::Image<double,3> T_Image;
  typedef T_Image::PointType T_Point;

  //typedef itk::MatrixOffsetTransformBase<double, 3,3> T_Transform;
  typedef itk::AffineTransform<double, 3> T_Transform;
  typedef itk::LandmarkBasedTransformInitializer< T_Transform, T_Image, T_Image > T_LandmarkBasedTransformInitializer;

  T_LandmarkBasedTransformInitializer::LandmarkPointContainer fixedLandmarks, movingLandmarks;

  T_Image::IndexType index[4];
  for( unsigned int ii=0; ii<4; ++ii ) {
    index[ii].Fill(0);
    if (ii > 0)
      index[ii][ii - 1] = 1;
  }

  //T_Image::SpacingType spacing;
  //spacing.Fill(1);
  //T_Image::SizeType size;
  //size.Fill(10);

  T_Image::RegionType region;
  region.SetSize(size);

  T_Image::Pointer img1 = T_Image::New();
  img1->SetOrigin(origin);
  img1->SetDirection(direction);
  img1->SetSpacing(spacing);
  img1->SetRegions(region);

  T_Image::Pointer img2 = T_Image::New();
  T_Image::PointType origin2;
  for(int ii=0; ii<3; ++ii)
    origin2[ii] = do_centre ? size[ii] * spacing[ii] * -0.5 : 0.0;
  img2->SetOrigin(origin2);
  img2->SetDirection(t_dir);
  img2->SetSpacing(spacing);
  img2->SetRegions(region);

  for( unsigned int ii=0; ii<4; ++ii ) {
    T_Point point;
    img1->TransformIndexToPhysicalPoint(index[ii], point);
    fixedLandmarks.push_back(point);
    img2->TransformIndexToPhysicalPoint(index[ii], point);
    movingLandmarks.push_back(point);
  }

  T_Transform::Pointer transform = T_Transform::New();

  T_LandmarkBasedTransformInitializer::Pointer tfm_maker = T_LandmarkBasedTransformInitializer::New();
  tfm_maker->SetFixedLandmarks(fixedLandmarks);
  tfm_maker->SetMovingLandmarks(movingLandmarks);
  tfm_maker->SetTransform(transform);
  tfm_maker->InitializeTransform();

  //Create the transform
//  typedef itk::MatrixOffsetTransformBase<double, 3,3> T_Transform;
//  T_Transform::Pointer transform = T_Transform::New();
//
//  T_Transform::InputPointType centre;
//  centre.Fill(0.0);
//  transform->SetCenter(centre);
//
//  T_Transform::OffsetType offset;
//  for(unsigned int ii=0; ii<3; ++ii)
//    offset[ii] = 0 - origin[ii];
//  transform->SetOffset(offset);
//
//  for(unsigned int ii=0; ii<3; ++ii)
//    offset[ii] = 0 - origin[ii];
//
//  T_Transform::MatrixType matrix;
//  //matrix = out_dir.GetInverse() * direction.GetVnlMatrix(); //Working
//  matrix = out_dir.GetVnlMatrix() * direction.GetInverse();
//  transform->SetMatrix(matrix);

  return transform;
}

}

}

