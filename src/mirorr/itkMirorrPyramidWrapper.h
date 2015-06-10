/*=========================================================================
Program: mirorr
Module: itkMirorrPyramidWrapper.h
Author: Nicholas Dowson
Created: 20 Feb 2009

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
// History:
// Revised: 13 Jul 2009 by N Dowson. Modification to set translation
//   of translation using set offset
// Revised: 17 Jun 2009 by N Dowson. Added ability to write out ITK
//   transform file

#ifndef __MirorrPyramidWrapper_h
#define __MirorrPyramidWrapper_h

#include <sstream>
#include <string>
#include <boost/timer/timer.hpp>
#include <itkTransformFileWriter.h>
#include <itkTransformFileReader.h>
#include <itkContinuousIndex.h>
#include <itkImageFileReader.h>

#if ITK_VERSION_MAJOR < 4
  #include <itkOrientedImage.h>
#else
  #include <itkImage.h>
#endif

#include <itkTransformFactory.h>
#include <itkTransformFactory.h>
#include <itkTransformBase.h>
#include <itkResampleImageFilter.h>
#include <vnl/vnl_cross.h>
#include "itkMirorrPyramidImplement.h"

namespace itk
{

/** \class A wrapper class for the MirorrPyramidImplement class.
 * This particular class adds some useful helper functions for
 * interacting with the file system and some other things:
 *   - reading and writing transform files
 *   - reading and writing images
 *   - performing an initial alignment of the two images
 *   - reading in a transformation file, and giving the user an example
 *     if this fails
 *   - saving the transformation to a file
 */
//template <unsigned int DIMENSION >
class MirorrPyramidWrapper
{
public:
  static const unsigned int DIMENSION = 3;

  //!Some useful type definitions
  typedef MirorrPyramidImplement/*<DIMENSION>*/                      MirorrType;
  typedef MirorrType::TransformType                                  TransformType;
  typedef itk::MatrixOffsetTransformBase<double,DIMENSION,DIMENSION>  InputTransformType;

  typedef MirorrType::ParametersType ParametersType;
  typedef MirorrType::ImageType      ImageType;
  typedef MirorrType::ImagePointer   ImagePointer;
  typedef ImageType::PixelType        PixelType;
  typedef ImageType::SpacingType      SpacingType;

#if ITK_VERSION_MAJOR < 4
  typedef itk::OrientedImage<unsigned char, DIMENSION> MaskType;
#else
  typedef itk::Image<unsigned char, DIMENSION> MaskType;
#endif
  
  typedef MaskType::Pointer                            MaskPointer;

  typedef itk::Vector<double, DIMENSION> VectorType;
  typedef itk::Size<DIMENSION>           SizeType;
  typedef itk::Index<DIMENSION>          IndexType;
  typedef itk::Point<double,DIMENSION>   PointType;

  MirorrPyramidWrapper()
  { //Defaults
    do_not_register = false;
    invert_input_transform = false;
    invert_output_transform = false;
    verbosity = 1;
    do_reorient_fixed = false;
    do_reorient_moving = false;
    program_name = "mirorr";
    do_force_resample_fixed = false;
    do_force_resample_moving = false;
    do_initialise_tfm_to_centre_images = true;
  }

  MirorrType & GetRegistrationPyramidObject() { return mirorr; }

  void ReadAndResampleImages();

  /**Constructor. This takes a MirorrPyramidWrapperInputsType as an
   * input because the class designed to run when created and then
   * instantly be destroyed.  There are too many inputs to be
   * specified individually*/
  void Update();

  //!Name of file storing fixed image
  void SetFixedName( const std::string & _ss ) {fixedName = _ss; }
  //!Name of file storing moving image
  void SetMovingName( const std::string & _ss ) {movingName = _ss; }
  const std::string & GetMovingName( ) const { return movingName; }
  //!Name of file storing fixed mask image
  void SetFixedMaskName( const std::string & _ss ) {fixedMaskName = _ss; }
  const std::string & GetFixedName( ) const { return fixedName; }
  //!Name of file storing moving mask image
  void SetMovingMaskName( const std::string & _ss ) {movingMaskName = _ss; }
  //! Do we attempt to reorient the image in the ARI direction?
  void SetDoReorientFixedInARI(bool in ) {do_reorient_fixed = in;}
  void SetDoReorientMovingInARI(bool in ) {do_reorient_moving = in;}

  //!Name of file string fixed image once it has been registered to moving
  //!image. If empty - no image is returned.
  const std::string & GetLastTransformedFixedName()
  { return lastTransformedFixedName; }
  void SetLastTransformedFixedName( const std::string & _ss, bool force_resample = false )
  { do_force_resample_fixed = force_resample; lastTransformedFixedName = _ss; }

  //!Name of file string moving image once fixed image been registered to moving
  //!image. If empty - no image is returned.
  const std::string & GetLastTransformedMovingName()
  { return lastTransformedMovingName; }
  void SetLastTransformedMovingName( const std::string & _ss, bool force_resample = false )
  { do_force_resample_moving = force_resample;  lastTransformedMovingName = _ss; }

  /** Name of file storing initial transform. If empty the image
   *  centres are aligned, with no rotation, and a file with the
   *  initial output is saved. */
  void SetInitialTransformName( const std::string & _ss ) {initialTransformName = _ss; }
  //!Name of file storing final transform out. "final.tfm" by default.
  void SetFinalTransformName( const std::string & _ss ) {finalTransformName = _ss; }
  const std::string & GetFinalTransformName( ) const { return finalTransformName; }
  void SetDoInitialiseTransformToCentreImages(bool in ) {do_initialise_tfm_to_centre_images = in;}

  /**Name of transform type based on corresponding ITK class names
   * e.g. Euler3DTransform, AffineTransform.  The names
   * "rigid","affine" and "translation" are converted to their
   * respective ITK classes */
  void SetTransformType( const std::string & _ss ) {transformType = _ss; }
  const std::string & GetTransformType( ) const { return transformType; }
  //!Output moving to fixed  image instead of vice versa
  void SetInvertInputTransform( bool in ) {invert_input_transform = in; }
  //!Output moving to fixed  image instead of vice versa
  void SetInvertOutputTransform( bool in ) {invert_output_transform = in; }

  void SetVerbosity( int in ) { verbosity = in; GetRegistrationPyramidObject().SetVerbosity(in); }
  int GetVerbosity() const { return verbosity; }

  //! Should we not actually register the images and just apply transform to resample the images?
  void SetDoNotRegister( bool in ) { do_not_register = in; }

  // For pretty display
  void SetProgramName(std::string program_name) {this->program_name = program_name;}

private:
  //!Write a transform to a file using standard ITK transform file writer
  void writeParametersUsingItkTransformFileWriter( std::string tfmName, //itk::TransformBase::Pointer
      InputTransformType::Pointer tfm,
      bool invert = false );

  //!Read  a transform from a file using standard ITK transform file writer
  int readParametersUsingItkTransformFileReader( std::string tfmName,
      InputTransformType::Pointer tfm,
      bool bInvertTransform = false );

  SpacingType GetMaxSpacingForResize( ImageType::Pointer image1,
      ImageType::Pointer image2,
      itk::Size<DIMENSION> size ) const;

  //! Create a transform of the correct type given an input string (rigid, quat, affine) */
  InputTransformType::Pointer CreateAppropriateTransform( std::string inputTransformType );
  void ReadInputTransform( InputTransformType::Pointer transform,
      std::string transformFileName, std::string _fixedName,
      bool _invert_output_transform );

  //! Method to apply a mask to an image
  //void ApplyMaskToImage( ImageType::Pointer image, MaskType::Pointer mask );

  //Attributes
  MirorrType mirorr;

  std::string fixedName; //!Name of file storing fixed image
  std::string movingName; //!Name of file storing moving image
  std::string fixedMaskName; //!Name of file storing fixed mask image
  std::string movingMaskName; //!Name of file storing moving mask image
  bool do_reorient_fixed; //! Do we attempt to reorient fixed image in the ARI direction?
  bool do_reorient_moving; //! Do we attempt to reorient moving image in the ARI direction?
  bool do_force_resample_fixed;
  bool do_force_resample_moving;

  bool do_initialise_tfm_to_centre_images;

  ImagePointer movingImage;
  ImagePointer fixedImage;
  MaskPointer movingMask;
  MaskPointer fixedMask;

  //!Name of file string fixed image once it has been registered to moving
  //!image. If empty - no image is returned.
  std::string lastTransformedFixedName;
  /** Name of file storing initial transform. If empty the image
   *  centres are aligned, with no rotation, and a file with the
   *  initial output is saved. */
  std::string lastTransformedMovingName;
  /** Name of file storing initial transform. If empty the image
   *  centres are aligned, with no rotation, and a file with the
   *  initial output is saved. */
  std::string initialTransformName; //Null if none supplied
  //!Name of file storing final transform out. "final.tfm" by default.
  std::string finalTransformName;
  /**Name of transform type based on corresponding ITK class names
   * e.g. Euler3DTransform, AffineTransform.  The names
   * "rigid","affine" and "translation" are converted to their
   * respective ITK classes */
  std::string transformType;
//  /** Name of the blockmatcher metric. Abbreviations are converted
//   * to full names as follows: "nc" to "normalized_correlation",
//   * "sd" to "sum_of_squared_differences", "cr" to "correlation_ratio",
//   * "mi" to "mutual_information", "npwmi" to "npw_mutual_information". */
//  std::string blockMetricType;
//
//  bool doSaveIntermediateImages; //! Should we save images after each iteration?
  bool invert_input_transform; //!Input moving to fixed  image instead of vice versa
  bool invert_output_transform; //!Output moving to fixed  image instead of vice versa
  int verbosity; //Verbosity of algorithm
  bool do_not_register; //Do not register the images. Just use transform to resample them
  std::string program_name;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMirorrPyramidWrapper.txx"
#endif

#endif
