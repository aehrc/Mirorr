/*=========================================================================
Program: mirorr
Module: MirorrPyramidWrapper.txx
Author: Nicholas Dowson
Created: 1 Oct 2010

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#ifndef ITKMIRORRPYRAMIDWRAPPER_TXX_
#define ITKMIRORRPYRAMIDWRAPPER_TXX_

#include "itkMirorrPyramidWrapper.h"
#include <itkOrientImageFilter.h>
#include <itkCenteredTransformInitializer.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include "itkIOUtils.h"
#include <iomanip>

namespace itk {

//An anonymous namespace to store important routines
namespace __MirorrPyramidWrapper {
//!Read an image file
template<typename TImageType>
typename TImageType::Pointer
ReadImage(std::string image_file_name) {
  //Read in the moving image
  typename itk::ImageFileReader<TImageType>::Pointer imageFileReader =
          itk::ImageFileReader<TImageType>::New();

  imageFileReader->SetFileName(image_file_name.c_str());
  try {
    imageFileReader->Update();
  } catch (itk::ExceptionObject &err) {
    std::cerr << "ExceptionObject caught:" << std::endl;
    std::cerr << err << std::endl;
    std::cout << err.what() << std::endl;
    exit(1);
  }

  typename TImageType::Pointer image = imageFileReader->GetOutput();

  return image;
}

//typedef itk::Image<double, 3> TImageType;
template< typename TImageType >
std::ostream & PrintImage( std::ostream & os, typename TImageType::Pointer image )
{
  typename TImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
  typename TImageType::PointType origin = image->GetOrigin();
  typename TImageType::SpacingType spacing = image->GetSpacing();
  typename TImageType::DirectionType dir = image->GetDirection();
  os<<"<Image> size=["<<size[0]<<" "<<size[1]<<" "<<size[2]<<"], "
    <<"spacing=["<<std::setprecision(3)<<spacing[0]<<" "<<spacing[1]<<" "<<spacing[2]<<"], "
    <<"origin=["<<std::setprecision(5)<<origin[0]<<" "<<origin[1]<<" "<<origin[2]<<"], "
    <<"dir=["<<std::setprecision(3)<<dir[0][0]<<" "<<dir[0][1]<<" "<<dir[0][2]<<"; "
    <<dir[1][0]<<" "<<dir[1][1]<<" "<<dir[1][2]<<"; "
    <<dir[2][0]<<" "<<dir[2][1]<<" "<<dir[2][2]<<"];";
  return os;
}

template<typename TImageType>
typename TImageType::Pointer
ReorientARIImage(typename TImageType::Pointer image) {
  typedef typename itk::OrientImageFilter<TImageType, TImageType> TOrientFilter;

  typename TOrientFilter::Pointer orienter = TOrientFilter::New();
  orienter->UseImageDirectionOn();
  orienter->SetDesiredCoordinateOrientation(
      itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARI);
  orienter->SetInput(image);
  orienter->Update();

  //=======================================================
  //ND 16 April 2015: Mod to reset direction to identity to
  //avoid adverse influence of interpolation artefacts when images are
  //almost but not quite aligned
  typename TImageType::Pointer out_image = orienter->GetOutput();

  typename TImageType::DirectionType dir;
  dir.Fill(0);
  dir[0][1] = 1;
  dir[1][0] = 1;
  dir[2][2] = 1;
  out_image->SetDirection(dir);

  typename TImageType::PointType origin;
  origin[0] = 0; origin[1] = 0; origin[2] = 0;
  out_image->SetOrigin(origin);
  return out_image;
  //=======================================================

  //return orienter->GetOutput();
}

template<typename TMatrix>
vnl_vector<double>
ComputeRotationVector(TMatrix &mat) {
  // see: http://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation
  double theta = acos((vnl_trace(mat.GetVnlMatrix()) - 1.) / 2.);
  const double EPS = 1e-16;

  vnl_matrix<double> N = (mat.GetVnlMatrix() - mat.GetTranspose()) / (2.0 * sin(theta) + EPS);

  vnl_vector<double> n(3);
  n(0) = N(2, 1);
  n(1) = N(0, 2);
  n(2) = N(1, 0);

  return theta * n;
}

vnl_matrix<double>
ComputeRotationMatrixFromRotationVector(vnl_vector<double> rotvec) {
  // see: http://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation
  double theta = rotvec.two_norm();
  vnl_vector<double> &w = rotvec.normalize();

  vnl_matrix<double> I(3, 3);
  I.set_identity();
  vnl_matrix<double> W(3, 3, 0.); //cross-product matrix of w
  W(0, 0) = 0.;
  W(0, 1) = -w(2);
  W(0, 2) = w(1);
  W(1, 0) = w(2);
  W(1, 1) = 0.;
  W(1, 2) = -w(0);
  W(2, 0) = -w(1);
  W(2, 1) = w(0);
  W(2, 2) = 0.;

  return I + W * sin(theta) + W * W * (1. - cos(theta));
}

template<typename TMatrix>
double
ComputeDirectionDifference(TMatrix &matA, TMatrix &matB) {
  vnl_vector<double> rotvecA = ComputeRotationVector(matA);
  vnl_vector<double> rotvecB = ComputeRotationVector(matB);
  //DRH: rotation vector will not convert back to the same matrix if the input matrix is not in SO3

  rotvecA -= rotvecB;
  return rotvecA.two_norm();
}

} //namespace


void
MirorrPyramidWrapper
::ReadAndResampleImages() {
  //Read in images
  ImagePointer movingImage_in = __MirorrPyramidWrapper::ReadImage<ImageType>(movingName);
  ImagePointer fixedImage_in = __MirorrPyramidWrapper::ReadImage<ImageType>(fixedName);
  MaskPointer movingMask_in;
  MaskPointer fixedMask_in;

  //Read or create masks
  if (!movingMaskName.empty()) {
    movingMask_in = __MirorrPyramidWrapper::ReadImage<MaskType>(movingMaskName);
    //ApplyMaskToImage( movingImage_in, movingMask_in );
  }
  else {
    movingMask_in = __MirorrPyramidWrapper::ReadImage<MaskType>(movingName);
    movingMask_in->FillBuffer(255);
  }

  if (!fixedMaskName.empty()) {
    fixedMask_in = __MirorrPyramidWrapper::ReadImage<MaskType>(fixedMaskName);
    //ApplyMaskToImage( fixedImage_in, fixedMask_in );
  }
  else {
    fixedMask_in = __MirorrPyramidWrapper::ReadImage<MaskType>(fixedName);
    fixedMask_in->FillBuffer(255);
  }

  movingImage = movingImage_in;
  fixedImage = fixedImage_in;
  movingMask = movingMask_in;
  fixedMask = fixedMask_in;

  // Try to re-orient the images
  if (this->do_reorient) {
    std::cerr << "#\n# Warning: Using the experimental --reorient feature.\n";
    std::cerr << "#          Carefully inspecting the result is strongly advised.\n#" << std::endl;

    movingImage = __MirorrPyramidWrapper::ReorientARIImage<ImageType>(movingImage);
    fixedImage = __MirorrPyramidWrapper::ReorientARIImage<ImageType>(fixedImage);
    movingMask = __MirorrPyramidWrapper::ReorientARIImage<MaskType>(movingMask);
    fixedMask = __MirorrPyramidWrapper::ReorientARIImage<MaskType>(fixedMask);
  }

  //Check if the two images have the same direction matrix and issue a warning if not
  const double DIRECTION_THRESHOLD = 0.017453; //== 1 deg
  const double delta = __MirorrPyramidWrapper::ComputeDirectionDifference(movingImage->GetDirection(),
                                                                          fixedImage->GetDirection());
  if (delta > DIRECTION_THRESHOLD) {
    std::cerr << "#\n# Warning: The moving and fixed images have different cosine direction matrices.\n";
    std::cerr << "#          Low quality results are to be expected. Consider using --reorient.\n";
    std::cerr << "#          Note: DELTA = " << delta << ".\n#" << std::endl;
  }
}

//!Write a transform to a file using standard ITK transform file writer
void
MirorrPyramidWrapper::
writeParametersUsingItkTransformFileWriter(
        std::string tfmName, //itk::TransformBase::Pointer
        InputTransformType::Pointer tfm,
        bool invert) {
  InputTransformType::Pointer tfm2 = tfm;
  //The inversion here means output is Rigid3DTransform_double_3_3 - which causes problems later
  if (invert) {
    InputTransformType::Pointer tfm_inverse =
            dynamic_cast<InputTransformType *>( tfm->GetInverseTransform().GetPointer());

    typedef AffineTransform<double, DIMENSION> AffineTransformType;
    AffineTransformType::Pointer tfm_affine = AffineTransformType::New();
    tfm_affine->SetTranslation(tfm_inverse->GetTranslation());
    tfm_affine->SetMatrix(tfm_inverse->GetMatrix());
    tfm2 = dynamic_cast<InputTransformType *>( tfm_affine.GetPointer());
  }

  //Write the transform to a file
  typedef itk::TransformFileWriter TransformWriterType;
  TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  transformWriter->SetFileName(tfmName.c_str());
  transformWriter->SetInput(tfm2);
  try {
    transformWriter->Update();
  } catch (itk::ExceptionObject &err) {
    std::cerr << "ExceptionObject caught:" << std::endl;
    std::cerr << err << std::endl;
    std::cout << err.what() << std::endl;
  }

}

//!Read  a transform from a file using standard ITK transform file writer
int
MirorrPyramidWrapper::
readParametersUsingItkTransformFileReader(
        std::string tfmName,
        InputTransformType::Pointer tfm,
        bool bInvertTransform) {
  //Write the transform to a file
  typedef itk::TransformFileReader TransformReaderType;
  TransformReaderType::Pointer transformReader = TransformReaderType::New();
  try {
    //Read the transform
    transformReader->SetFileName(tfmName.c_str());
    transformReader->Update();
    TransformReaderType::TransformPointer readTransform;
    if (transformReader->GetTransformList()->empty())
      return -1;
    readTransform = transformReader->GetTransformList()->front();

    InputTransformType::Pointer affineTransform;

    //Check input transform is affine or
    std::string className = readTransform->GetNameOfClass();
    if (className.find("Affine") != std::string::npos ||
        className.find("Euler") != std::string::npos) {
      affineTransform = dynamic_cast<InputTransformType *>( readTransform.GetPointer());
    }
    else {
      std::stringstream ss;
      ss << "Error: Input transform file '" << tfmName
                                               << "' does not contain a MatrixOffsetTransform type";
      throw(std::runtime_error(ss.str()));
    }

    if (bInvertTransform)
      affineTransform = dynamic_cast<InputTransformType *>(
              affineTransform->GetInverseTransform().GetPointer());

    //Output the transform
    //Ensure input is orthogonal for Euler transform
    if (className.find("Euler") != std::string::npos) {
      vnl_matrix_fixed<double, 3, 3> matrix =
              affineTransform->GetMatrix().GetVnlMatrix();
      matrix.set_row(2, vnl_cross_3d(
              matrix.get_row(0),
              matrix.get_row(1)));
      InputTransformType::MatrixType tmatrix(matrix);

      affineTransform->SetMatrix(tmatrix);
    }

    tfm->SetMatrix(affineTransform->GetMatrix());
    tfm->SetTranslation(affineTransform->GetTranslation());

  } catch (itk::ExceptionObject &err) {
    std::cerr << "Failed to read " << tfmName << " using as ITK transform file,"
                                                         "but failed, because " << err << "\nContinuing."
                                                 << std::endl;
    return -1;
  }
  return 0;
}

MirorrPyramidWrapper::SpacingType
MirorrPyramidWrapper::
GetMaxSpacingForResize(
        ImageType::Pointer image1,
        ImageType::Pointer image2,
        itk::Size<DIMENSION> size) const {
  SpacingType spacing1 = image1->GetSpacing();
  SpacingType spacing2 = image2->GetSpacing();
  itk::Size<DIMENSION> size1 = image1->GetLargestPossibleRegion().GetSize();
  itk::Size<DIMENSION> size2 = image2->GetLargestPossibleRegion().GetSize();

  for (unsigned int ii = 0; ii < DIMENSION; ++ii) {
    spacing1[ii] = spacing1[ii] * size1[ii] / size[ii];
    spacing2[ii] = spacing2[ii] * size2[ii] / size[ii];
    spacing1[ii] = std::max(spacing1[ii], spacing2[ii]);
  }

  return spacing1;
}

//* Create a transform of the correct type given an input string (rigid, quat, affine) */
MirorrPyramidWrapper::InputTransformType::Pointer
MirorrPyramidWrapper::
CreateAppropriateTransform(std::string inputTransformType) {
  //Create an initial transform with the initial translation
  std::string transformFullName;
  //Correct the name of the transformFullName
  //    if( 0 == transformFullName.compare("translation") )
  //      transformFullName = "TranslationTransform_double_";
  //    else
  if (0 == inputTransformType.compare("rigid"))
    transformFullName = "Euler3DTransform_double_";
  else if (0 == inputTransformType.compare("quat"))
    transformFullName = "QuaternionRigidTransform_double_";
  else if (0 == inputTransformType.compare("affine"))
    transformFullName = "AffineTransform_double_";
  else {
    std::cerr << "WARNING: Supplied Transform Type: '"
                 << inputTransformType << "' may not be supported."
                 << " Only translation, rigid and affine are explicitly "
                 << "allowed. Trying anyway."
                 << std::endl;
    transformFullName = inputTransformType;
    transformFullName += "_double_";
  }
  {
    std::stringstream ss;
    ss << DIMENSION << "_" << DIMENSION;
    transformFullName += ss.str();
  }

  //Create a transform factory
  itk::TransformFactoryBase::Pointer transformFactory;
  transformFactory = itk::TransformFactoryBase::New();
  transformFactory->RegisterDefaultTransforms();

  itk::LightObject::Pointer pTransform =
          transformFactory->CreateInstance(transformFullName.c_str());
  InputTransformType::Pointer transform =
          dynamic_cast<InputTransformType *>( pTransform.GetPointer());
  return transform;
}

void
MirorrPyramidWrapper::
ReadInputTransform(
        InputTransformType::Pointer transform,
        std::string transformFileName,
        std::string _fixedName,
        bool _invert_output_transform) {
  if (transformFileName.empty()) {
    //Extract bit of filename after last "/"
    transformFileName =
            _fixedName.substr(_fixedName.find_last_of("/") + 1);
    //Strip off two .dots (if any)
    transformFileName = transformFileName.substr(0,
                                                 transformFileName.find_last_of('.'));
    if (transformFileName.find('.') != std::string::npos)
      transformFileName = transformFileName.substr(0,
                                                   transformFileName.find_last_of('.'));

    transformFileName += "_tfmInitial.tfm";
    transformFileName = std::string("/tmp/") + transformFileName;

    if (verbosity >= 1)
      std::cout << "Setting Initial transform to have translation of "
                   << std::fixed << std::setprecision(3)
                   << PrettyPrint(transform->GetTranslation(), 7) << ". Tfm saved here: " << transformFileName
                   << std::endl;
    writeParametersUsingItkTransformFileWriter(transformFileName, transform);
  }
  else {
    if (readParametersUsingItkTransformFileReader(
            transformFileName, transform, _invert_output_transform)) {
      std::ostringstream name;
      name << "/tmp/example-" << transform->GetTransformTypeAsString() << ".tfm";

      std::cerr << "Could not read input transform file: "
                   << transformFileName << std::endl
                   << "An example transform file has been written here:\n\n"
                   << name.str() << "\n"
                   << "\nTry editing this file and using it as an input." << std::endl;
      writeParametersUsingItkTransformFileWriter(name.str(), transform);

      return;
    }
  }
}

void SplitFileName (const std::string iFileName, std::string &oBaseName, std::string &oFileExt, std::string &oFileExt2)
{
  //Strip folder names
  oBaseName = iFileName;
  if( oBaseName.find('/') != std::string::npos ) {
    oBaseName = oBaseName.substr(oBaseName.find_last_of("/") + 1);
  }

  //Strip extension * 2
  if( oBaseName.find('.') != std::string::npos ) {
    oFileExt2 = oBaseName.substr(oBaseName.find_last_of('.'));
    oBaseName = oBaseName.substr(0, oBaseName.find_last_of('.'));
  }

  if( oBaseName.find('.') != std::string::npos ) {
    oFileExt = oBaseName.substr(oBaseName.find_last_of('.'));
    oBaseName = oBaseName.substr(0, oBaseName.find_last_of('.'));
  }
}

/*
 * Run the registration
 */

void
MirorrPyramidWrapper::
Update()
{
  boost::timer::cpu_timer Mirorr_timer;

  //Create input structure for MirorrPyramidImplement and fill this in.
  ReadAndResampleImages();

  InputTransformType::Pointer transform = CreateAppropriateTransform( transformType );
  if( verbosity >= 1 )
    std::cout << "Transform Class: " << transform->GetNameOfClass() << std::endl;

  //Centre the transforms on the image centre
//  ImageType::PointType centre = movingImage->GetOrigin();
//  for(int ii=0; ii<3; ++ii )
//    centre[ii] += 0.5*movingImage->GetSpacing()[ii] * movingImage->GetLargestPossibleRegion().GetSize(ii);
//  transform->SetCenter( centre );
//  transform->SetIdentity();
  {
    typedef itk::CenteredTransformInitializer< TransformType, ImageType, ImageType > TransformInitializer;
    TransformInitializer::Pointer initialiser =  TransformInitializer::New();
    initialiser->SetTransform( transform );
    initialiser->SetFixedImage( movingImage );
    initialiser->SetMovingImage( fixedImage );
    initialiser->GeometryOn();
    initialiser->InitializeTransform();
  }
  ReadInputTransform( transform, initialTransformName, fixedName,
      invert_output_transform );

  //Display the input image


  //Register the two images
  std::cout<<"Loaded  FIXED: "; __MirorrPyramidWrapper::PrintImage<ImageType>(std::cout,fixedImage); std::cout<<"\n";
  std::cout<<"Loaded MOVING: "; __MirorrPyramidWrapper::PrintImage<ImageType>(std::cout,movingImage); std::cout<<"\n";
  mirorr.SetFixedImage( fixedImage );
  mirorr.SetMovingImage( movingImage );
  mirorr.SetFixedMask( fixedMask );
  mirorr.SetMovingMask( movingMask );

  if( verbosity >= 1 ) {
    std::cout
        << std::fixed << std::setprecision(3)
        << "Transformation parameters:       " << PrettyPrint(transform->GetParameters(), 7) << "\n"
        << "Transformation fixed parameters: " << PrettyPrint(transform->GetFixedParameters(), 7)
        << " (usually this is the rotation center)" << std::endl;
  }

  if( !do_not_register ) {
    mirorr.RunMirorr(dynamic_cast<TransformType *>( transform.GetPointer()));
  }
  else {
    std::cout << "\nINFO: Not running registration as requested." << std::endl;
  }

  //Save the parameters
  writeParametersUsingItkTransformFileWriter( finalTransformName,
                                              transform,
                                              invert_output_transform );

  // Save the shifted image if requested to by the user
  if( !lastTransformedFixedName.empty() )
  {
    ImagePointer resampledFixedImage;
    std::cout<<"Resampling FIXED: "; __MirorrPyramidWrapper::PrintImage<ImageType>(std::cout,fixedImage); std::cout<<"\n";
    if( transformType == "rigid" || transformType == "quat"  )
      resampledFixedImage =
          mirorr.GetReorientedImage( dynamic_cast<TransformType*>( transform.GetPointer() ) );
    else
      resampledFixedImage =
          mirorr.GetResampledImage( dynamic_cast<TransformType*>( transform.GetPointer() ) );

    typedef itk::ImageFileWriter< ImageType > WriterType;

    std::cout<<"Output FIXED: "; __MirorrPyramidWrapper::PrintImage<ImageType>(std::cout,resampledFixedImage); std::cout<<"\n";
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( resampledFixedImage );
    writer->SetFileName( lastTransformedFixedName );
    try {
      writer->Update();
    } catch( itk::ExceptionObject & err) {
      std::cerr << "ExceptionObject caught:" << std::endl;
      std::cerr << err << std::endl;
      std::cout << err.what() << std::endl;
    }

  }

  //===========================================================================

  //===========================================================================

  //And save the shifted image if requested to by the user
  if( !lastTransformedMovingName.empty() )
  {
    ImagePointer resampledMovingImage;
    std::cout<<"Resampling MOVING: "; __MirorrPyramidWrapper::PrintImage<ImageType>(std::cout,movingImage); std::cout<<"\n";

    if( transformType == "rigid" || transformType == "quat"  )
      resampledMovingImage =
          mirorr.GetReorientedImage( dynamic_cast<TransformType*>( transform.GetPointer() ), true );
    else
      resampledMovingImage =
          mirorr.GetResampledImage( dynamic_cast<TransformType*>( transform.GetPointer() ), true );
    std::cout<<"Output MOVING: "; __MirorrPyramidWrapper::PrintImage<ImageType>(std::cout,resampledMovingImage); std::cout<<"\n";

    typedef itk::ImageFileWriter< ImageType > WriterType;

    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( resampledMovingImage );
    writer->SetFileName( lastTransformedMovingName );
    try {
      writer->Update();
    } catch( itk::ExceptionObject & err) {
      std::cerr << "ExceptionObject caught:" << std::endl;
      std::cerr << err << std::endl;
      std::cout << err.what() << std::endl;
    }
  }


  if( verbosity >= 1 )
  {
    std::cout << "\nRegistration and all IO completed in: "
              << std::setprecision(4) << Mirorr_timer.elapsed().wall/1.E9 << "s" << std::endl;
    std::cout << "Output Tfm written here: " << finalTransformName << std::endl;

    if( verbosity >= 1 && !do_not_register )
    {
      // Print resampling instructions
      std::string fixedBase, fixedExt, fixedExt2;
      SplitFileName(fixedName, fixedBase, fixedExt, fixedExt2);

      std::string movingBase, movingExt, movingExt2;
      SplitFileName (movingName, movingBase, movingExt, movingExt2);

      std::cout << "\nTo get resampled FIXED image run:\n" << this->program_name << " --do-not-register";
      if( invert_output_transform ) {std::cout << " --invert-last-tfm";}
      std::cout << " -m " << movingName << " -f " << fixedName << " -l " << finalTransformName
                << " --save-fixed " << fixedBase << "_resampled" << fixedExt << fixedExt2  << std::endl;

      std::cout << "To get resampled MOVING image run:\n" << this->program_name << " --do-not-register";
      if( invert_output_transform ) {std::cout << " --invert-last-tfm";}
      std::cout << " -m " << movingName << " -f " << fixedName << " -l " << finalTransformName
                << " --save-moving " << movingBase << "_resampled" << movingExt << movingExt2  << std::endl;

      std::cout << "Note: it is possible to use --save-fixed and --save-moving at the same time" << std::endl;
    }
  }
}

}

#endif /* ITKMIRORRPYRAMIDWRAPPER_TXX_ */
