/*=========================================================================
Program: mirorr
Module: mirorr.cxx
Author: David Rivest-Henault and Nicholas Dowson
Created: Mon 11 Feb 2009

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

// File for running mirorr registration using an image pyramid. See
// comments below for more details.
// Revised: 09 Sep 2014 by D Rivest-Henault. Merged new improved GPU Block Matcher
// Revised: 26 May 2014 by D Rivest-Henault. More resampling space options
// Revised: 04 Apr 2014 by D Rivest-Henault. Handling of non-identity cosine matrix in image headers
// Revised: 15 Oct 2013 by D Rivest-Henault. Using intermediate space in classic mode
// Revised: 30 Jan 2013 by D Rivest-Henault. Added help for using the --use-gpu-bm on a multi GPUs system
// Revised: 15 Dec 2012 by D Rivest-Henault. Added --reg-mod and --nthreads options
// Revised: 10 Aug 2010 by N Dowson. Added explicit masking
// Revised: 13 Oct 2009 by N Dowson. Removed Transform as template parameter
// Revised: 12 Oct 2009 by N Dowson. Added equivalent of -py option and
//   ability to save resampled fixed image after registration.
// Revised: 17 Sep 2009 by N Dowson. Switched moving fixed in naming
//   convention to align with ITK, plus a host of little things to improve
//   usability and reduce annoyances.
// Revised: 13 Jul 2009 by N Dowson. Modification to set translation
//   of translation using set offset
// Revised: 12 Jun 2009 by N Dowson. Customise various things. Renamed
// Revised: 25 Mar 2009 by N Dowson. Added various transformations
// Revised: 13 May 2009 by N Dowson. Can read in and output transform files
// Revised: 20 Mar 2009 by N Dowson. Can now choose input files.
// Revised: 19 Mar 2009 by N Dowson. Generalised to 3D
// Revised: 19 Mar 2009 by N Dowson. First version.

#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <time.h>
#include <sys/types.h>

#if defined(WIN32)
  #include <process.h>
#endif

#include <itkMultiThreader.h>
#include <itkImageIOBase.h>

#include "itkMirorrPyramidWrapper.h"

#include "itkMirorrUtilities.h"

// By default, a rigid transformation is used.
//
// If no starting transformation is specified, the algorithm is
// initialised so that the two image centres overlap, and the initial
// transformation is saved to the name of the fixed image appended
// with tfmInitial.tfm. Likewise if no file for storing the final
// transform (after registration is specified, the fixed image is
// appended with tfmFinal.tfm.
//
// The algorithm relies on extracting blocks from the moving image and
// matching them to the most similar one of the nearby blocks of the
// same size (and resolution) in the fixed image. Only those blocks in
// the moving image with the greatest with their matching block in the
// fixed image are retained. Each of the remaining matches gives a
// suggested transformation. The global tranform agreeing with the
// majority of suggested translations is then found. Several
// iterations are performed at each image scale or resolution at which
// the images are being viewed.
//
// If the \emph{help} argument is used, a short help message is
// displayed, but nothing is done. If the \emph{verbose} argument is
// specified, the status of the algorithm is updated in the command
// line window. This of course may be piped to some test output file
// using pipe command \emph{e.g.} milxMirorr args $>$ output.txt.

/*!
 * \page CommandLine
 * \section mirorr
 * mirorr is a command line application to register
 * pairs of 3D images using linear transforms.
 *
 * Type mirorr --help for help.
 */
int main( int argc, char* argv[] )
{
  namespace po = boost::program_options;
  namespace fs = boost::filesystem;
  //Defaults ============================================================
  // Declare the supported options.
  po::options_description all_options(
  "mirorr\n\n"
  "Synopsis: The command line interface to a robust registration "
  "algorithm based on block matching. The algorithm seeks a transform "
  "from a 'moving' to a 'fixed' image that best aligns the two.\n"
  "Usage:\n"
  "  mirorr movingImage fixedImage [final.tfm] [start.tfm] ... \n\n"
  "Other examples:\n"
  "  To apply transform to moving image:\n"
  "    mirorr movingImage fixedImage last.tfm -R --save-reoriented-moving moving_tfmd.nii.gz\n"
  "  To apply inverse transform to fixed image:\n"
  "    mirorr movingImage fixedImage last.tfm -R --save-reoriented-fixed fixed_tfmd.nii.gz\n"
  "  To run fresh registration from a baseline, save moving inplace, and subsequently refine registration:\n"
  "    mirorr moving.nii.gz fixed.nii.gz last.tfm --fresh --reorient-moving --save-reoriented-moving moving.nii.gz\n"
  "    mirorr moving.nii.gz fixed.nii.gz last.tfm refined.tfm --reorient-moving --save-reoriented-moving moving.nii.gz\n"
  "  To reset fixed (or moving) image to identity orientation "
  "(Useful for baselining fixed image if reorient-fixed was used before):\n"
  "    # mirorr movingImage fixedImage last.tfm --reorient-fixed   ## From before\n"
  "    mirorr movingImage fixedImage --fresh --reorient-fixed -R -0 --save-reoriented-fixed baselineFixed.nii.gz\n"
  );

  po::options_description general_options(
      "General options to set input and output images and transforms");
  general_options.add_options()
    ("moving,m", po::value<std::string>(),
        "Specify Moving Image")
    ("fixed,f", po::value<std::string>(),
        "Specify Fixed Image")
    ("last-tfm,l", po::value<std::string>()->default_value( std::string("") ), //
        "File with optimised transformation from moving to fixed image.")
    ("start-tfm,s", po::value<std::string>()->default_value( std::string("") ),
        "File with initial transformation (will not be overwritten unless it is the same as last-tfm)")
    ("tfm-type,t", po::value<std::string>()->default_value("rigid"),
        "Type of transformation: rigid/affine")
    ("reg-mode", po::value<std::string>()->default_value("symmetric"),
            "Registration mode: classic/symmetric")
    ("do-not-register,R",
        "Do not run registration, just use transform file to resample the "
        "fixed and/or moving image.")
    ("save-reoriented-fixed", po::value<std::string>(),
     "Save fixed image (in moving image space) after registration. "
             "If rigid transform, only direction & origin updated, no resampling "
             "occurs.")
    ("save-reoriented-moving", po::value<std::string>(),
     "Save moving image (in fixed image space) after registration. "
             "If rigid transform, only direction & origin updated, no resampling "
             "occurs.")
    ("save-resampled-fixed", po::value<std::string>(),
     "Save fixed image (resampling on moving image's space and grid) "
             "after registration. ")
    ("save-resampled-moving", po::value<std::string>(),
     "Save moving image (resampling on fixed image's space and grid) "
             "after registration. ")
    ("final-interpolator", po::value<std::string>()->default_value("bspline"),
     "Type of interpolator used to save the final resampled images. Valid options: "
     "nn, linear, bspline, sinc")
    ("moving-mask,M", po::value<std::string>()->default_value(""),
        "Specify Moving Image")
    ("fixed-mask,F", po::value<std::string>()->default_value(""),
        "Specify Fixed Image")
    ("switch-images,S",
        "Switch images, so transform is from fixed to moving image. "
        "DEPRECATED - No longer needed with symmetric approach")
    ("invert-last-tfm,I",
        "Invert output transform. (Maybe use it with --invert-start-tfm)")
    ("invert-start-tfm,i",
        "Invert input transform. (Maybe use it with --invert-last-tfm)")
    ("nthreads", po::value<int>()->default_value(0),
        "Number of block matching threads. "
        "Default is 1 thread per available CPU core.")
    ("echo-cmd",
        "Echo the command line. Useful for logging.")
    ("fresh",
       "PREVENT use of output tfm as an input if it exists. By default "
        "if no input is provided and the output exists the output is "
        "used as an input.")
    ("help,h", "Produce help message")
    ("verbose,v", "Request extra verbose output")
    ("quiet,q", "Request quiet output")
  ;

  po::options_description algorithm_options(
      "Algorithm Options\n"
      "Control the registration algorithm. \nA hierarchical approach is "
      "used on an image pyramid. By convention pyramid level 1 indicates "
      "the moving image is resampled to a resolution of >=16^3. Each increase "
      "indicates a doubling of resolution. The maximum level ensures the "
      "moving image is not downsampled at all. The resolutions of the two "
      "images are matched to ensure approximately equal spacing. "
      "Zero gives default.\nThe number of block matching iterations may "
      "also be applied, as can the fraction of image blocks that are "
      "considered (Low image content blocks are ignored).");
  algorithm_options.add_options()
    ("pyr-start,a", po::value<int>()->default_value(1),
        "First pyramid level to be processed. (-ve numbers count back from max.)")
    ("pyr-end,c", po::value<int>()->default_value(1), //
        "Last pyramid level to be processed. (-ve numbers count back from max. "
        "0==max)")
    ("pyr-num,d", po::value<int>()->default_value(1024),
         "Number of levels in the pyramid (default is to create as many levels "
         "as maxDataSize > 32)")
    ("pyr-min-size,e", po::value<int>()->default_value(32),
         "Minimal image dimension at the first level of the pyramid")
    ("resampling-mode", po::value<std::string>()->default_value("middle"), //valid: basic, middle, max-resolution, max-size
          "Inner loop resampling mode: basic, middle, max-resolution, max-size, "
          "fixed, moving")
    ("crop", po::value<std::string>()->default_value(""),
     "Optionally crop the moving and fixed images on both sides by integer no. voxels using "
     "convention: \"movingX,movingY,fixedX,fixedY\"")
    ("iterations,n", po::value<int>()->default_value(0),
        "Number of iterations. Set negative to halve iterations each pyramid level")
    ("portion-kept", po::value<double>()->default_value(0.5),
        "Portion of image, by volume, to consider.")
#ifdef USE_OPENCL
    ("use-gpu-bm",
        "Use the GPU blockmatcher. This is much faster, but being tested."
        "When using this GPU implementation on a multi GPUs CUDA system, the "
        "CUDA_VISIBLE_DEVICES environment variable should be used to restrict the "
        "visibility of the GPUs to the ones you want to use. "
        "Specificaly on \"bragg-l.csiro.au\", it is advisable to place the following "
        "line in your script:"
        "\"export CUDA_VISIBLE_DEVICES=`gpu-which | sed 's/\"//g'`\"") //
#endif
    ("reorient-fixed",
        "Reorient fixed image in the RAI direction first, and set to reset "
        "position with identity directions and zero origin. This enables "
        "in-place image modification with re-use old transforms.")
    ("reorient-moving",
        "Reorient moving image in the RAI direction first, and set to reset "
        "position with identity directions and zero origin. This enables "
        "in-place image modification with re-use old transforms.")
    ("save-fixed", po::value<std::string>(),
     "Save fixed image (in moving image space) after registration. "
             "If rigid transform, only direction & origin updated, no resampling "
             "occurs. DEPRECATED - will be replaced by --save-resampled-fixed")
    ("save-moving", po::value<std::string>(),
     "Save moving image (in fixed image space) after registration. "
             "If rigid transform, only direction & origin updated, no resampling "
             "occurs. DEPRECATED - will be replaced by --save-resampled-moving")
    ("no-centering-init,0",
       "Turn off initialisation of transform to overlap image centres.")

    ("blockmetric", po::value<std::string>()->default_value("nc"),
        "metric used: normalized correlation (nc), sum of squared differences (sd), "
        "correlation ratio (cr), non-parametric window mutual information (mi)"
#ifdef USE_NPW
        ", np windows mutual information (npwmi)"
#endif
     )
#ifdef USE_NPW
     ("npw-bins", po::value<int>()->default_value(32),
        "the number of bins in the histogram used by the NPW MI metric" )
#ifndef USE_OPENCL
     ("npw-scale-factor", po::value<int>()->default_value(2),
        "render on a bigger framebuffer and scale this down to reduce render error "
        "this number specifies how many times bigger the framebuffer is" )
     ("npw-shape-expansion",
        "expands shapes when rendering reducing rasterization error "
        "if it is off (default), the histogram is normalized afterwards" )
#endif
#endif
     ("nhoodwidth", po::value<int>()->default_value(7),
         "Width of block matching neighbourhood (number of blocks tested along one edge of the neighbourhood [=7])")
     ("nhoodgap", po::value<int>()->default_value(3),
         "Gap between block matching neighbourhoods (origin BM space, in pixels [=3])")
     ("blockwidth", po::value<int>()->default_value(4),
         "Size of blocks (pixels [=4])")
     ("blockgap", po::value<int>()->default_value(1),
         "Gap between blocks in neighbourhood (search BM space, in pixels [=1])")
   ;

  all_options.add(general_options).add(algorithm_options);
  //f F h H l L m M n q R s S t v b

  po::positional_options_description positionalDescription;
  positionalDescription.add("moving", 1).add("fixed",1).add("last-tfm",1).add("start-tfm",1);

  po::variables_map variablesMap;
  po::store(po::command_line_parser(argc, argv).
            options(all_options).
            positional(positionalDescription).run(),
            variablesMap);
  po::notify(variablesMap);

  //Opening text ========================================================
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  std::cout << "------------------------------------------------------------------------\n"
            << "- CSIRO Mirorr (c) 2009-2015                                           -\n"
            << "------------------------------------------------------------------------\n"
            << "Multimodal Image Registration using blOck-matching and Robust Regression\n"
            << "By: David Rivest-Henault, Nicholas Dowson, and the AEHRC team\n"
            << "Build time stamp: <" << __DATE__ << "-" << __TIME__ << ">\n"
            << "Current time is:  " << asctime (timeinfo)
            << "------------------------------------------------------------------------\n";

  //Parse options =======================================================
  if (variablesMap.count("echo-cmd"))
  {
    #if defined(WIN32)
      std::cout<<"PID: " << _getpid() << std::endl;
    #else
      std::cout<<"PID: " << getpid() << std::endl;
    #endif
    for( int ii =0; ii<argc; ++ii )
      std::cout << argv[ii] << " " ;
    std::cout<< std::endl;
  }
  else
  {
    std::cout << "(use --echo-cmd to get more process information)" << std::endl;
  }
  std::cout << "------------------------------------------------------------------------\n"
            << std::endl;

  if (variablesMap.count("help"))
  {
  std::cout << all_options << std::endl << std::endl;
  return 1;
  }

  //The names of the moving and fixed files are mandatory
  if( !variablesMap.count("moving") )
  {
    std::cerr << "No moving image supplied. "
    "Type 'mirorr --help' for usage." <<std::endl;
    return -1;
  }
  if( !variablesMap.count("fixed") )
  {
    std::cerr << "No fixed image supplied. "
    "Type 'mirorr --help' for usage."<<std::endl;
    return -1;
  }

  //Check dimensionality of input images
  itk::ImageIOBase::Pointer imageIO =
      itk::ImageIOFactory::CreateImageIO(variablesMap["moving"].as<std::string>().c_str(),
          itk::ImageIOFactory::ReadMode);

  //Check file exists else imageIO falls over
  if( !fs::exists(variablesMap["moving"].as<std::string>()) )
    {
    std::cerr <<"ERROR: " << variablesMap["moving"].as<std::string>() << " not found. Exiting!"
    << std::endl;
    return -1;
    }
  if( !fs::exists(variablesMap["fixed"].as<std::string>()) )
    {
    std::cerr <<"ERROR: " << variablesMap["fixed"].as<std::string>() << " not found. Exiting!"
    << std::endl;
    return -1;
    }

  imageIO->SetFileName(variablesMap["moving"].as<std::string>().c_str());
  imageIO->ReadImageInformation();

  const int nDims = imageIO->GetNumberOfDimensions();

  //Return error if dimensionality is not 2 or 3
  //We have to template on dimension, so other dimensionalities
  //would require excess code
  //Case table for different dimensionalities and transformations
  //Although this is code verbose and results in large binaries - it results
  //in efficient code
  if( nDims != 3 )
    {
    std::cerr << "ERROR: Image dimension is " << nDims
    << "Only 3 dimensional registrations are currently supported."
    << std::endl;
    return -1;
    }
  itk::MirorrPyramidWrapper mirorr;

  //The registration mode need to be set first
  //***THIS CREATE A NEW (Sym)MirorrRegistrationMethod***
  if( variablesMap.count("reg-mode") ) {
    std::string mode = variablesMap["reg-mode"].as<std::string>();
    if( mode == "classic" ) {
      mirorr.GetRegistrationPyramidObject().SetRegistrationTypeClassic();
    } else if( mode == "symmetric" ) {
      mirorr.GetRegistrationPyramidObject().SetRegistrationTypeSymmetrical();
    } else {
      std::cerr << "ERROR: Unknown reg-mode. Available options are "
          "'classic' and 'symmetric' (default)." << std::endl;
      return -1;
    }
  } else {
    mirorr.GetRegistrationPyramidObject().SetRegistrationTypeSymmetrical();
  }
  mirorr.SetProgramName(argv[0]);

  //Get moving and fixed image file names from command line
  if( variablesMap.count("switch-images") )
    {
    mirorr.SetMovingName( variablesMap["fixed"].as<std::string>() );
    mirorr.SetFixedName( variablesMap["moving"].as<std::string>() );
    mirorr.SetMovingMaskName( variablesMap["fixed-mask"].as<std::string>() );
    mirorr.SetFixedMaskName( variablesMap["moving-mask"].as<std::string>() );
    }
  else
    {
    mirorr.SetMovingName( variablesMap["moving"].as<std::string>() );
    mirorr.SetFixedName( variablesMap["fixed"].as<std::string>() );
    mirorr.SetMovingMaskName( variablesMap["moving-mask"].as<std::string>() );
    mirorr.SetFixedMaskName( variablesMap["fixed-mask"].as<std::string>() );
    }
  mirorr.SetDoReorientFixedInARI( variablesMap.count("reorient-fixed") == true );
  mirorr.SetDoReorientMovingInARI( variablesMap.count("reorient-moving") == true );
  if( variablesMap.count("no-centering-init") )
    mirorr.SetDoInitialiseTransformToCentreImages(false);

  if( variablesMap.count("crop") )
  {
    std::vector<unsigned int> crop;
    itk::util::parseInputToList<unsigned int>(
            variablesMap["crop"].as<std::string>(), 4, crop);
    mirorr.SetCropMargins(crop[0], crop[1], crop[2], crop[3]);
  }

    //Should we bother applying registration
  bool do_not_register = variablesMap.count("do-not-register") > 0;
  mirorr.SetDoNotRegister( do_not_register );

  //Read in the name of the file for the final transformation or create one
  //from fixed image name
  std::string output_transform_file_name = variablesMap["last-tfm"].as<std::string>();
  std::string  input_transform_file_name = variablesMap["start-tfm"].as<std::string>();

  if( output_transform_file_name.empty() && ! do_not_register )
  {
    //Two basenames to catch .nii.gz
    output_transform_file_name =
        fs::basename( fs::basename( fs::path(mirorr.GetFixedName()).leaf() ) );
    output_transform_file_name += "_tfmFinal.tfm";
    std::cerr << "WARNING: No output transform supplied. Setting to: "
    << output_transform_file_name << std::endl;
    std::cerr << "INFO:    To avoid saving the last transform, use '--last-tfm ignore'" << std::endl;
  } else if( output_transform_file_name == "ignore" ) {
    output_transform_file_name = "";
  }
  mirorr.SetFinalTransformName( output_transform_file_name );

  //Check if output file exists
  if( input_transform_file_name.empty()
      && fs::exists(fs::path(output_transform_file_name)) )
    {
    std::cout << "INFO: Output transform file already exists ("
    << output_transform_file_name <<") and will be overwritten, ";
    if( variablesMap.count("fresh") )
      std::cout << "and it will NOT be used to initialise registration.\n" << std::endl;
    else
      {
      std::cout << "but it WILL be used to initialise registration.\n" << std::endl;
      input_transform_file_name = output_transform_file_name;
      }
    }

  //Read in name of the file for the initial transformation or assume
  //identity transform
  mirorr.SetInitialTransformName( input_transform_file_name );

  //Specify name of file to save fixed image to after registration
  if(variablesMap.count("save-fixed"))
    mirorr.SetLastTransformedFixedName( variablesMap["save-fixed"].as<std::string>(), true);
  if(variablesMap.count("save-moving"))
    mirorr.SetLastTransformedMovingName( variablesMap["save-moving"].as<std::string>(), true);
  if(variablesMap.count("save-reoriented-fixed"))
    mirorr.SetLastTransformedFixedName( variablesMap["save-reoriented-fixed"].as<std::string>() );
  if(variablesMap.count("save-reoriented-moving"))
    mirorr.SetLastTransformedMovingName( variablesMap["save-reoriented-moving"].as<std::string>() );

  if(variablesMap.count("save-resampled-fixed"))
  {
    if (!mirorr.GetLastTransformedFixedName().empty())
      throw std::runtime_error("Cannot use save-resampled-fixed and save-reoriented-fixed arguments at same time.");
    mirorr.SetLastTransformedFixedName(variablesMap["save-resampled-fixed"].as<std::string>(), true);
  }
  if(variablesMap.count("save-resampled-moving"))
  {
    if (!mirorr.GetLastTransformedMovingName().empty())
      throw std::runtime_error("Cannot use save-resampled-moving and save-reoriented-moving arguments at same time.");
    mirorr.SetLastTransformedMovingName(variablesMap["save-resampled-moving"].as<std::string>(), true);
  }

  if(variablesMap.count("final-interpolator"))
  {
      mirorr.SetFinalInterpolatorName(variablesMap["final-interpolator"].as<std::string>());
  }

  if( variablesMap.count("invert-start-tfm") )
    mirorr.SetInvertInputTransform( true );
  else
    mirorr.SetInvertInputTransform( false );

  if( variablesMap.count("invert-last-tfm") )
    mirorr.SetInvertOutputTransform( true );
  else
    mirorr.SetInvertOutputTransform( false );

  // Set the number threads
  if( variablesMap.count("nthreads") ) {
    int n = variablesMap["nthreads"].as<int>();
    if( mirorr.GetVerbosity() >= 2 ) {std::cout << "--nthread: " << n;}
    if (n==0) {
      if( mirorr.GetVerbosity() >= 2 ) {std::cout << ". Default multithread implementation" << std::endl;}
      mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->UseMultiThreadingOn();
      // Use default number of threads
    } else if (n<2) {
      if( mirorr.GetVerbosity() >= 2 ) {std::cout << ". Singlethreaded implementation" << std::endl;}
      mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->UseMultiThreadingOff();
    } else {
      if( mirorr.GetVerbosity() >= 2 ) {std::cout << ". Multithread implementation with " << n << " threads" << std::endl;}
      mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->UseMultiThreadingOn();
      mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->SetNumberOfBlockMatcherThreads(n);
      itk::MultiThreader::SetGlobalMaximumNumberOfThreads( n );
      
    }
  }

  //Store no. iterations
  mirorr.GetRegistrationPyramidObject().SetMaxIterations( variablesMap["iterations"].as<int>() );

  /*
  //Resampling for non-symmetric registration
  mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->SetUseSymResampling( variablesMap.count("basic-resampling") == 0 );

  //Resampling type for (semi)symmetric registration
  mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->SetUseMaxResampling( variablesMap.count("max-resampling") != 0 );
  */
  if( variablesMap.count("resampling-mode") ) {
    std::string mode = variablesMap["resampling-mode"].as<std::string>();
    if( mode == "basic" ) {
      mirorr.GetRegistrationPyramidObject().SetResamplingType(itk::MirorrPyramidWrapper::MirorrType::EResamplingBasic);
    } else if( mode == "middle" ) {
      mirorr.GetRegistrationPyramidObject().SetResamplingType(itk::MirorrPyramidWrapper::MirorrType::EResamplingMiddle);
    } else if( mode == "max-resolution" ) {
      mirorr.GetRegistrationPyramidObject().SetResamplingType(itk::MirorrPyramidWrapper::MirorrType::EResamplingMaxResolution);
    } else if( mode == "max-size" ) {
      mirorr.GetRegistrationPyramidObject().SetResamplingType(itk::MirorrPyramidWrapper::MirorrType::EResamplingMaxSize);
    } else if( mode == "fixed" ) {
      mirorr.GetRegistrationPyramidObject().SetResamplingType(itk::MirorrPyramidWrapper::MirorrType::EResamplingFixed);
    } else if( mode == "moving" ) {
      mirorr.GetRegistrationPyramidObject().SetResamplingType(itk::MirorrPyramidWrapper::MirorrType::EResamplingMoving);
    } else {
      std::cerr << "ERROR: Unknown resampling-mode. Available options are "
          "'classic' and 'symmetric' (default)." << std::endl;
      return -1;
    }
  } else {
    mirorr.GetRegistrationPyramidObject().SetResamplingType(itk::MirorrPyramidWrapper::MirorrType::EResamplingMiddle);
  }

  //Control pyramid levels that are used
  int aa = variablesMap["pyr-start"].as<int>();
  int cc = variablesMap["pyr-end"].as<int>();
  int maxLevel = variablesMap["pyr-num"].as<int>();
  int minBlockLength = variablesMap["pyr-min-size"].as<int>();
  mirorr.GetRegistrationPyramidObject().GetPyramidScheduleTuner()->SetMinLength(minBlockLength);
  mirorr.GetRegistrationPyramidObject().GetPyramidScheduleTuner()->SetLevelMin( aa );
  mirorr.GetRegistrationPyramidObject().GetPyramidScheduleTuner()->SetLevelMax( cc );

  mirorr.GetRegistrationPyramidObject().GetPyramidScheduleTuner()->SetMaxLevelNumber(maxLevel);

  //Read bottom level. This should be at least 0
  mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->SetPortionMatchesKept(
      variablesMap["portion-kept"].as<double>() );
#ifdef USE_OPENCL
  if(variablesMap.count("use-gpu-bm"))
    mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->SetUseGpuOn(true);
  else
    mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->SetUseGpuOn(false);
#endif

  //What level of verbosity is required
  mirorr.SetVerbosity(1);
  if( variablesMap.count("quiet") )
    mirorr.SetVerbosity(0);
  if( variablesMap.count("verbose") )
    mirorr.SetVerbosity(2);

  //Read in the default transform type
  mirorr.SetTransformType( variablesMap["tfm-type"].as<std::string>() );

  //Read in block metric type
  {
    std::map<std::string, std::string> metricMap;
    metricMap["nc"] = "normalized_correlation";
    metricMap["sd"] = "sum_of_squared_differences";
    metricMap["cr"] = "correlation_ratio";
    metricMap["mi"] = "mutual_information";
#ifdef USE_NPW
    metricMap["npwmi"] = "npw_mutual_information";
#endif
    std::string blockMetricInput = variablesMap["blockmetric"].as<std::string>();
    if ( metricMap.find(blockMetricInput) == metricMap.end() )
      {
      std::cerr << "ERROR: " << blockMetricInput << " is not a valid blockmetric type. Exiting!" << std::endl;
      std::cerr << "possible options are:" << std::endl; \
      for( std::map<std::string,std::string>::iterator it = metricMap.begin(); it != metricMap.end(); ++it )
        {
        std::cerr << "  " << it->first << ": " << it->second << std::endl;
        }
      return -1;
      }
    std::string blockMetricType = metricMap[ blockMetricInput ];
    mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->SetBlockMetricType( blockMetricType );
  }

  //Fill in deep algorithm parameters
  //Optionally fill in "deep" algorithm parameters
  {
    int blockWidth = variablesMap["blockwidth"].as<int>();
    int blockGap   = variablesMap["blockgap"].as<int>();
    int nhoodWidth = variablesMap["nhoodwidth"].as<int>();
    int nhoodGap   = variablesMap["nhoodgap"].as<int>();

    if( blockWidth>0 ) mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->SetBlockWidth( blockWidth );
    if( blockGap>0   ) mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->SetBlockGap(   blockGap );
    if( nhoodWidth>0 ) mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->SetNhoodWidth( nhoodWidth );
    if( nhoodGap>0   ) mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->SetNhoodGap(   nhoodGap );
#ifdef USE_NPW
    //Fill in npw algorithm parameters
    int NPWbins           = variablesMap["npw-bins"].as<int>();
    mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->SetNPWbins( NPWbins );
#ifndef USE_OPENCL
    int NPWscaleFactor    = variablesMap["npw-scale-factor"].as<int>();
    int NPWshapeExpansion = ( variablesMap.count("npw-shape-expansion") ? true : false );
    mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->SetNPWscaleFactor( NPWscaleFactor );
    mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->SetNPWshapeExpansion( NPWshapeExpansion );
#endif
#endif
  }

  //Display
  if( mirorr.GetVerbosity() >= 2 )
    {
    std::string switched = "";
    if(variablesMap.count("switch-images")) switched = " (switched)";
    std::cout << "Moving Image: " << mirorr.GetMovingName() << switched << std::endl;
    std::cout << " Fixed Image: " << mirorr.GetFixedName() << switched << std::endl;
    std::cout << "   Data Type: " << imageIO->GetComponentType() << std::endl;
    std::cout << "  Dimensions: " << nDims  << std::endl;
    std::cout << "   Transform: " << mirorr.GetTransformType() << std::endl;
    std::cout << "  Output Tfm: " << mirorr.GetFinalTransformName() << std::endl;
    std::cout << "Block Metric: " << mirorr.GetRegistrationPyramidObject().GetRegistrationObject()->GetBlockMetricType() << std::endl;
    }
  mirorr.Update();


  time_t rawtimeEnd; /* Seconds since the Epoch.  */
  time ( &rawtimeEnd );
  time_t elapsed = rawtimeEnd - rawtime;

  std::cout << "\nmirorr ran successfully" << std::endl
            << "Process started on:   " << asctime ( localtime ( &rawtime    ) );
  std::cout << "Process completed on: " << asctime ( localtime ( &rawtimeEnd ) )
            << "Elapsed time is ";
  if (elapsed/3600 > 0) { std::cout << elapsed/3600 << "h"; }
  if (elapsed/60   > 0) { std::cout << (elapsed/60)%60 << "m"; }
  std::cout << elapsed%60 << "s" << std::endl;

  return 0;
}
