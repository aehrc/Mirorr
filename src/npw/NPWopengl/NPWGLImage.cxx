/*=========================================================================
  Author:   Daan Broekhuizen
  Date:     May 2010
  Language: C++
=========================================================================*/

#include "NPWGLImage.h"

#ifdef USE_DEVIL
#include <IL/il.h>
#include <IL/ilu.h>
#endif

#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

NPWGLImage::NPWGLImage(const char * fileName, bool verbose)
{
    verbosity = verbose;

	imageDimenions = new int[3];
	imageDimenions[0] = 0;
	imageDimenions[1] = 0;
	imageDimenions[2] = 0;

	voxelCount = 0;
	data = 0;

	if (!LoadData(fileName))
	{
		std::cout << "NPWGLImage::NPWGLImage() Failed to load image from file: " << fileName << std::endl;

		if (data)
		{
			delete [] data;
			data = 0;
		}

		voxelCount = 0;
	}
}

NPWGLImage::NPWGLImage(int xDim, int yDim, int zDim, bool verbose)
{
    verbosity = verbose;
    if ( xDim < 1 || yDim < 1 || zDim < 1 )
    {
        std::cout << "NPWGLImage::NPWGLImage() *** ERROR *** : cannot create image with dimensions "
                  << xDim << "x" << yDim << "x" << zDim << std::endl;
        return;
    }

	imageDimenions = new int[3];
	imageDimenions[0] = xDim;
	imageDimenions[1] = yDim;
	imageDimenions[2] = zDim;

	voxelCount = xDim * yDim * zDim;

	data = new float[voxelCount];

	for (int i = 0; i < voxelCount; i++)
	{
		data[i] = 0.0f;
	}
}

NPWGLImage::NPWGLImage(float * imageData, int xDim, int yDim, int zDim, bool verbose)
{
    verbosity = verbose;
    imageDimenions = new int[3];
    imageDimenions[0] = xDim;
    imageDimenions[1] = yDim;
    imageDimenions[2] = zDim;

    voxelCount = xDim * yDim * zDim;

    data = new float[voxelCount];

    for ( int z = 0; z < zDim; ++z )
    {
        for ( int y = 0; y < yDim; ++y )
        {
            for ( int x = 0; x < xDim; ++x )
            {
                (*this)(x,y,z) = imageData[ x + xDim * ( y + z * yDim ) ];
            }
        }
    }
}

NPWGLImage::~NPWGLImage()
{
	if (data)
	{
		delete [] data;
		data = 0;
	}

	if (imageDimenions)
	{
		delete [] imageDimenions;
		imageDimenions = 0;
	}
}

/**
 * Attempts to load an image from file where the format is:
 * <int xDims> <int yDims> <int zDims>
 * <float intensity>
 * <float intensity>
 * <float intensity>
 * ...
 *
 * If this fails, it cleans up the data structures, setting the
 * voxel count to 0 and making sure data == 0
 */
bool NPWGLImage::LoadData(const char * fileName)
{
	if (!fileName)
	{
		std::cout << "NPWGLImage::LoadData() null file name" << std::endl;
		return false;
	}

//	if (!milxSimFileUtils::VerifyFileExists(fileName))
//	{
//		std::cout << "file not found: " << fileName << std::endl;
//		return false;
//	}

	if( verbosity )
	{
	    std::cout << "NPWGLImage::LoadData() loading image data from file "
	            << fileName << std::endl;
	}

	std::ifstream fileStream(fileName);

	std::string buffer;
    std::string temp;

	if ( !fileStream.good() || !fileStream.rdbuf()->is_open() )
	{
		std::cout << "NPWGLImage::LoadData() failed to open file: '"
		          << fileName <<  "'" << std::endl;
		return false;
	}

	/**
	 * First line contains dimensions
	 */
	if (!fileStream.eof())
	{
		getline(fileStream, buffer);
		std::istringstream line(buffer);

		line >> imageDimenions[0] >> imageDimenions[1] >> imageDimenions[2];

		voxelCount = imageDimenions[0] * imageDimenions[1] * imageDimenions[2];

		data = new float[voxelCount];

		if (!data)
		{
			std::cout << "NPWGLImage::LoadData() Failed to allocate memory for voxel count: " << voxelCount << std::endl;
			return false;
		}

		if( verbosity )
		{
		    std::cout << "NPWGLImage::LoadData() found image dimensions: [" <<
		            imageDimenions[0] << ", " << imageDimenions[1] << ", " << imageDimenions[2] <<
		            "], there are " << voxelCount << " voxels" << std::endl;
		}
	}

	int count = 0;
	float intensity = 0.0f;

	if( verbosity )
    {
        std::cout << "NPWGLImage::LoadData() loading image file" << std::endl;
    }

	while( !fileStream.eof() && count < voxelCount)
	{
		getline(fileStream, buffer);
		std::istringstream line(buffer);

		/*
		 * break we got an empty line
		 */
		if ( buffer.empty() )
		{
		    break;
		}

		line >> intensity;
		SetValueAtIndex(count, intensity);

		++count;
	}

	if (count != voxelCount )
	{
		std::cout << "NPWGLImage::LoadData() end of file was reached before all "
		          << voxelCount << " expected voxels were read in, only " << count
		          << " voxels were found in the file" << std::endl;

	}
	else if ( !fileStream.eof() )
	{
	    std::cout << "NPWGLImage::LoadData() all " << count
	              << " voxels were read before the end of the file was reached" << std::endl;
	}
	else
	{
	    if( verbosity )
	    {
	        std::cout << "NPWGLImage::LoadData() successfully loaded image" << std::endl;
	    }
    }

	fileStream.close();

	return true;
}

bool NPWGLImage::SaveRawDataToText( const char * fileName ) const
{
    /*
     * open stream
     */
    std::ofstream outputFile;

    outputFile.open( fileName );

    if ( !outputFile.is_open() )
    {
        std::cout << "NPWGLImage::SaveRawDataToText() *** ERROR *** : could not open file \""
                  << fileName << "\"" << std::endl;
        return false;
    }

    /*
     * output dimensions
     */
    outputFile << imageDimenions[0] << " "
               << imageDimenions[1] << " "
               << imageDimenions[2] << std::endl;

    /*
     * output data
     */
    for ( int i = 0; i < voxelCount; ++i )
    {
        outputFile << data[i] << std::endl;
    }

    /*
     * close stream
     */
    outputFile.close();

    if ( outputFile.is_open() )
    {
        std::cout << "NPWGLImage::SaveRawDataToText()  *** WARNING *** : could not close file \""
                  << fileName << "\"" << std::endl;
    }
    else
    {
        if( verbosity )
        {
            std::cout << "NPWGLImage::SaveRawDataToText() successfully saved file "
                    << fileName << std::endl;
        }
    }

    return true;

}

#ifdef USE_DEVIL
bool NPWGLImage::SaveToFile( const char * fileName,
        NormalizationMode normalizationMode, float imageMinimum, float imageMaximum ) const
{
    int pixelCount = imageDimenions[0] * imageDimenions[1];
    if ( imageDimenions[2] > 1 )
    {
        std::cout << "NPWGLImage::SaveToFile() WARNING: saving only the first of "
                  << imageDimenions[2] << " slices" << std::endl;
    }
    else if ( voxelCount != pixelCount )
    {
        std::cout << "NPWGLImage::SaveToFile() *** ERROR *** voxel count (" << voxelCount
                  << ") and pixel count (" << pixelCount << ") don't match" << std::endl;
    }

    /**
     * create a copy of the data
     */
    float * dataCopy = new float[pixelCount];
    std::copy( data, data + pixelCount, dataCopy );

    /**
     * normalize it if requested
     */
    switch ( normalizationMode )
    {
        case NO_NORMALIZATION:
        {
            break;
        }

        case LINEAR_SCALE:
        {
            if( verbosity )
            {
                std::cout << "NPWGLImage::SaveToFile() normalizing using a linear scale: [ "
                          << imageMinimum << " - " << imageMaximum << " ] "<< std::endl;
            }
            normaliseData( dataCopy, pixelCount, imageMinimum, imageMaximum, false );
            break;
        }

        case LOG_SCALE:
        {
            if( verbosity )
            {
                std::cout << "NPWGLImage::SaveToFile() normalizing using a log scale: [ "
                          << imageMinimum << " - " << imageMaximum << " ] "<< std::endl;
            }
            normaliseData( dataCopy, pixelCount, imageMinimum, imageMaximum, true );
            break;
        }

        default:
        {
            std::cout << "NPWGLImage::SaveToFile() invalid normalization mode, not normalizing"
                      << std::endl;
            break;
        }
    }

    /**
     * create a ILubyte image
     */
    if( verbosity )
    {
        std::cout << "NPWGLImage::SaveToFile() casting image to unsigned byte..." << std::endl;
    }
    ILuint image;
    ilGenImages(1, &image);
    ilBindImage(image);

    ilTexImage( (ILuint)imageDimenions[0], (ILuint) imageDimenions[1], (ILuint)1,
                (ILubyte)1, IL_LUMINANCE, IL_UNSIGNED_BYTE, 0 );

    ILubyte * imageData = ilGetData();

    for (int i = 0; i < pixelCount; ++i)
    {
        imageData[i] = (ILubyte) dataCopy[i];
    }


    /**
     * save it
     */
    if( verbosity )
    {
        std::cout << "NPWGLImage::SaveToFile() trying to save " << imageDimenions[0] << " x "
                  << imageDimenions[1] << " image \"" << fileName << "\"... " << std::flush;
    }
    ilEnable(IL_FILE_OVERWRITE);
    ilHint(IL_COMPRESSION_HINT, IL_NO_COMPRESSION);
    ilHint(IL_MEM_SPEED_HINT, IL_FASTEST);

    ilSaveImage(fileName);

    ilDeleteImages(1, &image);

    bool errored = false;

    ILenum code = ilGetError();

    if (code != IL_NO_ERROR)
    {
        errored = true;
        std::cout << "\nNPWGLImage::SaveToFile() error performing save operation, cause: "
                  << iluErrorString(code) << std::endl;
    }
    else if ( verbosity )
    {
        std::cout << "DONE!" << std::endl;
    }

    delete [] dataCopy;

    return !errored;
}
#endif

#ifdef USE_DEVIL
void NPWGLImage::normaliseData(
        float * data, const int dataSize,
        const float minReq, const float maxReq, const bool useLogScale ) const
{
    // sort the data
    float * normalisedData = new float[dataSize];
    std::copy( data, data + dataSize, normalisedData );
    std::sort( normalisedData, normalisedData + dataSize, std::greater<float>() );

    // location of first zero in the descending array
    float * firstZero = std::find( normalisedData, normalisedData + dataSize, 0.0 );

    // error checking
    if ( firstZero - normalisedData == 0 )
    {
        std::cout << "NPWGLImage::normaliseData() cannot normalise data because all values are 0.0" << std::endl;
        delete [] normalisedData;
        normalisedData = 0;
        return;
    }
    else if ( firstZero - normalisedData >= dataSize )
    {
        std::cout << "NPWGLImage::normaliseData() cannot normalise data because no 0.0 could be found in data" << std::endl;
        delete [] normalisedData;
        normalisedData = 0;
        return;
    }

    // minval is the first number greater than zero
    float minVal = *(firstZero - 1);
    float maxVal = *normalisedData;


    if ( useLogScale )
    {
        // scale data between minReq and maxReq,
        float alpha = 1e10; // some high value to ensure steep curve in low value range
        float beta  = 1.0f - minVal;
        float gamma = ( maxReq - minReq ) / (float)std::log( alpha * maxVal + beta );
        float delta = minReq;

        for ( int i = 0; i < dataSize; ++i )
        {
            if ( data[i] == 0.0f )
            {
                normalisedData[i] = 0.0f;
            }
            else
            {
                normalisedData[i] = gamma * std::log( data[i] * alpha + beta ) + delta;
            }
        }
    }
    else
    {
        float alpha = ( maxReq - minReq ) / ( maxVal - minVal );
        float beta  = ( maxVal * minReq - maxReq * minVal ) / (maxVal - minVal);

        for ( int i = 0; i < dataSize; ++i )
        {
            if ( data[i] == 0.0f )
            {
                normalisedData[i] = 0.0f;
            }
            else
            {
                normalisedData[i] = data[i] * alpha + beta;
            }
        }
    }

    std::copy( normalisedData, normalisedData + dataSize, data );
    delete [] normalisedData;
}
#endif
