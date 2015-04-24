
#include "NPWCudaHistogramCreator.h"

#include "NPWCudaConstants.h"

#include "MutualInformation.h"


#include <iostream>
#include <algorithm>

NPWCudaHistogramCreator::NPWCudaHistogramCreator() :
    geometryData ( 0 ),
    image1Data ( 0 ),
    image2Data ( 0 ),
    renderer ( 0 ),
    nBins ( 0 )
{
    renderer = new NPWCudaOpenGLRenderer();
}

NPWCudaHistogramCreator::~NPWCudaHistogramCreator()
{
    ClearDataStructures();

    if (renderer)
    {
        delete renderer;
        renderer = 0;
    }
}

void NPWCudaHistogramCreator::ClearDataStructures()
{
    if (geometryData)
    {
        delete geometryData;
        geometryData = 0;
    }

    if (image1Data)
    {
        delete image1Data;
        image1Data = 0;
    }

    if (image2Data)
    {
        delete image2Data;
        image2Data = 0;
    }
}

float * NPWCudaHistogramCreator::getHistogram(
                            float * image1, float * image2,
                            unsigned int dimX, unsigned int dimY, unsigned int dimZ,
                            int bins)
{
    imageDims[0] = dimX;
    imageDims[1] = dimY;
    imageDims[2] = dimZ;
    unsigned int imageSize = dimX * dimY * dimZ;

    if ( bins != -1 )
    {
//        nBins = bins;
        setBinSize1( image1, imageSize, nBins );
        setBinSize2( image2, imageSize, nBins );
    }
//    std::cout << "nBins: " << nBins << " binSize1: " << binSize1
//            << " binSize2: " << binSize2 << std::endl;

//    std::cout << "NPWCudaHistogramCreator::getHistogram() attempting to create CUDA images" << std::endl;

    image1Data = CreateImageData( image1, imageSize );
    image2Data = CreateImageData( image2, imageSize );

//    std::cout << "NPWCudaHistogramCreator::getHistogram() attempting to create CUDA geometry data" << std::endl;
    CreateGeometryData();

//    std::cout << "NPWCudaHistogramCreator::getHistogram() attempting to create host histogram" << std::endl;
    int histogramSize = nBins * nBins;
    float * histogram = new float [histogramSize];
    memset(histogram, 0, sizeof(float) * histogramSize);

//    std::cout << "NPWCudaHistogramCreator::getHistogram() created " << nBins << " x "
//              << nBins << " histogram and will now attempt to run kernels" << std::endl;

    bool result = RunHistogramKernels( histogram );
    if (!result)
    {
        delete [] histogram;
        histogram = 0;
    }

    ClearDataStructures();

//    double total = 0.0;
//    for (int i = 0; i < nBins*nBins; ++i)
//    {
//        total += (double)histogram[i];
//    }
//    std::cout << "NPWCudaHistogramCreator::getHistogram() total of histogram : " << total << std::endl;


    return histogram;
}


void NPWCudaHistogramCreator::CreateGeometryData()
{
    unsigned int n_geoUnits = (imageDims[0]-1) * (imageDims[1]-1) * (imageDims[2]-1) * 5;
    geometryData = new NPWCudaDataPointer<float> (n_geoUnits * FLOATS_PER_GEO_UNIT, 0, true, true);

    return;
}

NPWCudaDataPointer<float> * NPWCudaHistogramCreator::CreateImageData( float * data, unsigned int size )
{
//    std::cout << "NPWCudaHistogramCreator::CreateImageData() creating new cuda data pointer of "
//              << size << " floats" << std::endl;
    NPWCudaDataPointer<float> * cudaData =
            new NPWCudaDataPointer<float>( size, data, true, true);

    return cudaData;
}

bool NPWCudaHistogramCreator::RunHistogramKernels( float * histogram )
{
    if ( !geometryData )
    {
        std::cout << "NPWCudaHistogramCreator::RunHistogramKernels() geometry data null pointer!" << std::endl;
        return false;
    }
    if ( !image1Data   )
    {
        std::cout << "NPWCudaHistogramCreator::RunHistogramKernels() image 1 data null pointer!" << std::endl;
        return false;
    }
    if ( !image2Data   )
    {
        std::cout << "NPWCudaHistogramCreator::RunHistogramKernels() image 2 data null pointer!" << std::endl;
        return false;
    }

    // get geometry
    if ( !ExtractGeometryFromImagesCUDA(
            imageDims[0], imageDims[1], imageDims[2],
            geometryData->GetRawCudaPointer(),
            image1Data->GetRawCudaPointer(),
            image2Data->GetRawCudaPointer(),
            image1Min, image2Min,
            binSize1, binSize2) )
    {
        std::cout << "NPWCudaHistogramCreator::RunHistogramKernels() failed to run kernel to extract geometry" << std::endl;
        return false;
    }

    // render geometry
    bool result = renderer->renderGeometryOnResult(
                        histogram, nBins,
                        geometryData->GetRawCudaPointer(), geometryData->GetSize() );
    return result;
}


void NPWCudaHistogramCreator::setBinSize(
        const float * imageData, unsigned int dataSize, float & imageMin, float & binSize )
{
    if ( !imageData )
    {
        std::cout << "NPWCudaHistogramCreator::setBinSize() image data null pointer" << std::endl;
        return;
    }

    if (nBins < 2)
    {
        std::cout << "NPWCudaHistogramCreator::setBinSize() invalid number of bins specified: "
                  << nBins << std::endl;
        return;
    }

    imageMin       = *std::min_element( imageData, imageData + dataSize );

    float imageMax = *std::max_element( imageData, imageData + dataSize );

    binSize        = (imageMax - imageMin) / (float)nBins;

//    std::cout << "set bin size to " << binSize << " for range ["
//              << imageMin << " - " << imageMax << "] and " << nBins << "bins" << std::endl;

}

float NPWCudaHistogramCreator::getMutualInformation(
                      float * image1, float * image2,
                      unsigned int dimX, unsigned int dimY, unsigned int dimZ,
                      int bins)
{
//    std::cout << "NPWCudaHistogramCreator::getMutualInformation() the images are said to be "
//              << dimX << " x " << dimY << " x " << dimZ << std::endl;

    float * histogram = 0;
    histogram = getHistogram( image1, image2, dimX, dimY, dimZ, bins );

    if (!histogram)
    {
        std::cout << "NPWCudaHistogramCreator::getMutualInformation() error rendering histogram" << std::endl;
        return -1.0f;
    }

    float MI = NPW::getMutualInformationFromHistogram( histogram, nBins );

    if (histogram)
    {
        delete [] histogram;
        histogram = 0;
    }

    return MI;
}

void NPWCudaHistogramCreator::setBins( unsigned int bins )
{
    if ((int)bins < renderer->getFrameBufferBorderSize() )
    {
        std::cout << "cannot render on histograms smaller than "
                  << renderer->getFrameBufferBorderSize() << " x "
                  << renderer->getFrameBufferBorderSize() << ", setting bins to "
                  << 2*renderer->getFrameBufferBorderSize()
                  << " instead of " << bins << std::endl;
        bins = 2*renderer->getFrameBufferBorderSize();
    }
    nBins = bins;
}


