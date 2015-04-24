/*
 * NPWCudaOpenGLRenderer.cxx
 *
 *  Created on: 22/07/2010
 *      Author: bro86j
 */

#include "NPWCudaOpenGLRenderer.h"
#include "NPWCudaKernelPrototypes.h"
#include "NPWCudaDataPointer.h"
#include <iostream>
#include <cuda_gl_interop.h>

NPWCudaOpenGLRenderer::NPWCudaOpenGLRenderer() :
    window ( 0 )
{
    window = new NPWCudaWindow(800,600,false);

    cudaGLSetGLDevice(0);
}

NPWCudaOpenGLRenderer::~NPWCudaOpenGLRenderer()
{
    if ( window )
    {
        delete window;
        window = 0;
    }
}


bool NPWCudaOpenGLRenderer::renderGeometryOnResult(
            float * resultData, unsigned int bins,
            float * d_geoData   , unsigned int geoDataSize )
{
    window->initialize( bins, bins );

    float offset = (float)window->getFrameBufferBorderSize();

    NPWCudaDataPointer<float> * pointResultData =
            new NPWCudaDataPointer<float>(bins*bins, 0, true,true);

    if( !RenderGeometryDirectlyFromCUDA(
                d_geoData,
                geoDataSize,
                offset,
                pointResultData->GetRawCudaPointer(),
                bins
                ) )
    {
        window->finalize();
        return false;
    }

    // get the result of all points
    resultData = pointResultData->DownloadData(resultData);

    if ( pointResultData) { delete pointResultData; pointResultData = 0; }

    // add the OpenGL frame buffer containing all rendered triangles to the result
    window->readBackImage( resultData, bins, bins );

    window->finalize();

    return true;
}

