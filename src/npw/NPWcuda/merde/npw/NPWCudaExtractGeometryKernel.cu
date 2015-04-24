/*
 * NPWCudaExtractGeometryKernel.cu
 *
 *  Created on: 19/07/2010
 *      Author: bro86j
 */

#include <stdio.h>

#include "NPWCudaGeometryFunctions.h"
#include "NPWCudaShapeDetermination.h"

#include "NPWCudaDataPointer.h"


/*
 * extracts geometry data from two images
 */
__global__ void ExtractGeometryFromImagesKernel(
        float *      d_geometryData,
        float *      d_image1Data,
        float *      d_image2Data,
        unsigned int imageDimX,
        unsigned int imageDimY,
        unsigned int imageDimZ,
        float        innerWeight,
        float        outerWeight,
        float        image1Min,
        float        image2Min,
        float        invBinSize1,
        float        invBinSize2
        )
{
    // these are neighborhood indices,
    // corresponding to the first voxel of the neighborhood
    unsigned int x = blockIdx.x;
    unsigned int y = blockIdx.y;
    unsigned int z = threadIdx.x;

    // copy neighborhood to local memory, to speed things up
    float2 neighborhood_cache[8];
    for ( int i = 0; i < 8; ++i )
    {
        neighborhood_cache[i] = make_float2(0.0f,0.0f);
    }

    // check if this thread neighborhood is within
    if ( x < imageDimX-1 && y < imageDimY-1 && z < imageDimZ-1 )
    {
        // get 5 tets
        neighborhood_cache[0].x = d_image1Data[ getArrayIndex3D(x  ,y  ,z  ,imageDimX,imageDimY) ];
        neighborhood_cache[0].y = d_image2Data[ getArrayIndex3D(x  ,y  ,z  ,imageDimX,imageDimY) ];

        neighborhood_cache[1].x = d_image1Data[ getArrayIndex3D(x+1,y  ,z  ,imageDimX,imageDimY) ];
        neighborhood_cache[1].y = d_image2Data[ getArrayIndex3D(x+1,y  ,z  ,imageDimX,imageDimY) ];

        neighborhood_cache[2].x = d_image1Data[ getArrayIndex3D(x  ,y+1,z  ,imageDimX,imageDimY) ];
        neighborhood_cache[2].y = d_image2Data[ getArrayIndex3D(x  ,y+1,z  ,imageDimX,imageDimY) ];

        neighborhood_cache[3].x = d_image1Data[ getArrayIndex3D(x+1,y+1,z  ,imageDimX,imageDimY) ];
        neighborhood_cache[3].y = d_image2Data[ getArrayIndex3D(x+1,y+1,z  ,imageDimX,imageDimY) ];

        neighborhood_cache[4].x = d_image1Data[ getArrayIndex3D(x  ,y  ,z+1,imageDimX,imageDimY) ];
        neighborhood_cache[4].y = d_image2Data[ getArrayIndex3D(x  ,y  ,z+1,imageDimX,imageDimY) ];

        neighborhood_cache[5].x = d_image1Data[ getArrayIndex3D(x+1,y  ,z+1,imageDimX,imageDimY) ];
        neighborhood_cache[5].y = d_image2Data[ getArrayIndex3D(x+1,y  ,z+1,imageDimX,imageDimY) ];

        neighborhood_cache[6].x = d_image1Data[ getArrayIndex3D(x  ,y+1,z+1,imageDimX,imageDimY) ];
        neighborhood_cache[6].y = d_image2Data[ getArrayIndex3D(x  ,y+1,z+1,imageDimX,imageDimY) ];

        neighborhood_cache[7].x = d_image1Data[ getArrayIndex3D(x+1,y+1,z+1,imageDimX,imageDimY) ];
        neighborhood_cache[7].y = d_image2Data[ getArrayIndex3D(x+1,y+1,z+1,imageDimX,imageDimY) ];

        // transform image intensities to probability space
        for ( unsigned int i = 0; i < 8; ++i )
        {
            neighborhood_cache[i].x = ( neighborhood_cache[i].x - image1Min ) * invBinSize1;
            neighborhood_cache[i].y = ( neighborhood_cache[i].y - image2Min ) * invBinSize2;
        }

        // in the geometry data, get the n-hood index ...
        unsigned int nHoodIdx = getArrayIndex3D( x,y,z, imageDimX-1, imageDimY-1 );
        // ... and where we are actually going to write
        unsigned int geo_unit_offset = 5 * FLOATS_PER_GEO_UNIT * nHoodIdx;

        // extract geometry
        // alpha, beta, gamma, epsilon
        getGeometryFromVertices(
                neighborhood_cache[0], neighborhood_cache[1],
                neighborhood_cache[2], neighborhood_cache[4],
                d_geometryData, geo_unit_offset,
                outerWeight );
        geo_unit_offset += FLOATS_PER_GEO_UNIT;

        // beta, gamma, delta, theta
        getGeometryFromVertices(
                neighborhood_cache[1], neighborhood_cache[2],
                neighborhood_cache[3], neighborhood_cache[7],
                d_geometryData, geo_unit_offset,
                outerWeight );
        geo_unit_offset += FLOATS_PER_GEO_UNIT;

        // beta, epsilon, dzeta, theta
        getGeometryFromVertices(
                neighborhood_cache[1], neighborhood_cache[4],
                neighborhood_cache[5], neighborhood_cache[7],
                d_geometryData, geo_unit_offset,
                outerWeight );
        geo_unit_offset += FLOATS_PER_GEO_UNIT;

        // gamma, epsilon, eta, theta
        getGeometryFromVertices(
                neighborhood_cache[2], neighborhood_cache[4],
                neighborhood_cache[6], neighborhood_cache[7],
                d_geometryData, geo_unit_offset,
                outerWeight );
        geo_unit_offset += FLOATS_PER_GEO_UNIT;

        // beta, gamma, epsilon, dzeta
        getGeometryFromVertices(
                neighborhood_cache[1], neighborhood_cache[2],
                neighborhood_cache[4], neighborhood_cache[5],
                d_geometryData, geo_unit_offset,
                innerWeight );
    }

}

/*
 * if we don't care about the histogram, just the mutual information,
 * the joint histogram doesn't need to leave the gpu.
 * in that case, pass the resultData device pointer to render points on,
 * resulting in a faster execution.
 * Otherwise, pass a null-pointer (default)
 */
bool ExtractGeometryFromImagesCUDA(
        unsigned int xDim,
        unsigned int yDim,
        unsigned int zDim,
        float * d_geometryData,
        float * d_image1Data,
        float * d_image2Data,
        float   image1Min,
        float   image2Min,
        float   binSize1,
        float   binSize2
        )
{
    cudaError err;

    // every block is 1D, with as many elements as neighborhoods
    // in the z-dimension (of the image)
    dim3 dimBlock( zDim-1, 1     , 1 );
    // the grid is 2D, one block for every x-y neighborhood
    dim3 dimGrid ( xDim-1, yDim-1, 1 );

    //Precomputes the weightings of the inner tetrahedron and four outer ones and
    //the no. nhoods to do this.
    unsigned int  n_neighborhoods = (xDim-1) * (yDim-1) * (zDim-1);
    float         innerWeight     = 1.0f / ( ((float)n_neighborhoods) * 3.0f );
    float         outerWeight     = innerWeight * 0.5f;

    float invBinSize1 = 1.0f / binSize1;
    float invBinSize2 = 1.0f / binSize2;

    ExtractGeometryFromImagesKernel<<<dimGrid, dimBlock>>>(
            d_geometryData,
            d_image1Data,
            d_image2Data,
            xDim, yDim, zDim,
            innerWeight,
            outerWeight,
            image1Min,
            image2Min,
            invBinSize1,
            invBinSize2
            );

    //Ensures we have finished writing into d_geometryData (we call this after every kernel)
    cudaThreadSynchronize();

    if (((err = cudaGetLastError())) != cudaSuccess)
    {
        printf("Error Executing ExtractGeometryFromImagesKernel: %s\n", cudaGetErrorString(err));
        return false;
    }

    return true;
}


