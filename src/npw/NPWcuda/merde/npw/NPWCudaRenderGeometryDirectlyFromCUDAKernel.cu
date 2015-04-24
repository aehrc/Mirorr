/*
 * NPWCudaRenderGeometryDirectlyFromCUDA.cu
 *
 *  Created on: 28/07/2010
 *      Author: bro86j
 */

#include "NPWCudaConstants.h"
#include "NPWCudaGeometryFunctions.h"
#include "NPWCudaDataPointer.h"
#include "NPWCudaVertexBufferPointer.h"
#include <cudpp.h>
#include <stdio.h>


__global__ void GetTriangleCountKernel(
        int *   d_triangleCount,
        int *   d_vertexCount,
        float * d_geoData,
        int     geoUnitCount)
{
    int geoUnit_idx = threadIdx.x + blockIdx.x * blockDim.x;

    if ( geoUnit_idx < geoUnitCount )
    {

        int nTriangles = (int)d_geoData[ geoUnit_idx * FLOATS_PER_GEO_UNIT ];

        d_triangleCount[geoUnit_idx] = nTriangles;

        int nVertices;

        switch ( nTriangles )
        {
            case 0:
            {
                nVertices = 0;
                break;
            }
            case 1:
            {
                nVertices = 3;
                break;
            }
            case 2:
            {
                nVertices = 4;
                break;
            }
            case 3:
            {
                nVertices = 4;
                break;
            }
            case 4:
            {
                nVertices = 5;
                break;
            }

            default:
            {
                nVertices = 0;
                break;
            }
        }

        d_vertexCount[geoUnit_idx] = nVertices;
    }

}

__global__ void FillBuffersKernel(
        float * d_VertexBuffer,
        int *   d_IndexBuffer,
        int *   d_TriangleOffsets,
        int *   d_VertexOffsets,
        float * d_GeoData,
        int     geoUnitCount,
        float   renderOffset,
        float * d_pointResultData,
        int     resultDataDim
        )
{
    int geoUnit_idx = threadIdx.x + blockIdx.x * blockDim.x;

    if ( geoUnit_idx < geoUnitCount )
    {
        int geo_idx               = FLOATS_PER_GEO_UNIT * geoUnit_idx;
            // where in d_geoData the data for this unit begins
        int vertexBuffer_idx      = d_VertexOffsets[geoUnit_idx];
            // which vertex is the first vertex we are going to write
        int vertexfloatBuffer_idx = vertexBuffer_idx * 3;
            // but we're writing 3 floats per vertex, so multiply by 3
        int indexBuffer_idx       = d_TriangleOffsets[geoUnit_idx] * 3;
            // which triangle index is the first index we are going to write
            //   (multiplied by 3, cause there are 3 vertices per triangle)

        switch ( (int) d_GeoData[geo_idx] )
        {
            case 0:
            {
                int x = (int)d_GeoData[geo_idx+2];
                int y = (int)d_GeoData[geo_idx+3];

                int result_idx = getArrayIndex2D( x,y, resultDataDim,resultDataDim );
                atomicAdd( d_pointResultData + result_idx, d_GeoData[geo_idx+1] );

                break;
            }

            case 1:
            {
                d_VertexBuffer[vertexfloatBuffer_idx    ] = d_GeoData[geo_idx+2]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 1] = d_GeoData[geo_idx+3]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 2] = d_GeoData[geo_idx+1];

                d_VertexBuffer[vertexfloatBuffer_idx + 3] = d_GeoData[geo_idx+4]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 4] = d_GeoData[geo_idx+5]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 5] = 0.0f;

                d_VertexBuffer[vertexfloatBuffer_idx + 6] = d_GeoData[geo_idx+6]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 7] = d_GeoData[geo_idx+7]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 8] = 0.0f;

                d_IndexBuffer[indexBuffer_idx    ] = vertexBuffer_idx;
                d_IndexBuffer[indexBuffer_idx + 1] = vertexBuffer_idx+1;
                d_IndexBuffer[indexBuffer_idx + 2] = vertexBuffer_idx+2;

                break;
            }
            case 2:
            {
                d_VertexBuffer[vertexfloatBuffer_idx    ] = d_GeoData[geo_idx+2]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 1] = d_GeoData[geo_idx+3]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 2] = d_GeoData[geo_idx+1];

                d_VertexBuffer[vertexfloatBuffer_idx + 3] = d_GeoData[geo_idx+4]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 4] = d_GeoData[geo_idx+5]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 5] = 0.0f;

                d_VertexBuffer[vertexfloatBuffer_idx + 6] = d_GeoData[geo_idx+6]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 7] = d_GeoData[geo_idx+7]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 8] = 0.0f;

                d_VertexBuffer[vertexfloatBuffer_idx + 9] = d_GeoData[geo_idx+8]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx +10] = d_GeoData[geo_idx+9]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx +11] = 0.0f;

                d_IndexBuffer[indexBuffer_idx    ] = vertexBuffer_idx;
                d_IndexBuffer[indexBuffer_idx + 1] = vertexBuffer_idx+1;
                d_IndexBuffer[indexBuffer_idx + 2] = vertexBuffer_idx+2;

                d_IndexBuffer[indexBuffer_idx + 3] = vertexBuffer_idx;
                d_IndexBuffer[indexBuffer_idx + 4] = vertexBuffer_idx+1;
                d_IndexBuffer[indexBuffer_idx + 5] = vertexBuffer_idx+3;

                break;
            }
            case 3:
            {
                d_VertexBuffer[vertexfloatBuffer_idx    ] = d_GeoData[geo_idx+2]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 1] = d_GeoData[geo_idx+3]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 2] = d_GeoData[geo_idx+1];

                d_VertexBuffer[vertexfloatBuffer_idx + 3] = d_GeoData[geo_idx+4]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 4] = d_GeoData[geo_idx+5]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 5] = 0.0f;

                d_VertexBuffer[vertexfloatBuffer_idx + 6] = d_GeoData[geo_idx+6]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 7] = d_GeoData[geo_idx+7]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 8] = 0.0f;

                d_VertexBuffer[vertexfloatBuffer_idx + 9] = d_GeoData[geo_idx+8]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx +10] = d_GeoData[geo_idx+9]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx +11] = 0.0f;

                d_IndexBuffer[indexBuffer_idx    ] = vertexBuffer_idx;
                d_IndexBuffer[indexBuffer_idx + 1] = vertexBuffer_idx+1;
                d_IndexBuffer[indexBuffer_idx + 2] = vertexBuffer_idx+2;

                d_IndexBuffer[indexBuffer_idx + 3] = vertexBuffer_idx;
                d_IndexBuffer[indexBuffer_idx + 4] = vertexBuffer_idx+2;
                d_IndexBuffer[indexBuffer_idx + 5] = vertexBuffer_idx+3;

                d_IndexBuffer[indexBuffer_idx + 6] = vertexBuffer_idx;
                d_IndexBuffer[indexBuffer_idx + 7] = vertexBuffer_idx+3;
                d_IndexBuffer[indexBuffer_idx + 8] = vertexBuffer_idx+1;

                break;
            }
            case 4:
            {
                d_VertexBuffer[vertexfloatBuffer_idx     ] = d_GeoData[geo_idx + 2]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx +  1] = d_GeoData[geo_idx + 3]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx +  2] = d_GeoData[geo_idx + 1];

                d_VertexBuffer[vertexfloatBuffer_idx +  3] = d_GeoData[geo_idx + 4]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx +  4] = d_GeoData[geo_idx + 5]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx +  5] = 0.0f;

                d_VertexBuffer[vertexfloatBuffer_idx +  6] = d_GeoData[geo_idx + 6]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx +  7] = d_GeoData[geo_idx + 7]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx +  8] = 0.0f;

                d_VertexBuffer[vertexfloatBuffer_idx +  9] = d_GeoData[geo_idx + 8]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 10] = d_GeoData[geo_idx + 9]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 11] = 0.0f;

                d_VertexBuffer[vertexfloatBuffer_idx + 12] = d_GeoData[geo_idx + 10]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 13] = d_GeoData[geo_idx + 11]+renderOffset;
                d_VertexBuffer[vertexfloatBuffer_idx + 14] = 0.0f;

                d_IndexBuffer[indexBuffer_idx     ] = vertexBuffer_idx;
                d_IndexBuffer[indexBuffer_idx +  1] = vertexBuffer_idx + 1;
                d_IndexBuffer[indexBuffer_idx +  2] = vertexBuffer_idx + 2;

                d_IndexBuffer[indexBuffer_idx +  3] = vertexBuffer_idx;
                d_IndexBuffer[indexBuffer_idx +  4] = vertexBuffer_idx + 2;
                d_IndexBuffer[indexBuffer_idx +  5] = vertexBuffer_idx + 3;

                d_IndexBuffer[indexBuffer_idx +  6] = vertexBuffer_idx;
                d_IndexBuffer[indexBuffer_idx +  7] = vertexBuffer_idx + 3;
                d_IndexBuffer[indexBuffer_idx +  8] = vertexBuffer_idx + 4;

                d_IndexBuffer[indexBuffer_idx +  9] = vertexBuffer_idx;
                d_IndexBuffer[indexBuffer_idx + 10] = vertexBuffer_idx + 4;
                d_IndexBuffer[indexBuffer_idx + 11] = vertexBuffer_idx + 1;

                break;
            }

            default:
            {
                break;
            }
        }

    }
}

bool RenderGeometryDirectlyFromCUDA(
        float *       d_GeoData,
        int           geoDataSize,
        float         renderOffset,
        float *       d_pointResultData,
        int           resultDataDim
        )
{
    /*
     * create cuda array to store #tri's each shared data
     */
    int geoUnitCount = geoDataSize / FLOATS_PER_GEO_UNIT;

    NPWCudaDataPointer<int> * triangleCount =
            new NPWCudaDataPointer<int>( geoUnitCount , 0, true,true );

    NPWCudaDataPointer<int> * vertexCount =
            new NPWCudaDataPointer<int>( geoUnitCount , 0, true,true );

    dim3 dimBlock(128,1,1);
    dim3 dimGrid( (int)ceil((float)geoUnitCount/128.0f), 1, 1);

    GetTriangleCountKernel<<<dimGrid, dimBlock>>>(
            triangleCount->GetRawCudaPointer(),
            vertexCount->GetRawCudaPointer(),
            d_GeoData,
            geoUnitCount );

    cudaThreadSynchronize();

    cudaError err;

    if (((err = cudaGetLastError())) != cudaSuccess)
    {
        printf("Error Executing getTriangleCountKernel: %s\n", cudaGetErrorString(err));
        return false;
    }

    /*
     * get offsets for index- and vertexbuffer
     */
    NPWCudaDataPointer<int> * triangleOffsets =
                new NPWCudaDataPointer<int>( geoUnitCount, 0, true,true );

    NPWCudaDataPointer<int> * vertexOffsets =
                new NPWCudaDataPointer<int>( geoUnitCount, 0, true,true );

    CUDPPConfiguration config;
    config.op = CUDPP_ADD;
    config.datatype = CUDPP_INT;
    config.algorithm = CUDPP_SCAN;
    config.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_EXCLUSIVE;

    CUDPPHandle scanplan = 0;
    CUDPPResult result = cudppPlan(&scanplan, config, geoUnitCount, 1, 0);


    if (CUDPP_SUCCESS != result)
    {
        printf("Error creating CUDPPPlan\n");
        exit(-1);
    }

    /*
     * get total number of triangles
     */
    cudppScan( scanplan,
               triangleOffsets->GetRawCudaPointer(),
               triangleCount->GetRawCudaPointer(),
               geoUnitCount );

    cudaThreadSynchronize();

    if (((err = cudaGetLastError())) != cudaSuccess)
    {
        printf("Error Executing Prefix Sum Kernel cudppScan to get triangle offsets: %s\n", cudaGetErrorString(err));
        return false;
    }

    int lastTriangleOffset;
    cudaMemcpy( &lastTriangleOffset,
                triangleOffsets->GetRawCudaPointer() + geoUnitCount - 1 ,
                sizeof(int),
                cudaMemcpyDeviceToHost );

    int lastTriangleCount;
    cudaMemcpy( &lastTriangleCount,
                triangleCount->GetRawCudaPointer() + geoUnitCount - 1 ,
                sizeof(int),
                cudaMemcpyDeviceToHost );

    cudaThreadSynchronize();

    if (((err = cudaGetLastError())) != cudaSuccess)
    {
        printf("Error Copying triangle count from device: %s\n", cudaGetErrorString(err));
        return false;
    }

    int nAllTriangles = lastTriangleOffset + lastTriangleCount;

    /*
     * get total number of vertices
     */
    cudppScan( scanplan,
               vertexOffsets->GetRawCudaPointer(),
               vertexCount->GetRawCudaPointer(),
               geoUnitCount );

    cudaThreadSynchronize();

    if (((err = cudaGetLastError())) != cudaSuccess)
    {
        printf("Error Executing Prefix Sum Kernel cudppScan to get vertex offsets: %s\n", cudaGetErrorString(err));
        return false;
    }

    result = cudppDestroyPlan(scanplan);
    if (CUDPP_SUCCESS != result)
    {
        printf("Error destroying CUDPPPlan\n");
        exit(-1);
    }

    int lastVertexOffset = 0;
    cudaMemcpy( &lastVertexOffset,
                vertexOffsets->GetRawCudaPointer() + geoUnitCount - 1 ,
                sizeof(int),
                cudaMemcpyDeviceToHost );

    int lastVertexCount = 0;
    cudaMemcpy( &lastVertexCount,
                vertexCount->GetRawCudaPointer() + geoUnitCount - 1 ,
                sizeof(int),
                cudaMemcpyDeviceToHost );

    cudaThreadSynchronize();

    if (((err = cudaGetLastError())) != cudaSuccess)
    {
        printf("Error Copying vertex count from device: %s\n", cudaGetErrorString(err));
        return false;
    }

    int nAllVertices = lastVertexOffset + lastVertexCount;

    /*
     * create buffers
     */
    NPWCudaVertexBufferPointer<int> *   triangleIndexBuffer  = 0;
    NPWCudaVertexBufferPointer<float> * triangleVertexBuffer = 0;

    int *   d_triangleIndexBuffer = 0;
    float * d_vertexIndexBuffer   = 0;

    if ( nAllTriangles != 0 )
    {
        // 3 vertices/triangle
        triangleIndexBuffer =
                    new NPWCudaVertexBufferPointer<int>( nAllTriangles * 3 );
        d_triangleIndexBuffer = triangleIndexBuffer->GetBufferAsCudaPointer();
        // 3 floats/vertex
        triangleVertexBuffer =
                    new NPWCudaVertexBufferPointer<float>( nAllVertices * 3 );
        d_vertexIndexBuffer = triangleVertexBuffer->GetBufferAsCudaPointer();

        if ( !d_triangleIndexBuffer || !d_vertexIndexBuffer )
        {
            std::cout << "RenderGeometryDirectlyFromCUDA(): vertex buffer null pointer!" << std::endl;
            if (triangleIndexBuffer)    { triangleIndexBuffer->ReleaseBuffer(); }
            if (triangleVertexBuffer)   { triangleVertexBuffer->ReleaseBuffer(); }

            if (triangleCount)          { delete triangleCount;         triangleCount   = 0; }
            if (vertexCount)            { delete vertexCount;           vertexCount     = 0; }
            if (triangleOffsets)        { delete triangleOffsets;       triangleOffsets = 0; }
            if (vertexOffsets)          { delete vertexOffsets;         vertexOffsets   = 0; }
            if (triangleIndexBuffer)    { delete triangleIndexBuffer;   triangleIndexBuffer  = 0; }
            if (triangleVertexBuffer)   { delete triangleVertexBuffer;  triangleVertexBuffer = 0; }

            std::cout << "RenderGeometryDirectlyFromCUDA() cleaned up stuff" << std::endl;
            return false;
        }
    }

    FillBuffersKernel<<<dimGrid, dimBlock>>>(
            d_vertexIndexBuffer,
            d_triangleIndexBuffer,
            triangleOffsets->GetRawCudaPointer(),
            vertexOffsets->GetRawCudaPointer(),
            d_GeoData,
            geoUnitCount,
            renderOffset,
            d_pointResultData,
            resultDataDim
            );

    cudaThreadSynchronize();

    if (((err = cudaGetLastError())) != cudaSuccess)
    {
        printf("Error Executing FillBuffersKernel: %s\n", cudaGetErrorString(err));
        return false;
    }

    if ( triangleIndexBuffer && triangleVertexBuffer )
    {
        triangleIndexBuffer->ReleaseBuffer();
        triangleVertexBuffer->ReleaseBuffer();

        glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);

        // bind the vertex buffer
        glBindBuffer( GL_ARRAY_BUFFER, triangleVertexBuffer->GetVertexBufferID() );
        glVertexPointer( 3, GL_FLOAT, 3 * sizeof(float), (GLubyte*) (0));
        glEnableClientState(GL_VERTEX_ARRAY);

        // render all the triangles
        glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, triangleIndexBuffer->GetVertexBufferID() );
        glDrawElements( GL_TRIANGLES, 3 * nAllTriangles, GL_UNSIGNED_INT, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
        glBindBuffer(GL_ARRAY_BUFFER,0);
        glPopClientAttrib();

        glFinish();

        CheckOpenGLError();
    }

    /*
     * clean up
     */
    if (triangleCount)          { delete triangleCount;        triangleCount        = 0; }
    if (vertexCount)            { delete vertexCount;          vertexCount          = 0; }
    if (triangleOffsets)        { delete triangleOffsets;      triangleOffsets      = 0; }
    if (vertexOffsets)          { delete vertexOffsets;        vertexOffsets        = 0; }
    if (triangleIndexBuffer)    { delete triangleIndexBuffer;  triangleIndexBuffer  = 0; }
    if (triangleVertexBuffer)   { delete triangleVertexBuffer; triangleVertexBuffer = 0; }

    // done!
    return true;
}
