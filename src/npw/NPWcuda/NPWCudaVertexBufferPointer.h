/*
 * NPWCudaVertexBuffer.h
 *
 *  Created on: 28/07/2010
 *      Author: bro86j
 */

#ifndef NPW_CUDA_VERTEX_BUFFER_POINTER_H_
#define NPW_CUDA_VERTEX_BUFFER_POINTER_H_

#include <GL/glew.h>
#include <cuda_gl_interop.h>
#include <iostream>

bool CheckOpenGLError()
{
    GLenum errCode = glGetError();
    const GLubyte *errString;

    if (errCode != GL_NO_ERROR)
    {
        errString = gluErrorString(errCode);
        std::cout << "OpenGL error detected: " << errCode << " reason: " << errString << std::cout;
        return false;
    }
    else
    {
        return true;
    }
}

template <class ValueType> class NPWCudaVertexBufferPointer
{
    public:

        NPWCudaVertexBufferPointer(
                int size, ValueType * data=0 )
            : size(size)
        {
//            unsigned int freeMem = 0;
//            unsigned int totalMem = 0;
//            cuMemGetInfo( &freeMem, &totalMem );
//            std::cout << "NPWCudaVertexBufferPointer(): free memory before creation: " << freeMem << " out of a total of " << totalMem << "\t" << (float)freeMem*100.0f/(float)totalMem << " %" << std::endl;
//            std::cout << "NPWCudaVertexBufferPointer(): creating vertexbuffer of size " << size * sizeof(ValueType) << " bytes" << std::endl;


            glGenBuffers( 1, &vertexBufferID );

            glBindBuffer( GL_ARRAY_BUFFER, vertexBufferID );

            int sizeInBytes = size * sizeof(ValueType);

            glBufferData(
                    GL_ARRAY_BUFFER,
                    sizeInBytes,
                    data,
                    GL_DYNAMIC_DRAW );

            glBindBuffer( GL_ARRAY_BUFFER, 0 );

//            std::cout << "NPWCudaVertexBufferPointer() created vertex buffer: " << vertexBufferID << std::endl;

//            std::cout << "flushing..." << std::flush;

            glFinish();

//            std::cout << " DONE!" << std::endl;

            CheckOpenGLError();

//            cuMemGetInfo( &freeMem, &totalMem );

//            std::cout << "NPWCudaVertexBufferPointer(): free memory after creation: " << freeMem << " out of a total of " << totalMem << "\t" << (float)freeMem*100.0f/(float)totalMem << " %" << std::endl;
        }

        virtual ~NPWCudaVertexBufferPointer()
        {
//            std::cout << "~NPWCudaVertexBufferPointer() deleting vertex buffer: " << vertexBufferID << std::endl;
//            unsigned int freeMem = 0;
//            unsigned int totalMem = 0;
//            cuMemGetInfo( &freeMem, &totalMem );
//            std::cout << "~NPWCudaVertexBufferPointer(): free memory before deletion: " << freeMem << " out of a total of " << totalMem << "\t" << (float)freeMem*100.0f/(float)totalMem << " %" << std::endl;

            glDeleteBuffers(1, &vertexBufferID);
            vertexBufferID = 0;

            glFinish();

            CheckOpenGLError();

//            cuMemGetInfo( &freeMem, &totalMem );
//            std::cout << "~NPWCudaVertexBufferPointer(): free memory after deletion: " << freeMem << " out of a total of " << totalMem << "\t" << (float)freeMem*100.0f/(float)totalMem << " %" << std::endl;

        }

        ValueType * GetBufferAsCudaPointer()
        {
            cudaError err;

            void * result;
            // todo: eventually, change these deprecated calls
            err = cudaGLRegisterBufferObject(vertexBufferID);
            if ( err != cudaSuccess )
            {
                std::cout << "Quick Error registering openGL buffer (ID: " << vertexBufferID << ") as CUDA object: " << cudaGetErrorString(err) << std::endl;
                return 0;
            }

            if (((err = cudaGetLastError())) != cudaSuccess)
            {
                std::cout << "Error registering openGL buffer (ID: " << vertexBufferID << ") as CUDA object: " << cudaGetErrorString(err) << std::endl;
                return 0;
            }
//            cudaGraphicsGLRegisterBuffer(&cuda_vbo_resource, vertexBufferID, cudaGraphicsMapFlagsWriteDiscard);
            err = cudaGLMapBufferObject(&result, vertexBufferID);
            if ( err != cudaSuccess )
            {
                std::cout << "Quick Error binding openGL buffer (ID: " << vertexBufferID << ") as CUDA object: " << cudaGetErrorString(err) << std::endl;
                return 0;
            }
//            cudaGraphicsMapResources(1, vbo_resource, 0);
    
            if (((err = cudaGetLastError())) != cudaSuccess)
            {
                std::cout << "Error binding openGL buffer (ID: " << vertexBufferID << ") as CUDA object: " << cudaGetErrorString(err) << std::endl;
                return 0;
            }

            CheckOpenGLError();

            return (ValueType*)result;
        }

        void ReleaseBuffer()
        {
            // todo: eventually, change these deprecated calls
            cudaGLUnmapBufferObject(vertexBufferID);
//            cudaGraphicsUnmapResources(1, vbo_resource, 0)
            cudaGLUnregisterBufferObject(vertexBufferID);
//            cudaGraphicsUnregisterResource(vbo_res);

            CheckOpenGLError();
        }


        inline int          GetVertexBufferID()     { return vertexBufferID; }

        inline unsigned int GetSize(void)           { return size;}
        inline unsigned int GetSizeInBytes(void)    { return size * sizeof(ValueType);}

    protected:
        unsigned int size;
        GLuint       vertexBufferID;
};

#endif // NPW_CUDA_VERTEX_BUFFER_POINTER_H_
