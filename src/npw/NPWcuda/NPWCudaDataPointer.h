
#ifndef NPW_CUDA_DATA_POINTER_H
#define NPW_CUDA_DATA_POINTER_H


#include <cuda.h>
#include <cuda_runtime_api.h>

#include <vector_types.h>

#include <iostream>

#include "NPWCudaUtilities.h"

/**
 * \class NPWCudaDataPointer
 *
 * \brief Wraps memory allocated in cuda memory.
 *
 * This class is a wrapper around data allocated using CUDA. The class knows how
 * to allocate and deallocate itself.
 */
template <class ValueType> class NPWCudaDataPointer
{
    public:

        NPWCudaDataPointer( int size_n, ValueType * data=0,
                            bool allocate=true, bool free_n=true )
            : size(size_n), free(free_n), cudaPtr(0)
        {
//            std::cout << "NPWCudaDataPointer() entering constructor" << std::endl;
//            unsigned int freeMem = 0;
//            unsigned int totalMem = 0;
//            cuMemGetInfo( &freeMem, &totalMem );
//            std::cout << "NPWCudaDataPointer(): creating data pointer of size "
//                      << size * sizeof(ValueType) << " bytes" << std::endl;
//            std::cout << "NPWCudaDataPointer(): free memory before creation : "
//                      << freeMem << " out of a total of " << totalMem
//                      << "\t" << (float)freeMem*100.0f/(float)totalMem << " %" << std::endl;
            if (allocate)
            {
                if (size > 0)
                {
                    int byteSize = size * sizeof(ValueType);

//                    std::cout << "NPWCudaDataPointer() attempting to allocate " << byteSize << " bytes" << std::endl;

                    cudaPtr = (ValueType *) RawAllocateCUDAArray(byteSize, data);
                }
                else
                {
                    cudaPtr = 0;
                }
            }
            else
            {
                cudaPtr = data;
            }

//            cuMemGetInfo( &freeMem, &totalMem );

//            std::cout << "NPWCudaDataPointer(): free memory after creation: " << freeMem << " out of a total of " << totalMem << "\t" << (float)freeMem*100.0f/(float)totalMem << " %" << std::endl;
        }

        virtual ~NPWCudaDataPointer()
        {
//            unsigned int freeMem = 0;
//            unsigned int totalMem = 0;
//            cuMemGetInfo( &freeMem, &totalMem );
//            std::cout << "~NPWCudaDataPointer(): free memory before deletion: " << freeMem << " out of a total of " << totalMem << "\t" << (float)freeMem*100.0f/(float)totalMem << " %" << std::endl;

            if (free && cudaPtr)
            {
                FreeCUDAArray(cudaPtr);

                cudaPtr = 0;
            }

//            cuMemGetInfo( &freeMem, &totalMem );
//            std::cout << "~NPWCudaDataPointer(): free memory after deletion: " << freeMem << " out of a total of " << totalMem << "\t" << (float)freeMem*100.0f/(float)totalMem << " %" << std::endl;
        }

        void ClearData()
        {
            if (cudaPtr)
            {
                RawClearCUDAArray((void*) cudaPtr, size * sizeof(ValueType));
            }
        }

        void UploadData(ValueType * data)
        {
            if (cudaPtr && size > 0)
            {
                RawCopyCUDAArray((void*) data, (void*) cudaPtr,  size * sizeof(ValueType), true);
            }
        }

        ValueType * DownloadData(ValueType * data)
        {
            if (cudaPtr && size > 0)
            {
                if (!data)
                {
                    data = new ValueType[size];
                }

                RawCopyCUDAArray((void*) cudaPtr, (void*) data,  size * sizeof(ValueType), false);
            }

            return data;
        }

        inline unsigned int GetSize(void)           { return size;}
        inline unsigned int GetSizeInBytes(void)    { return size * sizeof(ValueType);}
        inline ValueType * GetRawCudaPointer(void)  { return cudaPtr;}

    protected:
        unsigned int size;
        bool         free;
        ValueType *  cudaPtr;
};


#endif // NPW_CUDA_DATA_POINTER_H
