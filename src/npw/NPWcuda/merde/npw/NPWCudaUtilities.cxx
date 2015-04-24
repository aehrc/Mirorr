
#include "NPWCudaUtilities.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

#include <vector_types.h>

#include <stdio.h>

extern "C" void * RawAllocateCUDAArray(int size, void * data)
{
    void * result = 0;
    {
//        fprintf(stderr,"RawAllocateCUDAArray(): trying to allocate %d bytes on CUDA\n",size);
    	cudaError err = cudaMalloc(&result, size);
//    	fprintf(stderr,"RawAllocateCUDAArray(): success?\n");
		if( cudaSuccess != err) {
			fprintf(stderr, "Cuda error allocating %d bytes on device - %s\n", size, cudaGetErrorString(err));
//			unsigned int freeMemory;
//			unsigned int totalMemory;
//			cuMemGetInfo(&freeMemory, &totalMemory);
//			fprintf(stderr, "Memory Usage %d free of %d total memory\n", freeMemory, totalMemory);
			abort();
		}

    }
    if (!result)
    {
    	fprintf(stderr,"Error - cudaMalloc returned null\n");
    	return result;
    }
    if (data)
    {
    	cudaError err = cudaMemcpy(result, data , size, cudaMemcpyHostToDevice);
		if( cudaSuccess != err)
		{
			fprintf(stderr, "Cuda error invoking copying %d bytes to device - source %p destination %p - %s\n", size, data, result, cudaGetErrorString(err));
//			unsigned int freeMemory;
//			unsigned int totalMemory;
//			cuMemGetInfo(&freeMemory, &totalMemory);
//			fprintf(stderr, "Memory Usage %d free of %d total memory\n", freeMemory, totalMemory);
			abort();
		}
		else
		{
    	    err = cudaThreadSynchronize();
    	    if( cudaSuccess != err)
    	    {
    			fprintf(stderr, "Cuda error copying %d bytes to device - source %p destination %p - %s\n", size, data, result, cudaGetErrorString(err));
//    			unsigned int freeMemory;
//    			unsigned int totalMemory;
//    			cuMemGetInfo(&freeMemory, &totalMemory);
//    			fprintf(stderr, "Memory Usage %d free of %d total memory\n", freeMemory, totalMemory);
    			abort();
    	    }
		}
    }
    else
    {
        cudaMemset(result, 0, size);
    }
    return result;
}

extern "C" void RawClearCUDAArray(void * ptr, int size)
{
    cudaMemset(ptr, 0, size);
}


template <class Type> Type * AllocateCUDAArray(int size, Type * data)
{
    int sizeInBytes = size * sizeof(Type);
    return (Type*) RawAllocateCUDAArray(sizeInBytes, (void*)data);
}

extern "C" void RawCopyCUDAArray(void * source, void * destination, int size, bool toDevice)
{
    if (toDevice)
    {
    	cudaError err = cudaMemcpy(destination, source , size, cudaMemcpyHostToDevice);
		if( cudaSuccess != err)
		{
			fprintf(stderr, "Cuda error invoking copying %d bytes to device - source %p destination %p - %s\n", size, source, destination, cudaGetErrorString( err));
//			unsigned int freeMemory;
//			unsigned int totalMemory;
//			cuMemGetInfo(&freeMemory, &totalMemory);
//			fprintf(stderr, "Memory Usage %d free of %d total memory\n", freeMemory, totalMemory);
			abort();
		}
		else
		{
    	    err = cudaThreadSynchronize();
    	    if( cudaSuccess != err)
    	    {
    			fprintf(stderr, "Cuda error copying %d bytes to device - source %p destination %p - %s\n", size, source, destination, cudaGetErrorString( err));
//    			unsigned int freeMemory;
//    			unsigned int totalMemory;
//    			cuMemGetInfo(&freeMemory, &totalMemory);
//    			fprintf(stderr, "Memory Usage %d free of %d total memory\n", freeMemory, totalMemory);
    			abort();
    	    }
		}
    }
    else
    {
    	cudaError err = cudaMemcpy(destination, source , size, cudaMemcpyDeviceToHost);
		if( cudaSuccess != err) {
			fprintf(stderr, "Cuda error invoking copying %d bytes to host - source %p destination %p - %s\n", size, source, destination, cudaGetErrorString( err));
//			unsigned int freeMemory;
//			unsigned int totalMemory;
//			cuMemGetInfo(&freeMemory, &totalMemory);
//			fprintf(stderr, "Memory Usage %d free of %d total memory\n", freeMemory, totalMemory);
			abort();
		}
		else
		{
    	    err = cudaThreadSynchronize();
    	    if( cudaSuccess != err)
    	    {
    			fprintf(stderr, "Cuda error copying %d bytes to host - source %p destination %p - %s\n", size, source, destination, cudaGetErrorString( err));
//    			unsigned int freeMemory;
//    			unsigned int totalMemory;
//    			cuMemGetInfo(&freeMemory, &totalMemory);
//    			fprintf(stderr, "Memory Usage %d free of %d total memory\n", freeMemory, totalMemory);
    			abort();
    	    }
		}

    }
}

template <class Type> void CopyCUDAArray(Type * source, Type * destination, int size, bool toDevice)
{
    int sizeInBytes = size * sizeof(Type);
    RawCopyCUDAArray((void *) source, (void *) destination, sizeInBytes, toDevice);
}

extern "C" void FreeCUDAArray(void * array)
{
    if(array)
    {
//        fprintf(stderr,"FreeCUDAArray(): freeing CUDA array\n");
        cudaFree(array);
        cudaThreadSynchronize();
    }
}

//extern "C" void CUDA_GetMemGetInfo(unsigned int * freeMemory, unsigned int * totalMemory)
extern "C" void CUDA_GetMemGetInfo(size_t * freeMemory, size_t * totalMemory)
{
    cuMemGetInfo(freeMemory, totalMemory);
}

