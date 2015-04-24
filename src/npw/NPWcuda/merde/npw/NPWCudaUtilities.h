/**
 * This file contains general utility functions for dealing with CUDA.
 */

#ifndef NPW_CUDA_UTILITIES_H
#define NPW_CUDA_UTILITIES_H

#include <cuda.h>
#include <cuda_runtime_api.h>

#include <vector_types.h>

//#define CUT_CHECK_ERROR(errorMessage)
//#define CUT_CHECK_ERROR_GL()
//#define CUT_CONDITION( val)
//#define CU_SAFE_CALL_NO_SYNC( call) call
//#define CU_SAFE_CALL( call) call
//#define CUDA_SAFE_CALL_NO_SYNC( call) call
//#define CUDA_SAFE_CALL( call) call
//#define CUT_SAFE_CALL( call) call
//#define CUFFT_SAFE_CALL( call) call

extern "C" void * RawAllocateCUDAArray(int size, void * data);

template <class Type> Type * AllocateCUDAArray(int size, Type * data);

//extern "C" void CUDA_GetMemGetInfo(unsigned int * freeMemory, unsigned int * totalMemory);
extern "C" void CUDA_GetMemGetInfo(size_t * freeMemory, size_t * totalMemory); //DRH g++ 4.8.2; CUDA 6.5

extern "C" void RawClearCUDAArray(void * ptr, int size);

extern "C" void RawCopyCUDAArray(void * source, void * destination, int size, bool toDevice);

template <class Type> void CopyCUDAArray(Type * source, Type * destination, int size, bool toDevice);

extern "C" void FreeCUDAArray(void * array);

#endif // NPW_CUDA_UTILITIES_H

