

#ifndef NPW_CUDA_HISTOGRAM_CREATOR_H
#define NPW_CUDA_HISTOGRAM_CREATOR_H

#define NPW_OPENGL_RENDER

//#define NPW_LOGGING

#include "NPWCudaKernelPrototypes.h"
#include "NPWCudaDataPointer.h"
#include "NPWCudaConstants.h"

#include "NPWCudaOpenGLRenderer.h"

#include <numeric>

/*
 * This class creates a joint-histogram based on two 3D images
 * and returns this or calculates mutual information
 *
 * depends on the cudpp library: http://code.google.com/p/cudpp/
 *
 */
class NPWCudaHistogramCreator
{
public:
  /*
   * constructor
   */
  NPWCudaHistogramCreator();

  /*
   * destructor
   */
  virtual ~NPWCudaHistogramCreator();

  /*
   * set the binsizes by specifying image data, image size and number of bins
   * user should call both methods with the same value of 'bins'
   * (i.e. histogram is square)
   *
   * inputs:
   *      data: pointer to image data
   *      size: size of image data
   *      bins: number of bins in histogram
   *              (n.b. ensure same for 1 and 2 or algorithm will screw up)
   */
  inline void setBinSize1( const float * data, unsigned int size, unsigned int bins )
  {
    setBins(bins);
    setBinSize( data, size, image1Min, binSize1 );
  }
  inline void setBinSize2( const float * data, unsigned int size, unsigned int bins )
  {
    setBins(bins);
    setBinSize( data, size, image2Min, binSize2 );
  }

  /*
   * Synonyms for setBinSize1 and setBinSize2
   */
  inline void setBaseBinSize( const float * data, unsigned int size, unsigned int bins )
  {
    setBinSize1( data, size, bins );
  }
  inline void setSearchBinSize( const float * data, unsigned int size, unsigned int bins )
  {
    setBinSize2( data, size, bins );
  }

  /*
   * renders the joint histogram of two images using CUDA
   *
   * if the binsize is specified here, the binsizes are recalculated,
   * otherwise they must have been set before by calling setBinsize1() and setBinsize2()
   *
   * inputs:
   *      image1 : first image data pointer
   *      image2 : second image data pointer
   *      dimX,Y,Z : dimensions of images (assumed to be same for both images)
   *      bins   : optional, number of bins to use on each axis in histogram
   */
  float * getHistogram( float * image1, float * image2,
      unsigned int dimX, unsigned int dimY, unsigned int dimZ,
      int bins = -1);

  /*
   * calculates mutual information between two images (all on the GPU)
   *
   * if the binsize is specified, the binsizes are (re)calculated,
   * otherwise they must have been set before by calling setBinsize1() and setBinsize2()
   *
   * this method calls getHistogram
   *
   * inputs:
   *      image1 : first image data pointer
   *      image2 : second image data pointer
   *      dimX,Y,Z : dimensions of images (assumed to be same for both images)
   *      bins   : optional, number of bins to use on each axis in histogram
   */
  float getMutualInformation(
      float * image1, float * image2,
      unsigned int dimX, unsigned int dimY, unsigned int dimZ,
      int bins = -1);


private:
  /*
   * CUDA data pointers to the geometry
   *
   * number of elements:   12 (floats per geo unit)
   *                     * 5 ( 5 tets/neighborhood )
   *                     * n (number of image neighborhoods)
   *
   * every geometry unit has 12 floats:
   *      # triangles
   *      height of first vertex (others are 0)
   *      x of first vertex
   *      y of first vertex
   *      and so on till fifth
   */
  NPWCudaDataPointer<float> * geometryData;

  /*
   * images on cuda memory as linear array
   */
  NPWCudaDataPointer<float> * image1Data;
  NPWCudaDataPointer<float> * image2Data;

  /*
   * the renderer renders all geometry (triangles AND points)
   */
  NPWCudaOpenGLRenderer     * renderer;

  /*
   * Clears the CUDA data structures defined above everytime
   * after a render is performed and in the destructor
   */
  void ClearDataStructures();

  /*
   * dimensions of both images
   */
  int imageDims[3];

  /*
   * number of bins in the output histogram (square)
   * setBins() ensure minimum number of bins is 8
   */
  void setBins( unsigned int bins );
  unsigned int nBins;

  /*
   * binsizes and image minima, this is needed to get coordinates in PDF space
   * we want to set these once, so we can keep them constant for several images
   */
  float binSize1; //i.e. units of bins / intensity levels in image 1
  float binSize2; //i.e. units of bins / intensity levels in image 2
  float image1Min; //Minimum value in image1 so that triangles are as close to 0,0, bin as possible
  float image2Min;

  /*
   * based on image data, extracts image minimum and calculates either binsize1 or binsize2
   */
  void setBinSize( const float * imageData, unsigned int dataSize,
      float & imageMin, float & binSize );

  /*
   * create image data on CUDA
   */
  NPWCudaDataPointer<float> * CreateImageData(
      int * imageDims, const char * fileName,
      float & imageMin, float & binSize );
  NPWCudaDataPointer<float> * CreateImageData( float * data, unsigned int size );

  /*
   * allocates memory on CUDA for geometry and fills it with 0.0f
   */
  void CreateGeometryData();

  /*
   * run a kernel to extract the  geometry
   * and lets the renderer render the geometry
   */
  bool RunHistogramKernels( float * histogram );

};

#endif // NPW_CUDA_H
