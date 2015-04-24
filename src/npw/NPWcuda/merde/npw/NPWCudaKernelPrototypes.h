
#ifndef NPW_CUDA_KERNEL_PROTOTYPES_H
#define NPW_CUDA_KERNEL_PROTOTYPES_H

/*
 * this will run a thread for every neighborhood,
 * extracting geometry and saving in geometryData
 *
 * Note : geometry size is implicit in image dimensions, so not needed here
 */
bool ExtractGeometryFromImagesCUDA(
        unsigned int  xDim,            //Image dimensions
        unsigned int  yDim,
        unsigned int  zDim,
        float *       d_geometryData,  //Device geometryData pointer
        float *       d_image1Data,    //Device image data pointer
        float *       d_image2Data,    //Device image data pointer
        float         image1Min,       //Minimum values
        float         image2Min,
        float         binSize1,        //Bin size in number of intensities
        float         binSize2
        );

/*
 * This renders the geometry onto a joint histogram (frame buffer)
 */
bool RenderGeometryDirectlyFromCUDA(
        float *       d_GeoData,        //Device geometry data pointer (triangles)
        int           geoDataSize,      //Size of data
        float         renderOffset,     //Thickness of border on frame buffer (to catch triangles that go slightly over the edge)
        float *       d_PointResultData,//Device pointer to result of rendering points
        int           resultDataDim     //Size of points result (bins*bins) as treated as linear array
        );


#endif // NPW_CUDA_KERNEL_PROTOTYPES_H
