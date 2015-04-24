/*=========================================================================
  Author:   Daan Broekhuizen
  Date:     May 2010
  Language: C++
=========================================================================*/

#include <algorithm>

#ifndef NPWGL_IMAGE_H
#define NPWGL_IMAGE_H

/**
 * \class NPWGLImage
 *
 * \brief Supports loading and manipulating 3D images as float arrays
 */
class NPWGLImage
{
  private:

    bool LoadData(const char * fileName);

    int *   imageDimenions;
    int     voxelCount;
    float * data;

    bool    verbosity;


  public:
    /**
     * constructors and destructor
     */
    NPWGLImage(const char * fileName, bool verbose = false);
    NPWGLImage(int xDim = 1, int yDim = 1, int zDim = 1, bool verbose = false);
    NPWGLImage(float * imageData, int xDim = 1, int yDim = 1, int zDim = 1, bool verbose = false);
    virtual ~NPWGLImage();

    /**
     * set the verbosity of this image
     */
    void SetVerbosity( bool verbose ) { verbosity = verbose; }

    /**
     * access to image values
     */
    inline float & operator()( int x, int y = 0, int z = 0 )    {
      return data[x + imageDimenions[0] * ( y + z * imageDimenions[1] )];
    }
    inline float & operator()( int x, int y = 0, int z = 0 ) const
    {
      return data[x + imageDimenions[0] * ( y + z * imageDimenions[1] )];
    }

    inline float & GetValueAtIndex( int index )                 { return data[index]; }
    inline float & GetValueAtIndex( int index ) const           { return data[index]; }

    void SetValueAtIndex(const int index, const float value)    { data[index] = value;  }

    /**
     * access to image statistics
     */
    void GetDimensions(int & x, int & y, int & z) const
    {
      x = imageDimenions[0];
      y = imageDimenions[1];
      z = imageDimenions[2];
    }
    inline int dim(int index) const                             { return imageDimenions[index]; }
    inline int GetVoxelCount() const                            { return voxelCount; }
    inline int GetNeighborhoods() const                         { return ( imageDimenions[0] - 1 ) *
                                                                         ( imageDimenions[1] - 1 ) *
                                                                         ( imageDimenions[2] - 1 ); }
    inline bool HasData() const                                 { return voxelCount > 0; }
    inline float * GetData() const                              { return data; }
    void clear()
    {
        for ( int i = 0; i < voxelCount; ++i )
        {
            data[i] = 0;
        }
    }

    inline float getMaximum() const { return *std::max_element( data, data+voxelCount ); }
    inline float getMinimum() const { return *std::min_element( data, data+voxelCount ); }

    /**
     * Save functionality
     */
    bool SaveRawDataToText( const char * fileName ) const;

#ifdef USE_DEVIL
    enum NormalizationMode {
        NO_NORMALIZATION,
        LINEAR_SCALE,
        LOG_SCALE,
    };
    bool SaveToFile(const char * fileName, NormalizationMode normalizationMode,
            float imageMinimum = 0.0f, float imageMaximum = 255.0f) const;
    void normaliseData(
            float * data, const int dataSize,
            const float minReq = 0.0f, const float maxReq = 255.0f,
            const bool useLogScale = false) const;
#endif

};

#endif //MILX_SIM_NPWGL_IMAGE_H
