/*
 * NPWCudaOpenGLRenderer.h
 *
 *  Created on: 22/07/2010
 *      Author: bro86j
 */

#include "NPWCudaConstants.h"
#include "NPWCudaWindow.h"

/*
 * renders all geometry onto the joint histogram
 *      triangles using an OpenGL frame buffer
 *      points on a CUDA data structure
 */
class NPWCudaOpenGLRenderer
{

public:

    /*
     * constructor & destructor
     */
    NPWCudaOpenGLRenderer();
    virtual ~NPWCudaOpenGLRenderer();

    /*
     * render geometry onto specified result
     */
    bool renderGeometryOnResult(
                float * resultData, unsigned int bins,
                float * geoData   , unsigned int geoDataSize );

    int getFrameBufferBorderSize()
    {
        return window->getFrameBufferBorderSize();
    }


private:

    inline int Abs( int a ) { return a < 0 ? -a : a; }

    inline int getArrayIndex2D(
           int x, int y, int xdim, int ydim )
    {
        int clamped_x = ( x < 0 ? Abs(x) : (x >= xdim ? (2*xdim-x-1) : x) );
        int clamped_y = ( y < 0 ? Abs(y) : (y >= ydim ? (2*ydim-y-1) : y) );

        return clamped_x + clamped_y * xdim;
    }

    /*
     * window that does all nasty OpenGL stuff
     */
    NPWCudaWindow * window;

};
