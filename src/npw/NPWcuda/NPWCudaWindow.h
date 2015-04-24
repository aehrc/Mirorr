/*=========================================================================
  Author:   Daan Broekhuizen
  Date:     May 2010
  Language: C++
=========================================================================*/

#ifndef NPWGL_WINDOW_H
#define NPWGL_WINDOW_H

#include <SDL/SDL.h>
#include <GL/glew.h>

/**
 * \class NPWCudaWindow
 *
 * \brief Handles creating and destroying a window and associated OpenGL objects
 * like a frame buffer and attached texture, shaders etc
 */

class NPWCudaWindow
{
	public:
		NPWCudaWindow(int windowWidth_n, int windowHeight_n, bool visible);
		virtual ~NPWCudaWindow();

		/**
		 * Initialize all OpenGL stuff
		 */
		void initialize( unsigned int width, unsigned int height );

		/**
		 * Restore all OpenGL stuff
		 */
		void finalize();

		/**
		 * Close the window
		 */
		void destroy();

		/**
		 * reads back the OpenGL frame buffer from the window class and adds this to result
		 */
        bool readBackImage( float * result, unsigned int width, unsigned int height );


	private:
        inline int Abs( int a )         { return a < 0 ? -a : a; }
        inline int Max( int a, int b )  { return ( a > b ? a : b ); }
        inline int Min( int a, int b )  { return ( a < b ? a : b ); }

        inline int getArrayIndex2D(
               int x, int y, int xdim, int ydim )
        {
            int clamped_x = ( x < 0 ? Abs(x) : (x >= xdim ? (2*xdim-x-1) : x) );
//            clamped_x = Max(0, clamped_x);
//            clamped_x = Min(xdim-1, clamped_x);
            int clamped_y = ( y < 0 ? Abs(y) : (y >= ydim ? (2*ydim-y-1) : y) );
//            clamped_y = Max(0, clamped_y);
//            clamped_y = Min(ydim-1, clamped_y);

            return clamped_x + clamped_y * xdim;
        }

	    /**
	     * Initialize glew
	     */
		bool initializeGlew();

		/**
		 * Create the window
		 */
		bool createWindow();

		/**
		 * Enter and leave orthographic projection mode
		 */
		void enterOrthoMode();
		void leaveOrthoMode();

		/**
		 * Frame buffer stuff
		 */
		GLuint frameBufferId;

		int frameBufferWidth;
		int frameBufferHeight;
		int frameBufferBorderSize;

		bool createFrameBuffer(int width, int height);
		void destroyFrameBuffer();
		void clearBuffer();

		void bindFrameBuffer();
	    void unBindFrameBuffer();

public:
	    /**
         * Access to frame buffer size
         */
        inline int getFrameBufferWidth()      { return frameBufferWidth;  }
        inline int getFrameBufferHeight()     { return frameBufferHeight; }
        inline int getFrameBufferBorderSize() { return frameBufferBorderSize; }
        inline int getFrameBufferCount()      { return 4 * frameBufferWidth * frameBufferHeight; }

private:
	    /**
	     * Texture stuff
	     */
		GLuint textureId;

		void bindTexture();
	    void unbindTexture();

	    GLfloat * readBackTexture();


		/**
		 * Shader stuff
		 */
		SDL_Surface * sdlSurface;

	    GLuint vertexShader;
	    GLuint fragmentShader;
	    GLuint program;

	    const char * fragmentText;
	    const char * vertexText;

	    bool compileShader();
	    void associateShader();
	    void disassociateShader();


	    /**
	     * Size of the window
	     */
		int windowWidth;
		int windowHeight;

		/**
		 * Should the window be hidden upon creation
		 */
		bool visible;

		/**
		 * was a framebuffer created?
		 */
		bool hasCreatedFrameBuffer;
};

#endif // NPWGL_WINDOW_H

