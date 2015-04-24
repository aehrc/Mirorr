/*=========================================================================
  Author:   Daan Broekhuizen
  Date:     May 2010
  Language: C++
=========================================================================*/

#ifndef NPWGL_WINDOW_H
#define NPWGL_WINDOW_H

#include "NPWGLImage.h"
#include "NPWOpenGLRenderer.h"

#include <SDL/SDL.h>
#include <GL/glew.h>

/**
 * \class NPWWindow
 *
 * \brief Handles creating and destroying a window and associated OpenGL objects
 * like a frame buffer and attached texture, shaders etc
 */

class NPWOpenGLRenderer;

class NPWWindow
{
    private:
        /**
         * The renderer
         */
        NPWOpenGLRenderer * renderer;

	public:
		NPWWindow(int windowWidth_n, int windowHeight_n, bool visible);
		virtual ~NPWWindow();

		/**
		 * Access to the renderer
		 */
		NPWOpenGLRenderer * getRenderer() { return renderer; }

		/**
		 * Initialize all OpenGL stuff
		 */
		void initialize( NPWGLImage * histogram );

		/**
		 * Restore all OpenGL stuff
		 */
		void finalize();

		/**
		 * Close the window
		 */
		void destroy();

		/**
		 * returns the texture as an NPWGLImage
		 */
        bool readBackImage( NPWGLImage * result );


	private:
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

		int frameBufferHeight;
		int frameBufferWidth;
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
