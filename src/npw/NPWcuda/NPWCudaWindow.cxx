/*=========================================================================
  Author:   Daan Broekhuizen
  Date:     May 2010
  Language: C++
=========================================================================*/

#include <stdexcept>
#include <iostream>
#include <GL/glew.h>

#include "NPWCudaWindow.h"
#include "NPWStopwatch.h"
#include "NPWCudaConstants.h"

//#define DB_STATUS

NPWCudaWindow::NPWCudaWindow(int windowWidth_n, int windowHeight_n, bool visible_n) :
    frameBufferId(0), frameBufferWidth(0), frameBufferHeight(0),
    frameBufferBorderSize ( (int)(MAXIMUM_VERTEX_CORRECTION+1.0f) ),
    textureId(0), sdlSurface(0),
    windowWidth(windowWidth_n), windowHeight(windowHeight_n), visible(visible_n),
    hasCreatedFrameBuffer(false)
{
    createWindow();

    if ( !GLEW_VERSION_2_0 )
    {
        printf("Your video card (%s %s) supports OpenGL version \"%s\", but NPW rendering requires 2.0.\n",
                glGetString(GL_VENDOR),
                glGetString(GL_RENDERER),
                glGetString(GL_VERSION) );
        destroy();
        throw std::runtime_error("Quit program because OpenGL 2.0 is not supported.");
    }

    vertexText =
       "varying float height; \
        void main(void) \
        { \
            height = gl_Vertex.z; \
            gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xy, 0.0, 1.0); \
        }";

    fragmentText =
       "varying float height; \
        void main(void) \
        { \
            gl_FragColor = vec4(height, 0.0, 0.0, 0.0); \
        }";

    compileShader();
}

NPWCudaWindow::~NPWCudaWindow()
{
    destroy();
}

void NPWCudaWindow::destroy()
{
    destroyFrameBuffer();
    sdlSurface = 0;
    SDL_Quit(); // whoa check this
}

bool NPWCudaWindow::createWindow()
{
    /**
     * Initialize SDL
     */
    if (SDL_Init(SDL_INIT_VIDEO) != 0)
    {
        std::cout << "unable to initialize SDL: " << SDL_GetError() << ", exiting" << std::endl;
        return false;
    }

    SDL_WM_SetCaption("NPW Renderer", "NPW Renderer");

    int flags = SDL_OPENGL;

    SDL_WM_GrabInput( SDL_GRAB_OFF);



    /**
     * Set up double buffering and vertical synch to prevent tearing
     */
    SDL_GL_SetAttribute(SDL_GL_SWAP_CONTROL, 1);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

    /**
     * 24 bit depth buffer
     */
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

    /**
     * Create the window
     */
    sdlSurface = SDL_SetVideoMode(windowWidth, windowHeight, 24, flags);

    if (!sdlSurface)
    {
        std::cout << "unable to set video mode: " << SDL_GetError() << std::endl;
        return false;
    }

    if (!initializeGlew())
    {
        std::cout << "failed to initialize GLEW!" << std::endl;
        return false;
    }

    glClearDepth(1.0); // Enables Clearing Of The Depth Buffer
    glDepthFunc(GL_LESS); // The Type Of Depth Test To Do
    glEnable(GL_DEPTH_TEST); // Enables Depth Testing
    glShadeModel(GL_SMOOTH); // Enables Smooth Color Shading
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

    if (!visible)
    {
        SDL_WM_IconifyWindow();
    }

    return true;
}

bool NPWCudaWindow::createFrameBuffer(int width, int height)
{
    if ( frameBufferWidth  == width  &&
         frameBufferHeight == height &&
         hasCreatedFrameBuffer )
    {
        // return true, since the required frame buffer is already there
        return true;
    }

    frameBufferWidth  = width;
    frameBufferHeight = height;

    destroyFrameBuffer();

    /**
     * Create a new texture to render the frame buffer into
     */
    glGenTextures(1, &textureId);

    glBindTexture(GL_TEXTURE_2D, textureId);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA_FLOAT32_ATI, width, height, 0,
            GL_RGBA, GL_FLOAT, NULL);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    /**
     * Create a frame buffer and attach the texture to it
     */
    glGenFramebuffersEXT(1, &frameBufferId);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frameBufferId);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, textureId, 0);

    hasCreatedFrameBuffer = true;

    /**
     * Check the status of the frame buffer
     */
    GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);


    bool result = true;

    switch (status)
    {
        case GL_FRAMEBUFFER_COMPLETE_EXT:
        {
            break;
        }
        case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
        {
            std::cout << "Frame buffer is not supported!" << std::endl;
            result = false;
            break;
        }
        default:
        {
            std::cout << "Frame buffer programming error, check configuration: " << status << std::endl;
            result = false;
        }
    }

    /**
     * Unbind the buffer
     */
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

    if (!result)
    {
        std::cout << "Failed to create frame buffer" << std::endl;
    }

    // Check opengl error

    return result;
}

bool NPWCudaWindow::initializeGlew()
{
    GLenum err = glewInit();

    if (GLEW_OK != err)
    {
        std::cout << "failed to initialize GLEW: " << glewGetErrorString(err) << std::endl;
        return false;
    }

    return true;
}

void NPWCudaWindow::clearBuffer()
{
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frameBufferId);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

void NPWCudaWindow::bindFrameBuffer()
{
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frameBufferId);
}

void NPWCudaWindow::unBindFrameBuffer()
{
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

void NPWCudaWindow::bindTexture()
{
    glEnable(GL_TEXTURE_2D);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, textureId);
}

void NPWCudaWindow::unbindTexture()
{
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, 0);
}

GLfloat * NPWCudaWindow::readBackTexture()
{
    int frameBufferCount = getFrameBufferCount();
//    std::cout << "NPWCudaWindow::readBackTexture() creating resultbuffer with " << frameBufferCount << " elements" << std::endl;
    GLfloat * resultBuffer = new GLfloat [frameBufferCount];

    bindTexture();

    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, resultBuffer);

    unbindTexture();

    return resultBuffer;
}

bool NPWCudaWindow::readBackImage( float * result, unsigned int width, unsigned int height )
{
    GLfloat * data = readBackTexture();

    if ( !data )
    {
        return false;
    }

    /*
     * add the frame buffer to the result while folding the edges in
     * add, because points are rendered directly onto this result by the openGL renderer
     */
    int outputX    = 0;
    int outputY    = 0;
    int result_idx = 0;
    int data_idx   = 0;

//    std::cout << "NPWCudaWindow::readBackImage() result should be " << width << " x " << height << std::endl;
//
//    std::cout << "  x\t  y\t  idx\n" << std::endl;

    for ( int y = 0; y < frameBufferHeight; ++y )
    {
        for ( int x = 0; x < frameBufferWidth; ++x )
        {
            outputX = x - frameBufferBorderSize;
            outputY = y - frameBufferBorderSize;

            result_idx = getArrayIndex2D(outputX,outputY, width, height);
            data_idx   = 4 * ( x + frameBufferWidth * y );
//            std::cout << "  " << x << "\t  " << y << "\t  " << result_idx
//                      << ((result_idx < 0 || result_idx >= (int)(width*height))?"   <---":"")
//                      << std::endl;
//            std::cout << "NPWCudaWindow::readBackImage() result idx : " << result_idx << std::endl;
//            std::cout << "NPWCudaWindow::readBackImage() data   idx : " << data_idx << std::endl;

            result[result_idx] += data[data_idx];
        }
    }

    /*
     * delete the data array, which we created
     */
    if ( data )
    {
        delete [] data;
        data = 0;
    }

    return true;
}

void NPWCudaWindow::enterOrthoMode()
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glFrontFace(GL_CCW);
    glDisable(GL_CULL_FACE);
    glDisable(GL_POLYGON_OFFSET_FILL);
    glDepthMask(GL_FALSE);
    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT, GL_FILL);
    glPolygonMode(GL_BACK, GL_FILL);

    glBlendEquation( GL_FUNC_ADD );
    glEnable( GL_BLEND );
    glBlendFunc( GL_ONE, GL_ONE );

    glDisable(GL_DEPTH_TEST);


    glViewport(0, 0, frameBufferWidth, frameBufferHeight);

    // Push OpenGL state
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();

    // Set up projection for 2D display
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0f, frameBufferWidth, 0.0f, frameBufferHeight);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void NPWCudaWindow::leaveOrthoMode()
{
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    // Pop OpenGL state
    glPopAttrib();
}

bool NPWCudaWindow::compileShader()
{
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

    glShaderSource(vertexShader, 1, (const char **) &vertexText, NULL);
    glShaderSource(fragmentShader, 1, (const char **) &fragmentText, NULL);

    glCompileShader(vertexShader);
    glCompileShader(fragmentShader);

    program = glCreateProgram();

    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);

    glLinkProgram(program);

    return true;
}

void NPWCudaWindow::associateShader()
{
    glUseProgram(program);
}

void NPWCudaWindow::disassociateShader()
{
    glUseProgram(0);
}

void NPWCudaWindow::initialize( unsigned int width, unsigned int height )
{
    createFrameBuffer( width  + 2 * frameBufferBorderSize,
                       height + 2 * frameBufferBorderSize );
    clearBuffer();
    bindFrameBuffer();
    associateShader();
    enterOrthoMode();
}

void NPWCudaWindow::finalize()
{
    leaveOrthoMode();
    disassociateShader();
    unBindFrameBuffer();
}

void NPWCudaWindow::destroyFrameBuffer()
{
    if (frameBufferId)
    {
//        std::cout << "NPWCudaWindow::destroyFrameBuffer() Deleting frame buffer: " << frameBufferId << std::endl;
        glDeleteFramebuffersEXT(1, &frameBufferId);
        frameBufferId = 0;
    }

    if (textureId)
    {
//        std::cout << "NPWCudaWindow::destroyFrameBuffer() Deleting texture: " << textureId << std::endl;
        glDeleteTextures(1, &textureId);
        textureId = 0;
    }

    hasCreatedFrameBuffer = false;



}
