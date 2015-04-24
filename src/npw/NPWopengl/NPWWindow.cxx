/*=========================================================================
  Author:   Daan Broekhuizen
  Date:     May 2010
  Language: C++
=========================================================================*/

#include <stdexcept>
#include <iostream>
#include <GL/glew.h>

#include "NPWWindow.h"
#include "NPWStopwatch.h"
#include "NPWGeometry.h"

//#define DB_STATUS

NPWWindow::NPWWindow(int windowWidth_n, int windowHeight_n, bool visible_n) :
    frameBufferId(0), textureId(0), sdlSurface(0),
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
    // was: vec4(gl_Vertex.xy, 0.0, 1.0)

    fragmentText =
       "varying float height; \
        void main(void) \
        { \
            gl_FragColor = vec4(height, 0.0, 0.0, 0.0); \
        }";
    // was vec4(1.0, 0.0, 0.0, height)

    compileShader();

    renderer = new NPWOpenGLRenderer( this );

#ifdef DB_STATUS
    std::cout << "NPW Window created\n"
              << "   status messages ON" << std::endl;
#endif
}

NPWWindow::~NPWWindow()
{
#ifdef DB_STATUS
    std::cout << "destroying NPW window" << std::endl;
#endif

    destroy();

    if ( renderer )
    {
        delete renderer;
    }
}

void NPWWindow::destroy()
{
    sdlSurface = 0;
    SDL_Quit(); // whoa check this
}

bool NPWWindow::createWindow()
{
#ifdef DB_STATUS
    std::cout << "Creating a window" << std::endl;
#endif

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

#ifdef DB_STATUS
    std::cout << "setting rendering window properties" << std::endl;
#endif

    glClearDepth(1.0); // Enables Clearing Of The Depth Buffer
    glDepthFunc(GL_LESS); // The Type Of Depth Test To Do
    glEnable(GL_DEPTH_TEST); // Enables Depth Testing
    glShadeModel(GL_SMOOTH); // Enables Smooth Color Shading
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

    if (!visible)
    {
#ifdef DB_STATUS
        std::cout << "hiding window" << std::endl;
#endif
        SDL_WM_IconifyWindow();
    }

#ifdef DB_STATUS
    std::cout << "window initialization complete" << std::endl;
#endif

    return true;
}

bool NPWWindow::createFrameBuffer(int width, int height)
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

#ifdef DB_STATUS
    std::cout << "Got new frame buffer id: " << frameBufferId << " width: " << width << " height: " << height
              << " border size : " << frameBufferBorderSize << std::endl;
#endif

    /**
     * Check the status of the frame buffer
     */
    GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);


    bool result = true;

    switch (status)
    {
        case GL_FRAMEBUFFER_COMPLETE_EXT:
        {
#ifdef DB_STATUS
            std::cout << "Frame buffer creation was successful" << std::endl;
#endif
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
    else
    {
#ifdef DB_STATUS
        std::cout << "Successfully initialised buffer!" << std::endl;
#endif
    }

    return result;
}

bool NPWWindow::initializeGlew()
{
    GLenum err = glewInit();

    if (GLEW_OK != err)
    {
        std::cout << "failed to initialize GLEW: " << glewGetErrorString(err) << std::endl;
        return false;
    }

#ifdef DB_STATUS
    std::cout << "GLEW initialization complete" << std::endl;
#endif

    return true;
}

void NPWWindow::clearBuffer()
{
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frameBufferId);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

void NPWWindow::bindFrameBuffer()
{
#ifdef DB_STATUS
    std::cout << "Binding buffer" << std::endl;
#endif

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frameBufferId);
}

void NPWWindow::unBindFrameBuffer()
{
#ifdef DB_STATUS
    std::cout << "Unbinding buffer" << std::endl;
#endif
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

void NPWWindow::bindTexture()
{
    glEnable(GL_TEXTURE_2D);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, textureId);
}

void NPWWindow::unbindTexture()
{
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, 0);
}

GLfloat * NPWWindow::readBackTexture()
{
#ifdef DB_STATUS
    std::cout << "Reading data from OpenGL" << std::endl;
#endif

    int frameBufferCount = getFrameBufferCount();
    GLfloat * resultBuffer = new GLfloat[frameBufferCount];

    bindTexture();

    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, resultBuffer);

    unbindTexture();

    return resultBuffer;
}

bool NPWWindow::readBackImage(NPWGLImage * result)
{
    NPWStopwatch stopwatch;

    GLfloat * data = readBackTexture();

#ifdef DB_STATUS
    std::cout << "reading back from frame buffer took "
              << stopwatch.reset() << " s" << std::endl;
#endif

    /*
     * add the frame buffer to the histogram while folding the edges in
     * add, because points are rendered directly onto this result by the openGL renderer
     */

    int i = 0;
    for ( int y = 0; y < frameBufferHeight; ++y )
    {
        for ( int x = 0; x < frameBufferWidth; ++x )
        {
            int outputX = x - frameBufferBorderSize;
            int outputY = y - frameBufferBorderSize;

            /*
             * fold x
             */
            if ( outputX < 0 )
            {
                outputX = -outputX - 1;
            }
            else if ( outputX >= result->dim(0) )
            {
                outputX = 2 * result->dim(0) - outputX - 1;
            }

            /*
             * fold y
             */
            if ( outputY < 0 )
            {
                outputY = -outputY - 1;
            }
            else if ( outputY >= result->dim(1) )
            {
                outputY = 2 * result->dim(1) - outputY - 1;
            }

            /*
             * add frame buffer to result
             */
            i = x + frameBufferWidth * y;

            (*result)(outputX,outputY) += data[i*4];
        }
    }

    /*
     * delete the data array, which we created
     */
    delete [] data;
    data = 0;

    #ifdef DB_STATUS
    std::cout << "converted GLFloat array to NPW Image in "
              << stopwatch.reset() << " s" << std::endl;
#endif

    return true;
}

void NPWWindow::enterOrthoMode()
{
#ifdef DB_STATUS
    std::cout << "Entering orthographic projection mode and setting OpenGL blending modes" << std::endl;
#endif

    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glFrontFace(GL_CCW);
    glDisable(GL_CULL_FACE);
    glDisable(GL_POLYGON_OFFSET_FILL);
    glDepthMask(GL_FALSE);
    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT, GL_FILL);
    glPolygonMode(GL_BACK, GL_FILL);

//    glBlendEquationSeparate( GL_FUNC_ADD, GL_FUNC_ADD );
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

void NPWWindow::leaveOrthoMode()
{
#ifdef DB_STATUS
    std::cout << "Leaving orthographic projection mode" << std::endl;
#endif

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    // Pop OpenGL state
    glPopAttrib();
}

bool NPWWindow::compileShader()
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

#ifdef DB_STATUS
    std::cout << "Compiled the shader" << std::endl;
#endif

    return true;
}

void NPWWindow::associateShader()
{
#ifdef DB_STATUS
    std::cout << "Associating GLSL shader" << std::endl;
#endif
    glUseProgram(program);
}

void NPWWindow::disassociateShader()
{
#ifdef DB_STATUS
    std::cout << "Disassociating GLSL shader" << std::endl;
#endif
    glUseProgram(0);
}

void NPWWindow::initialize( NPWGLImage * histogram )
{
    frameBufferBorderSize = NPWGeometryMath::Ceil( MAXIMUM_VERTEX_CORRECTION );

    createFrameBuffer( histogram->dim(0) + 2 * frameBufferBorderSize,
                       histogram->dim(1) + 2 * frameBufferBorderSize);
    clearBuffer();
    bindFrameBuffer();
    associateShader();
    enterOrthoMode();
    renderer->setHistogram( histogram );
}

void NPWWindow::finalize()
{
    leaveOrthoMode();
    disassociateShader();
    unBindFrameBuffer();
}
