/*=========================================================================
  Program: MILX milxSimPlugins 1.2
  Module: milxSimNPWTetrahedron.h
  Author:
  Modified by:
  Language: C++
  Created:

  Copyright: (c) CSIRO, Australia

  This software is protected by international copyright laws.
  Any unauthorised copying, distribution or reverse engineering is prohibited.

  Licence:
  All rights in this Software are reserved to CSIRO. You are only permitted
  to have this Software in your possession and to make use of it if you have
  agreed to a Software License with CSIRO.

  BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
=========================================================================*/

#ifndef NPW_OPENGL_RENDERER_H_
#define NPW_OPENGL_RENDERER_H_

#include <GL/glew.h>

#include "Vector2.h"
#include "NPWGLImage.h"
#include "NPWVertices.h"


using namespace npw;

class NPWWindow;
class NPWVertices;

class NPWOpenGLRenderer {

    public:

		/*
		 * constructor
		 */
        NPWOpenGLRenderer( NPWWindow * parentWindow_n );

        /**
         * destructor
         */
        virtual ~NPWOpenGLRenderer();

    	/**
    	 * renders a tetrahedron with a specified weight
    	 */
    	void renderTetrahedron(
    			const Vector2 & v1, const Vector2 & v2, const Vector2 & v3, const Vector2 & v4,
    			double weight_n );

    	/**
    	 * set the histogram for direct point rendering
    	 */
    	void setHistogram( NPWGLImage * histogram_n );

    	/**
    	 * turn shape expansion on and off
    	 */
    	void turnShapeExpansionOn();
    	void turnShapeExpansionOff();
    	bool isShapeExpansionOn();

#ifdef DB_COUNT
        /**
         * print counters of vertices
         */
        void printCounters();
#endif


    private:

        /*
         * renders the current render list in vertices
         */
        void render();

        /*
         * calculate the height of the fan center based on the area
         */
        void calculateHeight();

        /**
         * the histogram image that will eventually hold the histogram
         * render points directly onto this
         */
        NPWGLImage * histogram;

        /**
         * object to hold the vertices
         * why does this need to be a pointer??? ASK NICK!
         */
        NPWVertices * vertices;

        /**
         * the border size of the frame buffer
         * this is set to getMaximumVertexCorrection() in the constructor,
         * though the parent window may set a different border size, we know it does not.
         * if we assume this, we don't need to ask the parent window for it for every shape
         */
        float offset;

        /**
         * shape properties
         */
        double weight;
        float  height;
        double area;

        /**
         * index for a for loop, saves computation time by declaring it here?
         */
        unsigned int i;

        /**
         * store the window for access to the frame buffer size
         */
        NPWWindow * parentWindow;


};

#endif // NPW_OPENGL_RENDERER_H_
