/*=========================================================================
  Program: MILX milxSimPlugins 1.2
  Module: NPWOpenGLRenderer.cxx
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

#include "NPWOpenGLRenderer.h"

#include "NPWGeometry.h"
#include "NPWWindow.h"
#include "NPWConstants.h"

#include <algorithm>


NPWOpenGLRenderer::NPWOpenGLRenderer( NPWWindow * parentWindow_n ) :
        histogram( 0 ),
        parentWindow ( parentWindow_n )
{
    vertices = new NPWVertices();

    vertices->setParentWindow( parentWindow_n );

    offset = std::ceil( MAXIMUM_VERTEX_CORRECTION );

}


NPWOpenGLRenderer::~NPWOpenGLRenderer()
{
    delete vertices;
}


/**
 *  Calculates the height of the first vertex.
 *  weight must be set and area calculated.
 */
void NPWOpenGLRenderer::calculateHeight()
{
    switch ( vertices->n_renderList_size )
    {
        /*
         * points
         */
        case 1:
        {
            height = (float)weight;
            break;
        }

        /*
         * triangles with a vertex pair and unsafe borders
         */
        case 3:
        {
            area = 0.5 * NPWGeometryMath::Abs(
                              ( (double)(vertices->renderList[1]->x - vertices->renderList[0]->x) *
                                (double)(vertices->renderList[2]->y - vertices->renderList[0]->y) ) -
                              ( (double)(vertices->renderList[2]->x - vertices->renderList[0]->x) *
                                (double)(vertices->renderList[1]->y - vertices->renderList[0]->y) ) );
            height = (float)( 3.0 * weight / NPWGeometryMath::Max(area, MINIMUM_AREA) );
            break;
            break;
        }

        /*
         * triangles with three subtriangles
         * and
         * triangles with a vertex on the border and unsafe borders
         */
        case 4:
        case 5:
        {
            //    area = 0.5f * abs( ( (x1-x0)*(y2-y0) ) - ( (x2-x0)*(y1-y0) ) );
            area = 0.5 * NPWGeometryMath::Abs(
                              ( (double)(vertices->renderList[2]->x - vertices->renderList[1]->x) *
                                (double)(vertices->renderList[3]->y - vertices->renderList[1]->y) ) -
                              ( (double)(vertices->renderList[3]->x - vertices->renderList[1]->x) *
                                (double)(vertices->renderList[2]->y - vertices->renderList[1]->y) ) );
            height = (float)( 3.0 * weight / NPWGeometryMath::Max(area, MINIMUM_AREA) );
            break;
        }

        /*
         * quads
         */
        case 6:
        {
            // volume = area[base] * height / 3
            //        = weight
            // area = 0.5 * | p x q | , where p and q are the diagonals
            //     note: z coordinate is 0.0
            //     ref: http://mathworld.wolfram.com/Quadrilateral.html
            // Vector2 p = *renderList[3] - *renderList[1];
            // Vector2 q = *renderList[4] - *renderList[2];
            // float abs_cross = abs( p.x * q.y - p.y * q.x );
            // area = 0.5f * abs_cross;

            // do it in one fast computation:
            area = 0.5 * NPWGeometryMath::Abs(
                              ( (double)(vertices->renderList[3]->x - vertices->renderList[1]->x) *
                                (double)(vertices->renderList[4]->y - vertices->renderList[2]->y) ) -
                              ( (double)(vertices->renderList[3]->y - vertices->renderList[1]->y) *
                                (double)(vertices->renderList[4]->x - vertices->renderList[2]->x) ) );
            height = (float)( 3.0 * weight / NPWGeometryMath::Max(area, MINIMUM_AREA) );
            break;
        }

        default:
        {
            height = 0.0f;
            break;
        }
    }
}


/**
 *  Renders the projection of the tetrahedron
 */
void NPWOpenGLRenderer::render()
{
    switch ( vertices->n_renderList_size )
    {

        /*
         * points
         */
        case 1:
        {
            /*
             * note: no offset, since we are rendering directly to the histogram
             */
            int x = NPWGeometryMath::Floor( vertices->renderList[0]->x );
            int y = NPWGeometryMath::Floor( vertices->renderList[0]->y );
            x = NPWGeometryMath::Max( 0, NPWGeometryMath::Min( histogram->dim(0) - 1, x ) );
            y = NPWGeometryMath::Max( 0, NPWGeometryMath::Min( histogram->dim(1) - 1, y ) );
            (*histogram)(x,y) += weight;
            break;
        }

        /*
         * triangles with a vertex pair, with unsafe borders
         */
        case 3:
        {
            glBegin( GL_TRIANGLE_FAN );

                // only the first vertex has non-zero height
                glVertex3f( vertices->renderList[0]->x + offset,
                            vertices->renderList[0]->y + offset, height );

                // all others have height of zero
                glVertex3f( vertices->renderList[1]->x + offset,
                            vertices->renderList[1]->y + offset, 0.0f );
                glVertex3f( vertices->renderList[2]->x + offset,
                            vertices->renderList[2]->y + offset, 0.0f );

            glEnd();
            break;
        }

        /*
         * triangles with a vertex on the border, with unsafe borders
         */
        case 4:
        {
            glBegin( GL_TRIANGLE_FAN );

                // only the first vertex has non-zero height
                glVertex3f( vertices->renderList[0]->x + offset,
                            vertices->renderList[0]->y + offset, height );

                // all others have height of zero
                glVertex3f( vertices->renderList[1]->x + offset,
                            vertices->renderList[1]->y + offset, 0.0f );
                glVertex3f( vertices->renderList[2]->x + offset,
                            vertices->renderList[2]->y + offset, 0.0f );
                glVertex3f( vertices->renderList[3]->x + offset,
                            vertices->renderList[3]->y + offset, 0.0f );

            glEnd();
            break;
        }

        /*
         * triangles
         */
        case 5:
        {
            glBegin( GL_TRIANGLE_FAN );

                // only the first vertex has non-zero height
                glVertex3f( vertices->renderList[0]->x + offset,
                            vertices->renderList[0]->y + offset, height );

                // all others have height of zero
                glVertex3f( vertices->renderList[1]->x + offset,
                            vertices->renderList[1]->y + offset, 0.0f );
                glVertex3f( vertices->renderList[2]->x + offset,
                            vertices->renderList[2]->y + offset, 0.0f );
                glVertex3f( vertices->renderList[3]->x + offset,
                            vertices->renderList[3]->y + offset, 0.0f );
                glVertex3f( vertices->renderList[4]->x + offset,
                            vertices->renderList[4]->y + offset, 0.0f );

            glEnd();
            break;
        }

        /*
         * quads
         */
        case 6:
        {
            glBegin( GL_TRIANGLE_FAN );

                // only the first vertex has non-zero height
                glVertex3f( vertices->renderList[0]->x + offset,
                            vertices->renderList[0]->y + offset, height );

                // all others have height of zero
                glVertex3f( vertices->renderList[1]->x + offset,
                            vertices->renderList[1]->y + offset, 0.0f );
                glVertex3f( vertices->renderList[2]->x + offset,
                            vertices->renderList[2]->y + offset, 0.0f );
                glVertex3f( vertices->renderList[3]->x + offset,
                            vertices->renderList[3]->y + offset, 0.0f );
                glVertex3f( vertices->renderList[4]->x + offset,
                            vertices->renderList[4]->y + offset, 0.0f );
                glVertex3f( vertices->renderList[5]->x + offset,
                            vertices->renderList[5]->y + offset, 0.0f );

            glEnd();
            break;
        }

        default:
        {
            break;
        }
    }
}



void NPWOpenGLRenderer::renderTetrahedron(
		const Vector2 & v1, const Vector2 & v2, const Vector2 & v3, const Vector2 & v4,
		double weight_n )
{
	/*
	 * set the vertices and weight
	 */
	weight = weight_n;

	vertices->setVerticesAndDetectShape( v1,v2,v3,v4 );

	/*
	 * set the render list
	 */
	vertices->setRenderList();

	/*
     * calculate height
     */
    calculateHeight();

    /*
     * render the shape
     */
	render();

}

void NPWOpenGLRenderer::setHistogram( NPWGLImage * histogram_n )
{
    histogram = histogram_n;
}

#ifdef DB_COUNT
void NPWOpenGLRenderer::printCounters()
{
    if ( vertices )
    {
        vertices->printCounters();
    }
}
#endif


void NPWOpenGLRenderer::turnShapeExpansionOn()
{
    vertices->turnShapeExpansionOn();
}


void NPWOpenGLRenderer::turnShapeExpansionOff()
{
    vertices->turnShapeExpansionOff();
}

bool NPWOpenGLRenderer::isShapeExpansionOn()
{
    return vertices->isShapeExpansionOn();
}

