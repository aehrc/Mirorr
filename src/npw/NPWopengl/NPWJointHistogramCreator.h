/*
 * NPWJointHistogramCreator.h
 *
 *  Created on: 14/05/2010
 *      Author: bro86j
 */


#ifndef NPW_JOINT_HISTOGRAM_CREATOR_H
#define NPW_JOINT_HISTOGRAM_CREATOR_H

#include "NPWOpenGLRenderer.h"
#include "NPWGLImage.h"
#include "NPWWindow.h"
#include "Vector2.h"

class NPWJointHistogramCreator
{
public:

	/**
	 * constructor
	 */
	NPWJointHistogramCreator();

	/**
	 * destructor
	 */
	virtual ~NPWJointHistogramCreator();

	/**
	 * fills a joint histogram of image1 and image2 using NPW on OpenGL
	 */
	int getHistogram(
	        const NPWGLImage * baseImage, const NPWGLImage * searchImage, NPWGLImage * hist );

   /**
	 * When registering images by splitting it up in small pieces (ie blockmatching),
	 * use the same bin size for all blocks.
	 * Therefore these needs to be set by the blockmatcher
	 */
	void setBaseBinSize  ( const float min, const float max, const int bins );
	void setSearchBinSize( const float min, const float max, const int bins );

private:

	/**
	 * The window
	 */
	NPWWindow * window;

	/**
	 * Creates a hidden window which allows OpenGL state. Uses SDL
	 */
	void createWindow();

	/**
	 * assumes the dimensions of output are scaleFactor smaller than input
	 * scales input down and updates output
	 */
	void scaleDown( const NPWGLImage * input, NPWGLImage * output );

	/**
	 * normalizes the data, this is necessary if the render error is large
	 */
	void normalizeData( NPWGLImage * data );

	/**
	 * When registering images by splitting it up in small pieces (ie. blockmatching),
	 * use the same bin size for all blocks.
	 * Therefore these needs to be set by the blockmatcher
	 */
	float baseBinSize;
	float searchBinSize;
	int   baseBins;
	int   searchBins;
	float baseMin;
	float searchMin;

	/**
	 * Vertices may be moved around a bit, this decreases accuracy. To correct this, the histogram
	 * may be rendered onto a larger frame buffer which is subsequently scaled down.
	 * This factor specifies how many times larger the frame buffer is in comparison to the histogram
	 * It is clamped between 1 and the maximum vertex correction.
	 */
	int scaleFactor;

	/**
	 * The eight voxels in a neighbourhood
	 */
	npw::Vector2 alpha, beta, gamma, delta, epsilon, dzeta, eta, theta;


public:

	/*
	 * returns the renderer
	 */
	NPWOpenGLRenderer * getRenderer() { return window->getRenderer(); }

	/*
	 * turn shape expansion on or off
	 */
	/**
	 * if shape expansion is turned on, every point of maximum intensity in every shape
	 * is guaranteed to be a certain distance from all shape edges, minimizing
	 * rasterizing errors.
	 *
	 * if this is turned off, more errors are introduced and the histogram is normalized in
	 * the end to correct this.
	 */
	void turnShapeExpansionOn()   { window->getRenderer()->turnShapeExpansionOn(); }
	void turnShapeExpansionOff()  { window->getRenderer()->turnShapeExpansionOff(); }
	bool isShapeExpansionOn()     { return window->getRenderer()->isShapeExpansionOn(); }

	/*
	 * set the scale factor, this is clipped to the maximum vertex correction.
	 * scaling down further does not increase accuracy
	 */
	void setScaleFactor( int s );
	int  getScaleFactor()   { return scaleFactor; };

};

#endif // NPW_JOINT_HISTOGRAM_CREATOR_H
