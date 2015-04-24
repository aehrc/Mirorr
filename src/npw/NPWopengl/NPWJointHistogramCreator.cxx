/*
 * NPWJointHistogramCreator.cxx
 *
 *  Created on: 14/05/2010
 *      Author: bro86j
 */

#include "NPWJointHistogramCreator.h"
#include "NPWStopwatch.h"
#include "NPWGeometry.h"

#include <iostream>
#include <stdexcept>

//#define DB_STATUS
//#define DB_SAFE_BUT_SLOW


NPWJointHistogramCreator::NPWJointHistogramCreator() :
    baseBinSize(-1.0f), searchBinSize(-1.0f),
    baseBins(-1), searchBins(-1),
    baseMin(-1.0f), searchMin(-1.0f),
    scaleFactor ( 2 )
{
#ifdef DB_STATUS
    std::cout << "constructing NPWJointHistogramCreator object"
              << "    status messages ON" << std::endl;
#endif

    window = new NPWWindow(800, 600, false);

}

NPWJointHistogramCreator::~NPWJointHistogramCreator()
{
	if (window)
	{
	    window->destroy();
	    delete window;
	    window = 0;
	}

#ifdef DB_STATUS
	std::cout << "destroying NPWJointHistogramCreator object" << std::endl;
#endif
}


int NPWJointHistogramCreator::getHistogram(
		const NPWGLImage * baseImage, const NPWGLImage * searchImage, NPWGLImage * hist )
{
#ifdef DB_SAFE_BUT_SLOW
	if ( !baseImage || !searchImage )
	{
		std::cout << "*** ERROR ***: image null pointer exception" << std::endl;
		return -1;
	}
	if ( !hist )
	{
		std::cout << "*** ERROR ***: histogram null pointer exception" << std::endl;
		return -1;
	}

	if ( baseImage->dim(0) != searchImage->dim(0) ||
		 baseImage->dim(1) != searchImage->dim(1) ||
		 baseImage->dim(2) != searchImage->dim(2) )
	{
		std::cout << "*** ERROR ***: images are of unequal dimensions" << std::endl;
		return -2;
	}

	if ( baseBinSize < 0.0f || searchBinSize < 0.0f )
	{
	    std::cout << "*** ERROR ***: bin sizes are not set (or smaller than 0.0)" << std::endl;
	    return -2;
	}

	if ( baseBins != hist->dim(0) || searchBins != hist->dim(1) )
	{
	    std::cout << "*** ERROR ***: dimensions of supplied histogram do not match "
	              << "those used to calculate bin sizes: "
	              << "first dimension should correspond to base image and "
	              << "second dimension to search image" << std::endl;
	    return -2;
	}
#endif

#ifdef DB_STATUS
	NPWStopwatch stopwatch;
#endif

	NPWGLImage * renderHistogram = hist;

	/*
	 * if we use scaling, we need to create a new histogram to render on
	 */
	if ( scaleFactor != 1 )
	{
	    renderHistogram = new NPWGLImage( baseBins, searchBins );
	}

	window->initialize( renderHistogram );

#ifdef DB_STATUS
	std::cout << "Rendering joint histogram...\t ( "
	          << "shape expansion " << (window->getRenderer()->isShapeExpansionOn()?"ON":"OFF")
	          << ",\tscale factor " << scaleFactor << ",\t"
	          << baseBins/scaleFactor << "x" << searchBins/scaleFactor << " bins )"
	          << std::endl;
#endif

	int maxX = baseImage->dim(0) - 1;
	int maxY = baseImage->dim(1) - 1;
	int maxZ = baseImage->dim(2) - 1;

	double outerWeight = 1.0 / static_cast<double>( 6 * baseImage->GetNeighborhoods() );
	double innerWeight = 1.0 / static_cast<double>( 3 * baseImage->GetNeighborhoods() );

	for ( int z = 0; z < maxZ; ++z )
	{
		for ( int y = 0; y < maxY; ++y )
		{
			for ( int x = 0; x < maxX; ++x )
			{
				// get tetrahedron vertices in bin coordinates
				alpha.x   = ( (*baseImage)(x  ,y  ,z  ) - baseMin ) / baseBinSize;
				beta.x    = ( (*baseImage)(x+1,y  ,z  ) - baseMin ) / baseBinSize;
				gamma.x   = ( (*baseImage)(x  ,y+1,z  ) - baseMin ) / baseBinSize;
				delta.x   = ( (*baseImage)(x+1,y+1,z  ) - baseMin ) / baseBinSize;
				epsilon.x = ( (*baseImage)(x  ,y  ,z+1) - baseMin ) / baseBinSize;
				dzeta.x   = ( (*baseImage)(x+1,y  ,z+1) - baseMin ) / baseBinSize;
				eta.x     = ( (*baseImage)(x  ,y+1,z+1) - baseMin ) / baseBinSize;
				theta.x   = ( (*baseImage)(x+1,y+1,z+1) - baseMin ) / baseBinSize;

				alpha.y   = ( (*searchImage)(x  ,y  ,z  ) - searchMin ) / searchBinSize;
				beta.y    = ( (*searchImage)(x+1,y  ,z  ) - searchMin ) / searchBinSize;
				gamma.y   = ( (*searchImage)(x  ,y+1,z  ) - searchMin ) / searchBinSize;
				delta.y   = ( (*searchImage)(x+1,y+1,z  ) - searchMin ) / searchBinSize;
				epsilon.y = ( (*searchImage)(x  ,y  ,z+1) - searchMin ) / searchBinSize;
				dzeta.y   = ( (*searchImage)(x+1,y  ,z+1) - searchMin ) / searchBinSize;
				eta.y     = ( (*searchImage)(x  ,y+1,z+1) - searchMin ) / searchBinSize;
				theta.y   = ( (*searchImage)(x+1,y+1,z+1) - searchMin ) / searchBinSize;

				/*
				 * outer tetrahedra
				 */
                window->getRenderer()->renderTetrahedron(alpha, beta   , gamma, epsilon, outerWeight);
                window->getRenderer()->renderTetrahedron(beta , gamma  , delta, theta  , outerWeight);
                window->getRenderer()->renderTetrahedron(beta , epsilon, dzeta, theta  , outerWeight);
                window->getRenderer()->renderTetrahedron(gamma, epsilon, eta  , theta  , outerWeight);

                /*
                 * inner tetrahedron
                 */
                window->getRenderer()->renderTetrahedron(beta, gamma, epsilon, dzeta, innerWeight);

			}
		}
	}


#ifdef DB_TIME
//	window->getRenderer()->printCounters();
#endif
#ifdef DB_STATUS
	std::cout << "reading back texture data (lap time is " << stopwatch.lap() << " s)"
	          << std::endl;
#endif

	if (!window->readBackImage(renderHistogram))
	{
		std::cout << "*** ERROR ***: could not read back texture data" << std::endl;
		return -3;
	}

	/*
	 * if we use scaling, we have to scale the rendered histogram down
	 * and normalize it
	 */
	if ( scaleFactor != 1 )
	{
	    scaleDown( renderHistogram, hist );

	    delete renderHistogram;
	}

	if ( !window->getRenderer()->isShapeExpansionOn() )
	{
	    normalizeData( hist );
	}

#ifdef DB_STATUS
	std::cout << "finished rendering joint histogram in "
			  << stopwatch.lap() << " s" << std::endl;
#endif

	return 1;
}

void NPWJointHistogramCreator::setBaseBinSize(
        const float min, const float max, const int bins )
{
    if ( min == max )
    {
        // at least one image contains only one intensity, cannot render histogram
        std::cout << "*** ERROR *** cannot set bin sizes because base image is mono-valued"
                  << std::endl;
        throw std::runtime_error("Quit because base bin size could not be set");
    }

    baseBins = bins * scaleFactor;
    baseMin = min;

    baseBinSize = (max - baseMin) / static_cast<float>( baseBins );
}

void NPWJointHistogramCreator::setSearchBinSize(
        const float min, const float max, const int bins )
{
    if ( min == max )
    {
        // at least one image contains only one intensity, cannot render histogram
        std::cout << "*** ERROR *** cannot set bin sizes because search image is mono-valued"
                  << std::endl;
        throw std::runtime_error("Quit because search bin size could not be set");
    }

    searchBins = bins * scaleFactor;
    searchMin  = min;

    searchBinSize = (max - searchMin) / static_cast<float>( searchBins );
}

void NPWJointHistogramCreator::setScaleFactor( int s )
{
    int oldScaleFactor = scaleFactor;
    scaleFactor = NPWGeometryMath::Min(
                       NPWGeometryMath::Max(s, 1 ),
                       NPWGeometryMath::Ceil(npw::MAXIMUM_VERTEX_CORRECTION) );

    searchBins    *= scaleFactor / oldScaleFactor;
    baseBins      *= scaleFactor / oldScaleFactor;

    searchBinSize *= (float)( oldScaleFactor / scaleFactor );
    baseBinSize   *= (float)( oldScaleFactor / scaleFactor );

}


void NPWJointHistogramCreator::scaleDown( const NPWGLImage * input, NPWGLImage * output )
{
    if ( scaleFactor == 1 )                                 { return; }

#ifdef DB_SAFE_BUT_SLOW
    if ( !input || !output )                                { return; }

    if ( output->dim(0) * scaleFactor != input->dim(0) ||
         output->dim(1) * scaleFactor != input->dim(1) )    { return; }
#endif

    float nTotal;

    for ( int y = 0; y < output->dim(0); ++y )
    {
        for ( int x = 0; x < output->dim(1); ++x )
        {
            /*
             * get the total of the neighborhoods
             */
            nTotal = 0.0f;
            for ( int yn = 0; yn < scaleFactor; ++yn )
            {
                for (int xn = 0; xn < scaleFactor; ++ xn )
                {
                    nTotal += (*input)( x*scaleFactor+xn, y*scaleFactor+yn );
                }
            }

            /*
             * add the total to the output
             * don't average it, because the
             * input and output should both sum to one
             */
            (*output)(x,y) = nTotal;
        }
    }

}


void NPWJointHistogramCreator::normalizeData( NPWGLImage * data )
{
    /*
     * get the total
     */
    float total = 0.0f;

    for ( int y = 0; y < data->dim(1); ++y )
    {
        for ( int x = 0; x < data->dim(0); ++x )
        {
            total += (*data)(x,y);
        }
    }

    float invt = 1.0f / total;

    /*
     * normalize data
     */
    for ( int y = 0; y < data->dim(1); ++y )
    {
        for ( int x = 0; x < data->dim(0); ++x )
        {
            (*data)(x,y) *= invt;
        }
    }

}

