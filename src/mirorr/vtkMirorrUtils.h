/*=========================================================================
Program: mirorr
Module: vtkMirorrUtils.h
Author: David Rivest-Henault
Created: 12 Dec 2012

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*//

#ifndef VTKMIRORRUTILS_H_
#define VTKMIRORRUTILS_H_


#include <vector>
#include <string>
#include "itkPoint.h"


namespace vtk
{

class vtkMirorrUtils
{
public:

  typedef std::vector<itk::Point<double,3> > PointListType;

  /**
   * Save the vector defined by {originPt, endPt - originPt} in
   * a VTK unstructured grid file.
   *
   * Recommended file extension is .vtu
   */
  static bool saveVtkVectorFile(std::string filename,
      const PointListType &originPt, const PointListType &endPt);

private:
  vtkMirorrUtils() {};
  virtual
  ~vtkMirorrUtils() {};
};

} /* namespace milx */
#endif /* VTKMIRORRUTILS_H_ */
