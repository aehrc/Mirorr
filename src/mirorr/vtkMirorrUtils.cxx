/*=========================================================================
Program: mirorr
Module: vtkMirorrUtils.cxx
Author: David Rivest-Henault
Created: 12 Dec 2012

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include "vtkMirorrUtils.h"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkLine.h"
#include "vtkCellType.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"


namespace vtk
{

bool vtkMirorrUtils::saveVtkVectorFile(std::string filename,
      const PointListType &originPt, const PointListType &endPt)
{
  const size_t N = originPt.size();
  if (endPt.size() != N) {
    return false;
  }

  vtkSmartPointer<vtkUnstructuredGrid> grid =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkPoints> points =
      vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkDoubleArray> scalars =
        vtkSmartPointer<vtkDoubleArray>::New();
  scalars->SetNumberOfComponents(1);

  vtkSmartPointer<vtkDoubleArray> vects =
      vtkSmartPointer<vtkDoubleArray>::New();
  vects->SetNumberOfComponents(3);

  for (size_t i=0; i<N; i++) {
    points->InsertPoint(i, originPt[i].GetDataPointer());
    //points->InsertPoint(i, endPt[i].GetDataPointer());

    itk::Vector<double,3> vec = endPt[i] - originPt[i];
    vects->InsertTuple(i, vec.GetDataPointer());

    double v = vec.GetNorm();
    scalars->InsertTuple(i, &v);
  }

  grid->SetPoints(points);

  scalars->SetName("BMVectorMagnitude");
  grid->GetPointData()->SetScalars(scalars);

  vects->SetName("BMVector");
  grid->GetPointData()->SetVectors(vects);

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  writer->SetFileName(filename.c_str());
  writer->SetInput(grid);
  writer->Update();

  return true;
}

} /* namespace milx */
