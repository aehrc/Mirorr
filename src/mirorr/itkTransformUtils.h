/*=========================================================================
Program: mirorr
Module: itkTransformUtils.h
Author: David Rivest-Henault
Created: 26 Movember 2012

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#ifndef ITKTRANSFORMUTILS_H_
#define ITKTRANSFORMUTILS_H_

namespace itk {

template <typename TransformTypePointer>
void setCenterToOrigin(TransformTypePointer transform) {
  typedef typename TransformTypePointer::ObjectType TransformType;
  typedef typename TransformType::CenterType PointType;
  typedef typename TransformType::OutputVectorType VectorType;

  PointType C0;
  for (unsigned int i=0; i<3; ++i)
    C0[i]=0;
  
  VectorType offset = transform->GetOffset();
  transform->SetCenter(C0);
  transform->SetTranslation(offset);
}


/**
 * Update the transformation centre from C1 to C2 while preserving
 * the relation y = T(x)
 */
template <typename TransformTypePointer, typename CenterType>
void updateCenter(TransformTypePointer transform, CenterType C2) {
  typedef typename TransformTypePointer::ObjectType TransformType;
  typedef typename TransformType::OutputVectorType VectorType;

  setCenterToOrigin(transform);

  CenterType iC2;
  iC2[0] = -C2[0]; iC2[1] = -C2[1]; iC2[2] = -C2[2];
  transform->SetCenter(iC2);
  VectorType off = transform->GetOffset();

  transform->SetCenter(C2);
  transform->SetTranslation(off);
  transform->GetParameters();
}


/**
 * Given an itk transform T centred around C1, so that y = T_{C1} x,
 * compute the itk transform U centred around C2, so that U_{C2} y = x.
 *
 * The reason this function exists is because the so-called 'centered'
 * transforms add a bit of complexity. In pure matrix notation, we have:
 *
 * y = C1 T C1^-1 x, and
 * C2 U C2^-1 y = x.
 *
 * Thus U is not exactly the inverse of T if C1 != C2.
 *
 * ---
 *
 * Template TransformTypePointer is assumed to be an itk::SmartPointer<>
 * templated over an itk transform class, and CenterType should
 * correspond to TransformType::CenterType.
 */
template <typename TransformTypePointer, typename CenterType>
void invertAndChangeCenter(TransformTypePointer transform, CenterType C2) {
  typedef typename TransformTypePointer::ObjectType TransformType;

  TransformTypePointer temp =
        dynamic_cast<TransformType*>(transform->CreateAnother().GetPointer());
  temp->SetCenter(C2);
  transform->GetInverse(temp);

  transform->SetParameters(temp->GetParameters());
  transform->SetFixedParameters(temp->GetFixedParameters());
}

}

#endif /* ITKTRANSFORMUTILS_H_ */
