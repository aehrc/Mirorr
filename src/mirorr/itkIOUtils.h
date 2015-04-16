//
// Created by riv019 on 7/04/15.
//

#ifndef MIRORRTOP_ITKIOUTILS_H
#define MIRORRTOP_ITKIOUTILS_H

#include <string>
#include <iomanip>
#include <itkArray.h>
#include <itkFixedArray.h>
#include <vnl/vnl_vector.h>

/**
 * This file implement pretty printing for object of type itk::Array.
 *
 * This is necessary (as far as pretty printing is) since ITK purposely (why?)
 * prevent user-control of precision for output stream. See:
 *   "itkArrayOutputSpecialization.cxx"
 *   "itkNumberToString.h"
 *   "itkNumberToString.hxx"
 *
 * Use case:
 * itk::Array< double > vec = this->m_Transform->GetParameters();
 * std::cout << std::scientific;
 * std::cout.precision(2);
 * std::cout << "Nicely printed array: " << PrettyPrint(vec, 10) << std::endl;
 * std::cout << "Nicely printed array: " << PrettyPrint(this->m_Transform->GetParameters(), 10) << std::endl;
 */

namespace itk {

template<typename TValue> unsigned int GetArrayLength(const itk::Array<TValue> &arr) {return arr.size();}
template<typename TValue> unsigned int GetArrayLength(const itk::FixedArray<TValue> &arr) {return arr.Size();}
template<typename TValue> unsigned int GetArrayLength(const vnl_vector<TValue> &arr) {return arr.size();}


/** \class PrettyPrinter
 *  \brief
 */
template<typename TArray>
class PrettyPrinter
{
public:
    PrettyPrinter(const TArray &arr, unsigned int width=0, const std::string delim="; ") {
        this->arr = &arr;
        this->width = width;
        this->delim = delim;
    }

private:
    const TArray * arr;
    unsigned int width;
    std::string delim;

public:
    template<typename TArray_> friend std::ostream & operator<< (std::ostream & os, const PrettyPrinter<TArray_> & pp);
};


template<typename TArray>
PrettyPrinter<TArray> PrettyPrint(const TArray &arr, unsigned int width=0, const std::string delim="; ") {
    return PrettyPrinter<TArray>(arr, width, delim);
}


template<typename TArray>
std::ostream & operator<< (std::ostream & os, const PrettyPrinter<TArray> & pp)
{
    const TArray & arr = *(pp.arr);
    os << "[";
    const unsigned int length = GetArrayLength(arr);
    if ( length >= 1 )
    {
        const unsigned int   last   = length - 1;

        if (pp.width > 0) {
            for (unsigned int i = 0; i < last; ++i) {
                os << std::setw(pp.width) << arr[i] << pp.delim;
            }
            os << std::setw(pp.width) << arr[last];
        } else {
            for (unsigned int i = 0; i < last; ++i) {
                os << arr[i] << pp.delim;
            }
            os << arr[last];
        }
    }
    os << "]";
    return os;
}

}
#endif //MIRORRTOP_ITKIOUTILS_H
