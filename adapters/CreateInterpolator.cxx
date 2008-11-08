#include "CreateInterpolator.h"
#include "itkInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"

template <class TPixel, unsigned int VDim>
typename CreateInterpolator<TPixel, VDim>::InterpolatorType *
CreateInterpolator<TPixel, VDim>
::operator() ()
{
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double> NNInterpolatorType;
  typedef itk::LinearInterpolateImageFunction<ImageType,double> LinearInterpolatorType;
  typedef itk::BSplineInterpolateImageFunction<ImageType,double> CubicInterpolatorType;

  // Create an interpolating function
  if(c->m_Interpolation == "NearestNeighbor")
    m_Interp = NNInterpolatorType::New();
  else if(c->m_Interpolation == "Linear")
    m_Interp = LinearInterpolatorType::New();
  else if(c->m_Interpolation == "Cubic")
    m_Interp = CubicInterpolatorType::New();
  else throw string("Unknown interpolator type");

  return m_Interp;
}

// Invocations
template class CreateInterpolator<double, 2>;
template class CreateInterpolator<double, 3>;
