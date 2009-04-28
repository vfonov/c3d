#ifndef __SimpleElasticRegistration_h_
#define __SimpleElasticRegistration_h_

#include "ConvertAdapter.h"
#include "itkGaussianInterpolateImageFunction.h"

template<class TPixel, unsigned int VDim>
class SimpleElasticRegistration : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  // Other typedefs
  typedef itk::FixedArray<float, VDim>  VectorType;
  typedef itk::Image<VectorType, VDim> VectorImageType;
  typedef typename VectorImageType::Pointer VectorImagePointer;
  typedef typename itk::ImageRegionIteratorWithIndex<VectorImageType> VectorIterator;

  typedef itk::GaussianInterpolateImageFunction<ImageType, double> GaussInterpolator;

  SimpleElasticRegistration(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

  double ComputeObjective();

  // Stuff for optimization
  ImagePointer iref, imov;
  VectorImagePointer v;

  // Interpolators
  typename GaussInterpolator::Pointer giref, gimov;

};

#endif

