#include "ApplyMetric.h"
#include "itkMutualInformationHistogramImageToImageMetric.h"
#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkAffineTransform.h"


template <class TPixel, unsigned int VDim>
void
ApplyMetric<TPixel, VDim>
::operator() (const char *metric_name)
{
  // Two images must be on a stack
  if(c->m_ImageStack.size() < 2)
    throw ConvertException("Two images required for metric computation");

  // Get the images
  ImagePointer iref = c->m_ImageStack[c->m_ImageStack.size() - 2];
  ImagePointer imov = c->m_ImageStack[c->m_ImageStack.size() - 1];

  // Create the appropriate metric
  typedef itk::ImageToImageMetric<ImageType, ImageType> MetricType;
  typename MetricType::Pointer metric;
  if(!strcmp(metric_name,"MI"))
    {
    metric = 
      itk::MutualInformationHistogramImageToImageMetric<
        ImageType,ImageType>::New();
    }
  else if(!strcmp(metric_name,"NMI"))
    {
    metric = 
      itk::NormalizedMutualInformationHistogramImageToImageMetric<
        ImageType,ImageType>::New();
    }
  else throw ConvertException("Unknown metric %s", metric_name);

  // Configure the identity transform
  typedef itk::AffineTransform<double, VDim> TransformType;
  typename TransformType::Pointer tran = TransformType::New();
  tran->SetIdentity();

  // Interpolator
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double> NNInterpolatorType;
  metric->SetInterpolator(NNInterpolatorType::New());

  // Configure the metric
  metric->SetMovingImage(imov);
  metric->SetFixedImage(iref);
  metric->SetTransform(tran);
  metric->SetFixedImageRegion(iref->GetBufferedRegion());
  metric->Initialize();

  double mvalue = metric->GetValue(tran->GetParameters());

  // Print the value
  cout << metric_name << " = " << mvalue << endl;

  // Get image from stack
  // ImagePointer img = c->m_ImageStack.back();

  // Do some processing ...
  // ImagePointer result = ...;
  
  // Put result on stack
  // c->m_ImageStack.pop_back();
  // c->m_ImageStack.push_back(result);
}

// Invocations
template class ApplyMetric<double, 2>;
template class ApplyMetric<double, 3>;
