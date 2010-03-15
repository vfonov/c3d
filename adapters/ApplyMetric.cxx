#include <string>
#include <iostream>
#include "ApplyMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkMutualInformationHistogramImageToImageMetric.h"
#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkAffineTransform.h"
#include <itkTransformFileReader.h>
#include <itkTransformFactory.h>
#include "itkResampleImageFilter.h"
#include "gsGSAffine3DTransform.h"


template<unsigned int VDim>
void ReadMatrix(const char *fname, itk::Matrix<double,VDim+1,VDim+1> &mat)
  {
  ifstream fin(fname);
  for(size_t i = 0; i < VDim+1; i++)
    for(size_t j = 0; j < VDim+1; j++)
      if(fin.good())
        {
        fin >> mat[i][j];
        }
      else
        {
        throw ConvertException("Unable to read matrix %s", fname);
        }
  fin.close();
  }


template <class TPixel, unsigned int VDim>
void
ApplyMetric<TPixel, VDim>
::operator() (const char *metric_name, const char *fn_tran)
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
  else if(!strcmp(metric_name,"MMI"))
    {
    metric = 
      itk::MattesMutualInformationImageToImageMetric<
        ImageType,ImageType>::New();
    }
  else if(!strcmp(metric_name,"MSQ"))
    {
    metric = 
      itk::MeanSquaresImageToImageMetric<
        ImageType,ImageType>::New();
    }
  else if(!strcmp(metric_name,"NCOR"))
    {
    metric = 
      itk::NormalizedCorrelationImageToImageMetric<
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

  if (!strcmp(fn_tran,"none"))
    tran->SetIdentity();
  else
    {
    // Read the matrix
    itk::Matrix<double,VDim+1,VDim+1> matrix;
    itk::Matrix<double,VDim,VDim> amat;
    itk::Vector<double, VDim> aoff;

    ReadMatrix<VDim>(fn_tran, matrix);
    amat.GetVnlMatrix().update(
      matrix.GetVnlMatrix().extract(VDim, VDim));
    aoff.GetVnlVector().update(
      matrix.GetVnlMatrix().get_column(VDim).extract(VDim));

    // Extrernal matrices are assumed to be RAS to RAS, so we must convert to LPS to LPS
    vnl_vector<double> v_lps_to_ras(VDim, 1.0);
    v_lps_to_ras[0] = v_lps_to_ras[1] = -1.0;
    vnl_diag_matrix<double> m_lps_to_ras(v_lps_to_ras);
    vnl_matrix<double> mold = amat.GetVnlMatrix();
    amat.GetVnlMatrix().update(m_lps_to_ras * mold * m_lps_to_ras);
    aoff.GetVnlVector().update(m_lps_to_ras * aoff.GetVnlVector());

    // Put the values in the transform
    tran->SetMatrix(amat);
    tran->SetOffset(aoff);
    }

  // Interpolator
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double> NNInterpolatorType;
  //metric->DebugOn();
  metric->SetInterpolator(NNInterpolatorType::New());

  // Configure the metric
  metric->SetMovingImage(imov);
  metric->SetFixedImage(iref);
  metric->SetTransform(tran);
  metric->SetFixedImageRegion(iref->GetBufferedRegion());
  metric->Initialize();

  //std::cout << "metric transform parameters: " << metric->GetTransform()->GetParameters() << endl;
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
