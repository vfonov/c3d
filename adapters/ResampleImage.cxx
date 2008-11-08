#include "ResampleImage.h"
#include "itkResampleImageFilter.h"

template <class TPixel, unsigned int VDim>
void
ResampleImage<TPixel, VDim>
::operator() (SizeType &sz)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();
  
  // Build the resampling filter
  typedef itk::ResampleImageFilter<ImageType,ImageType> ResampleFilterType;
  typename ResampleFilterType::Pointer fltSample = ResampleFilterType::New();

  // Initialize the resampling filter with an identity transform
  fltSample->SetInput(input);
  fltSample->SetTransform(itk::IdentityTransform<double,VDim>::New());
  fltSample->SetInterpolator(c->GetInterpolator());

  // Compute the spacing of the new image
  typename ImageType::SpacingType spc_pre = input->GetSpacing();
  typename ImageType::SpacingType spc_post = spc_pre;
  for(size_t i = 0; i < VDim; i++)
    spc_post[i] *= input->GetBufferedRegion().GetSize()[i] * 1.0 / sz[i];

  // Get the bounding box of the input image
  typename ImageType::PointType origin_pre = input->GetOrigin();

  // Recalculate the origin. The origin describes the center of voxel 0,0,0
  // so that as the voxel size changes, the origin will change as well.
  typename ImageType::SpacingType off_pre = (input->GetDirection() * spc_pre) * 0.5;
  typename ImageType::SpacingType off_post = (input->GetDirection() * spc_post) * 0.5;
  typename ImageType::PointType origin_post = origin_pre - off_pre + off_post;
  
  // Set the image sizes and spacing.
  fltSample->SetSize(sz);
  fltSample->SetOutputSpacing(spc_post);
  fltSample->SetOutputOrigin(origin_post);
  fltSample->SetOutputDirection(input->GetDirection());

  // Set the unknown intensity to positive value
  fltSample->SetDefaultPixelValue(c->m_Background);

  // Describe what we are doing
  *c->verbose << "Resampling #" << c->m_ImageStack.size() << " to have" << sz << " voxels." << endl;
  *c->verbose << "  Interpolation method: " << c->m_Interpolation << endl;
  *c->verbose << "  Background intensity: " << c->m_Background << endl;
  *c->verbose << "  Input spacing: " << spc_pre << endl;
  *c->verbose << "  Input origin: " << origin_pre << endl;
  *c->verbose << "  Output spacing: " << spc_post << endl;
  *c->verbose << "  Output origin: " << origin_post << endl;

  // Perform resampling
  fltSample->UpdateLargestPossibleRegion();
    
  // Change the source to the output 
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(fltSample->GetOutput());
}

// Invocations
template class ResampleImage<double, 2>;
template class ResampleImage<double, 3>;
