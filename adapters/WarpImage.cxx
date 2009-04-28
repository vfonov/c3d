#include "WarpImage.h"
#include "itkWarpImageFilter.h"
#include "itkImageToVectorImageFilter.h"
#include "CreateInterpolator.h"

template <class TPixel, unsigned int VDim>
void
WarpImage<TPixel, VDim>
::operator() ()
{
  // Check input availability
  if(c->m_ImageStack.size() < VDim + 1)
    {
    cerr << "Warp operation requires " << VDim+1 << " images on the stack" << endl;
    throw -1;
    }

  // Write something 
  *c->verbose << "Warping image #" << c->m_ImageStack.size() << endl;

  // Get the image to warp
  ImagePointer isrc = c->m_ImageStack[c->m_ImageStack.size() - 1];

  // Store the direction temporarily
  // itk::Matrix<double, VDim, VDim> dirin = isrc->GetDirection(), dirtemp;
  // dirtemp.SetIdentity();
  // isrc->SetDirection(dirtemp);

  // Create a deformation field
  typedef itk::ImageToVectorImageFilter<ImageType> VectorFilter;
  typename VectorFilter::Pointer fltVec = VectorFilter::New();

  // Get the array of warps
  for(size_t d = 0; d < VDim; d++)
    {
    size_t k = c->m_ImageStack.size() + d - (VDim + 1); 
    ImagePointer iwarp = c->m_ImageStack[k];
    fltVec->SetInput(d, iwarp);
    }

  // Compute the field
  fltVec->Update();
  typedef typename VectorFilter::OutputImageType FieldType;
  typename FieldType::Pointer field = fltVec->GetOutput();

  // Create the warp filter
  typedef itk::WarpImageFilter<ImageType, ImageType, FieldType> WarpType;
  typename WarpType::Pointer fltWarp = WarpType::New();
  fltWarp->SetInput(isrc);
  fltWarp->SetDeformationField(field);

  // Create interpolator
  fltWarp->SetInterpolator(c->GetInterpolator());

  // Update the warp fileter
  fltWarp->SetOutputSpacing(field->GetSpacing());
  fltWarp->SetOutputOrigin(field->GetOrigin());
  fltWarp->SetOutputDirection(field->GetDirection());
  fltWarp->SetEdgePaddingValue(c->m_Background);
  fltWarp->Update();

  // Update the output image
  //
  ImagePointer imgout = fltWarp->GetOutput();

  // Drop image and warps from stack
  c->m_ImageStack.pop_back();
  for(size_t d = 0; d < VDim; d++)
    c->m_ImageStack.pop_back();

  // Add the warped image to the stack
  c->m_ImageStack.push_back(imgout);
}

// Invocations
template class WarpImage<double, 2>;
template class WarpImage<double, 3>;
