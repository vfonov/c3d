#include "WarpLabelImage.h"
#include "itkWarpImageFilter.h"
#include "itkImageToVectorImageFilter.h"
#include "CreateInterpolator.h"
#include "ThresholdImage.h"
#include "SmoothImage.h"
#include <set>

template <class TPixel, unsigned int VDim>
void
WarpLabelImage<TPixel, VDim>
::operator() (RealVector &stdev)
{
  // Check input availability
  if(c->m_ImageStack.size() < VDim + 1)
    {
    cerr << "Warp operation requires " << VDim+1 << " images on the stack" << endl;
    throw -1;
    }

  // Write something 
  *c->verbose << "Warping image label-wise #" << c->m_ImageStack.size() << endl;

  // Get the image to warp
  ImagePointer isrc = c->m_ImageStack[c->m_ImageStack.size() - 1];

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

  // Get the origin and spacing of the warp field
  typename ImageType::PointType warp_origin = fltVec->GetInput(0)->GetOrigin();
  typename ImageType::SpacingType warp_spacing = fltVec->GetInput(0)->GetSpacing();

  // Set the origin to zero for both images (why?)
  typename ImageType::PointType zero_origin;
  zero_origin.Fill(0.0);
  isrc->SetOrigin(zero_origin);
  field->SetOrigin(zero_origin);

  // Create the warp filter
  typedef itk::WarpImageFilter<ImageType, ImageType, FieldType> WarpType;
  typename WarpType::Pointer fltWarp = WarpType::New();
  fltWarp->SetDeformationField(field);

  // Create interpolator
  CreateInterpolator<TPixel, VDim> interp(c);
  fltWarp->SetInterpolator(interp());

  // Update the warp fileter
  fltWarp->SetOutputSpacing(field->GetSpacing());
  fltWarp->SetOutputOrigin(field->GetOrigin());
  fltWarp->SetEdgePaddingValue(c->m_Background);

  // Find all unique labels in the input image
  std::set<TPixel> labels;
  size_t n = isrc->GetBufferedRegion().GetNumberOfPixels();
  for(size_t i = 0; i < n; i++)
    labels.insert(isrc->GetBufferPointer()[i]);

  // Create the output image
  ImagePointer ilabel = ImageType::New();
  ilabel->SetRegions(field->GetBufferedRegion());
  ilabel->SetOrigin(warp_origin);
  ilabel->SetSpacing(warp_spacing);
  ilabel->Allocate();

  // Create the 'score' image
  ImagePointer iscore = ImageType::New();
  iscore->SetRegions(field->GetBufferedRegion());
  iscore->Allocate();
  iscore->FillBuffer(0.0);

  // For each unique label in the input image, threshold that image
  for(typename std::set<TPixel>::iterator it = labels.begin(); it!=labels.end(); ++it)
    {
    // Current label value
    TPixel val = *it;
    
    // Threshold the image (get only current label)
    ThresholdImage<TPixel, VDim> thresh(c);
    thresh(val, val, 1, 0);

    // Smooth the image
    SmoothImage<TPixel, VDim> smooth(c);
    smooth(stdev);

    // Apply the warp to the smoothed image
    fltWarp->SetInput(c->m_ImageStack.back());
    fltWarp->Update();

    // Iterate over voxels, update the score and label images
    ImagePointer iwarped = fltWarp->GetOutput();
    size_t nvox = iwarped->GetBufferedRegion().GetNumberOfPixels();
    TPixel *pscore = iscore->GetBufferPointer(), 
           *pwarped = iwarped->GetBufferPointer(), 
           *plabel = ilabel->GetBufferPointer();
    for(size_t i = 0; i < nvox; i++)
      {
      if(pscore[i] < pwarped[i])
        {
        pscore[i] = pwarped[i];
        plabel[i] = val;
        }
      }

    // Replace everything on the stack
    c->m_ImageStack.pop_back();
    c->m_ImageStack.push_back(isrc);
    }

  // Drop image and warps from stack
  c->m_ImageStack.pop_back();
  for(size_t d = 0; d < VDim; d++)
    c->m_ImageStack.pop_back();

  // Store the result on the stack
  c->m_ImageStack.push_back(ilabel);
}

// Invocations
template class WarpLabelImage<double, 2>;
template class WarpLabelImage<double, 3>;
