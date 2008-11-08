#include "LabelStatistics.h"
#include "itkLabelStatisticsImageFilter.h"
#include <set>

template <class TPixel, unsigned int VDim>
void
LabelStatistics<TPixel, VDim>
::operator() ()
{
  // Get images from stack
  size_t n = c->m_ImageStack.size();
  if(n < 2)
    throw ConvertException("Label statistics require two image inputs");
  ImagePointer label = c->m_ImageStack[n-1];
  ImagePointer image = c->m_ImageStack[n-2];

  // Create a short image for the labels
  typedef itk::Image<short, VDim> LabelImageType;
  typename LabelImageType::Pointer slab = LabelImageType::New();

  // Allocate the image
  slab->SetRegions(label->GetBufferedRegion());
  slab->Allocate();

  // Accumulator of label values    
  std::set<short> sval;

  // Round off doubles to create labels
  size_t nv = label->GetBufferedRegion().GetNumberOfPixels();
  for(size_t i = 0; i < nv; i++)
    {
    slab->GetBufferPointer()[i] = (short) (label->GetBufferPointer()[i] + 0.5);
    sval.insert(slab->GetBufferPointer()[i]);
    }

  // Create the label statistics filter
  typedef itk::LabelStatisticsImageFilter<ImageType, LabelImageType> StatFilter;
  typename StatFilter::Pointer fltStat = StatFilter::New();
  
  // Set the inputs
  fltStat->SetInput(image);
  fltStat->SetLabelInput(slab);

  // Update the filter
  fltStat->Update();

  // Get the number of labels 
  printf("LabelID        Mean        StdD         Max         Min       Count\n");
  for(set<short>::iterator it = sval.begin(); it != sval.end(); ++it)
    {
    // printf("xxxxx    xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx");
    printf("%5i    %10.2f  %10.2f  %10.2f  %10.2f  %10lu\n",
      (int) *it, 
      fltStat->GetMean(*it), 
      fltStat->GetSigma(*it), 
      fltStat->GetMaximum(*it), 
      fltStat->GetMinimum(*it), 
      (long unsigned) fltStat->GetCount(*it));
    }
}

// Invocations
template class LabelStatistics<double, 2>;
template class LabelStatistics<double, 3>;

