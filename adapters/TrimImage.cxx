#include "TrimImage.h"
#include "ExtractRegion.h"

template<unsigned int VDim>
void ExpandRegion(itk::ImageRegion<VDim> &region, const itk::Index<VDim> &idx)
{
  if(region.GetNumberOfPixels() == 0)
    {
    region.SetIndex(idx);
    for(size_t i = 0; i < VDim; i++)
      region.SetSize(i, 1);
    }
  else {
    for(size_t i = 0; i < VDim; i++)
      {
      if(region.GetIndex(i) > idx[i])
        {
        region.SetSize(i, region.GetSize(i) + (region.GetIndex(i) - idx[i]));
        region.SetIndex(i, idx[i]);
        }
      else if(region.GetIndex(i) + (long) region.GetSize(i) <= idx[i]) {
        region.SetSize(i, 1 + idx[i] - region.GetIndex(i));
        }
      }
  }
}

template <class TPixel, unsigned int VDim>
void
TrimImage<TPixel, VDim>
::operator() (const RealVector &margin)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();

  // Debugging info
  *c->verbose << "Trimming #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Wrapping non-background voxels with margin of " << margin << " mm." << endl;

  // Initialize the bounding box
  RegionType bbox;

  // Find the extent of the non-background region of the image
  Iterator it(input, input->GetBufferedRegion());
  for( ; !it.IsAtEnd(); ++it)
    if(it.Value() != c->m_Background)
      ExpandRegion(bbox, it.GetIndex());

  // Pad the region by radius specified by user
  typename ImageType::SizeType radius;
  for(size_t i = 0; i < VDim; i++)
    radius[i] = (int) ceil(margin[i] / input->GetSpacing()[i]);
  bbox.PadByRadius(radius);

  // Use the extract region code for the rest
  ExtractRegion<TPixel, VDim> extract(c);
  extract(bbox);
}

// Invocations
template class TrimImage<double, 2>;
template class TrimImage<double, 3>;
