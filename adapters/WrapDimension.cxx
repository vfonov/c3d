#include "WrapDimension.h"
#include "itkWrapPadImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "vnl/vnl_finite.h"

template <class TPixel, unsigned int VDim>
void
WrapDimension<TPixel, VDim>
::operator() (const IndexType &xWrap)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // We will wrap the image in place. We must first calculate the offset
  // in the flattened image that corresponds to the wrap
  int k = 0, n = 1;
  for(size_t i = 0; i < VDim; i++)
    {
    k += xWrap[i] * n;
    n *= img->GetBufferedRegion().GetSize(i);
    }

  // Compute the GCD of number of voxels, offset
  int gcd = vnl_finite_int<0>::gcd(abs(k), n);
  int m = n / gcd;
 
  // Explain what we are doing
  *c->verbose << "Wrapping image #" << c->m_ImageStack.size() << " by " << xWrap << endl;
  *c->verbose << "  Rotate in memory by " << k << " bytes." << endl;
  *c->verbose << "  GCD(" << n << "," << k << ") = " << gcd << endl;

  // Get data pointer
  TPixel *ptr = img->GetBufferPointer();

  // There are gcd different starting points for the wrapping
  for(int i = 0; i < gcd; i++)
    {
    TPixel q = ptr[i];
    int p = i;
    for(int j = 1; j < m; j++)
      {
      int pnext = (p + k) % n;
      ptr[p] = ptr[pnext];
      p = pnext;
      if(p < 0)
        p = p + n;
      }
    ptr[p] = q;
    }

  // Change the origin. Wrapping by (dx,dy,dz) means that voxel (0,0,0) is changed
  // to the value of voxel (dx,dy,dz). So, the coordinate of (dx,dy,dz) should be
  // the new origin.
  itk::Point<double, VDim> org;
  img->TransformIndexToPhysicalPoint(xWrap, org);
  img->SetOrigin(org);

  // Flag our image as modified
  img->Modified();
}

// Invocations
template class WrapDimension<double, 2>;
template class WrapDimension<double, 3>;
