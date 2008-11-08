#include "ReciprocalImage.h"

template <class TPixel, unsigned int VDim>
void
ReciprocalImage<TPixel, VDim>
::operator() ()
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Verbose
  *c->verbose << "Taking the reciprocal of #" << c->m_ImageStack.size() << endl;

  // Simply go through and divide
  size_t n = img->GetBufferedRegion().GetNumberOfPixels();
  TPixel *buffer = img->GetBufferPointer();
  for(size_t i = 0; i < n; i++)
    buffer[i] = 1.0 / buffer[i];
  img->Modified();
}

// Invocations
template class ReciprocalImage<double, 2>;
template class ReciprocalImage<double, 3>;
