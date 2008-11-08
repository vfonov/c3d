#ifndef __TrimImage_h_
#define __TrimImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class TrimImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  TrimImage(Converter *c) : c(c) {}

  void operator() (const RealVector &margin);

private:
  Converter *c;

};

#endif

