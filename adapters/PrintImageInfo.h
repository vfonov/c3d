#ifndef __PrintImageInfo_h_
#define __PrintImageInfo_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class PrintImageInfo : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  PrintImageInfo(Converter *c) : c(c) {}

  void operator() (bool full);

private:
  Converter *c;

};

#endif

