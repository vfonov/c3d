#ifndef __ApplyMetric_h_
#define __ApplyMetric_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ApplyMetric : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ApplyMetric(Converter *c) : c(c) {}

  void operator() (const char *metric);

private:
  Converter *c;

};

#endif

