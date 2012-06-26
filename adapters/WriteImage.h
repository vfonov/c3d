#ifndef __WriteImage_h_
#define __WriteImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class WriteImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  WriteImage(Converter *c) : c(c) {}

  void operator() (const char *file, bool force, int pos=-1,const char *history=NULL);

  void WriteMultiComponent(const char *file, int ncomp,const char *history=NULL);

private:
  Converter *c;

  template <class TOutPixel> 
    void TemplatedWriteImage(const char *file, double xRoundFactor, int pos,const char *history=NULL);

  template <class TOutPixel> 
    void TemplatedWriteMultiComponentImage(const char *file, double xRoundFactor, int pstart,const char *history=NULL);

  void MakeSPMOriginFix(itk::ImageBase<VDim> *input, itk::ImageBase<VDim> *output);

};

#endif

