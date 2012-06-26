#include "ConvertImageND.h"
#include "itkVoxBoCUBImageIOFactory.h"
#include "itkPovRayDF3ImageIOFactory.h"

#ifdef HAVE_MINC4ITK
#include "itkMincImageIOFactory.h"
#include <time_stamp.h>    // for creating minc style history entry
#endif //HAVE_MINC4ITK


int main(int argc, char *argv[])
{
  // Load the ITK factories
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());
  itk::ObjectFactoryBase::RegisterFactory(itk::PovRayDF3ImageIOFactory::New());
  const char *history=NULL;
#ifdef HAVE_MINC4ITK 
  itk::RegisterMincIO();
  history = time_stamp(argc, argv); 
#endif //HAVE_MINC4ITK

  ImageConverter<double, 3> convert;
  int r=convert.ProcessCommandLine(argc, argv,history);
  
  if(history) free((void*)history);
    
  return r;
}
