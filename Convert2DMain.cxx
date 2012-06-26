#include "ConvertImageND.h"
#include "itkVoxBoCUBImageIOFactory.h"
#include "itkPovRayDF3ImageIOFactory.h"
#ifdef HAVE_MINC4ITK
#include "itkMincImageIOFactory.h"
#endif //HAVE_MINC4ITK

int main(int argc, char *argv[])
{
  // Load the ITK factories
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());
  itk::ObjectFactoryBase::RegisterFactory(itk::PovRayDF3ImageIOFactory::New());
#ifdef HAVE_MINC4ITK 
  itk::RegisterMincIO();
#endif //HAVE_MINC4ITK

  ImageConverter<double, 2> convert;
  return convert.ProcessCommandLine(argc, argv);
}

