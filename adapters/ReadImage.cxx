#include "ReadImage.h"
#include "itkIOCommon.h"
#include "itkImageFileReader.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"

template <class TPixel, unsigned int VDim>
void
ReadImage<TPixel, VDim>
::operator() (const char *file)
{
  // Set up the reader
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(file);

  // Report
  *c->verbose << "Reading #" << (1 + c->m_ImageStack.size()) << " from " << file << endl;
  
  // Try reading this file
  try { reader->Update(); }
  catch(itk::ExceptionObject &exc)
    {
    cerr << "Error reading image: " << file << endl;
    cerr << "ITK Exception: " << exc << endl;
    throw -1;
    }
  ImagePointer image = reader->GetOutput();

  // TODO: get rid of this code, nifti is more reliable
  string ext = itksys::SystemTools::GetFilenameExtension(file);
  if((ext == ".hdr" || ext == ".img.gz" || ext == ".img") && c->m_FlagSPM)
    {
    string temp;
    if(itk::ExposeMetaData<std::string>(
        image->GetMetaDataDictionary(), itk::ITK_FileOriginator, temp))
      {
      typename ImageType::PointType oitk = image->GetOrigin();
      typename ImageType::SpacingType sitk = image->GetSpacing();
      
      // Read the SPM-style origin
      *c->verbose << "  Applying SPM origin :";
      for(size_t i=0; i < VDim; i++)
        {
        short xspm = (temp[2*i+1] << 8) + temp[2*i];
        oitk[i] = -sitk[i] * xspm;
        *c->verbose << xspm << " ";
        }
      *c->verbose << endl;

      image->SetOrigin(oitk);
      }
    }




  // Push the image into the stack
  c->m_ImageStack.push_back(image);
}

// Invocations
template class ReadImage<double, 2>;
template class ReadImage<double, 3>;
