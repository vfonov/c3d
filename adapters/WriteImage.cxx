#include "WriteImage.h"
#include "itksys/SystemTools.hxx"
#include "itkImageFileWriter.h"
#include "itkIOCommon.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"

template<class TPixel, unsigned int VDim>
template<class TOutPixel>
void 
WriteImage<TPixel, VDim>
::TemplatedWriteImage(const char *file, double xRoundFactor)
{
  // Get the input image
  if(c->m_ImageStack.size() == 0)
    { cerr << "No data has been generated! Can't write to " << file << endl; throw -1; }
  ImagePointer input = c->m_ImageStack.back();
  
  // Create the output image 
  typedef itk::OrientedRASImage<TOutPixel, VDim> OutputImageType;
  typename OutputImageType::Pointer output = OutputImageType::New();
  output->SetRegions(input->GetBufferedRegion());
  output->SetSpacing(input->GetSpacing());
  output->SetOrigin(input->GetOrigin());
  output->SetDirection(input->GetDirection());
  output->SetMetaDataDictionary(input->GetMetaDataDictionary());
  output->Allocate();

  // Describe what we are doing
  *c->verbose << "Writing #" << c->m_ImageStack.size() << " to file " << file << endl;
  *c->verbose << "  Output voxel type: " << c->m_TypeId << "[" << typeid(TOutPixel).name() << "]" << endl;
  *c->verbose << "  Rounding off: " << (xRoundFactor == 0.0 ? "Disabled" : "Enabled") << endl;
  
  // Set the SPM originator header
  if(c->m_FlagSPM)
    {
    size_t i;
    string originator;
    originator.resize(VDim * 2);
    
    // Compute the SPM-style origin of the image
    *c->verbose << "  Setting SPM origin field to:";
    for(i = 0; i < VDim; i++)
      {
      short ospm = (short)(0.5 - input->GetOrigin()[i] / input->GetSpacing()[i]);
      originator[2*i] = (char)(ospm & 0x00ff);
      originator[2*i+1] = (char)(ospm >> 8);
      *c->verbose << ospm << " ";
      }
    originator[2*i] = 0;
    *c->verbose << endl;

    itk::EncapsulateMetaData<string>(
      output->GetMetaDataDictionary(),itk::ITK_FileOriginator,originator);
    }

  // Copy everything, rounding if the pixel type is integer
  size_t n = input->GetBufferedRegion().GetNumberOfPixels();
  for(size_t i = 0; i < n; i++)
    output->GetBufferPointer()[i] = (TOutPixel) (input->GetBufferPointer()[i] + xRoundFactor);

  // Write the image out
  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(output);
  writer->SetFileName(file);
  try { writer->Update(); }
  catch (itk::ExceptionObject &exc) {
    cerr << "Error writing image to " << file << endl;
    cerr << "ITK Exception: " << exc << endl;
    throw -1;
  }
}

template <class TPixel, unsigned int VDim>
void
WriteImage<TPixel, VDim>
::operator() (const char *file, bool force)
{
  // Unless in 'force' mode, check if the image already exists
  if(!force && itksys::SystemTools::FileExists(file))
    {
    cerr << "File " << file << " already exists. Use -o option to override!" << endl;
    throw -1;
    }

  if(c->m_TypeId == "char" || c->m_TypeId == "byte")
    TemplatedWriteImage<char>(file, c->m_RoundFactor);
  if(c->m_TypeId == "uchar" || c->m_TypeId == "ubyte")
    TemplatedWriteImage<unsigned char>(file, c->m_RoundFactor);
  
  if(c->m_TypeId == "short") 
    TemplatedWriteImage<short>(file, c->m_RoundFactor);
  if(c->m_TypeId == "ushort")
    TemplatedWriteImage<unsigned short>(file, c->m_RoundFactor);

  if(c->m_TypeId == "int") 
    TemplatedWriteImage<int>(file, c->m_RoundFactor);
  if(c->m_TypeId == "uint")
    TemplatedWriteImage<unsigned int>(file, c->m_RoundFactor);

  if(c->m_TypeId == "float") 
    TemplatedWriteImage<float>(file, 0.0);
  if(c->m_TypeId == "double")
    TemplatedWriteImage<double>(file, 0.0);

}

// Invocations
template class WriteImage<double, 2>;
template class WriteImage<double, 3>;
