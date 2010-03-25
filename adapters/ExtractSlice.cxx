#include "ExtractSlice.h"
#include <string>
#include <iostream>
#include "itkExtractImageFilter.h"
#include "itkImageFileWriter.h"


template <class TPixel, unsigned int VDim>
void
ExtractSlice<TPixel, VDim>
::operator() (string axis, char* pos, char* filename)
{
  // Check input availability
  if(c->m_ImageStack.size() < 1)
    {
    cerr << "No images on the stack" << endl;
    throw 0;
    }

  // Get the image
  ImagePointer image = c->m_ImageStack.back();
  typename ImageType::SpacingType spacing = image->GetSpacing();
  typename ImageType::RegionType inputRegion = image->GetLargestPossibleRegion();
  typename ImageType::SizeType size = inputRegion.GetSize();
  typedef itk::Image<TPixel, 2> SliceType;

  typename SliceType::SpacingType slicespacing;
  for (int ii=0; ii<2; ii++) slicespacing[ii]=spacing[ii];
  unsigned int slicedir;
  if (!axis.compare("x"))
    {
    slicedir = 0;
    slicespacing[0] = spacing[1];
    slicespacing[1] = spacing[2];
    }
  else if (!axis.compare("y"))
    {
    slicedir = 1;
    slicespacing[0] = spacing[0];
    slicespacing[1] = spacing[2];
    }
  else if (!axis.compare("z"))
    {
    slicedir = 2;
    slicespacing[0] = spacing[0];
    slicespacing[1] = spacing[1];
    }
  else
    {
    cerr << "axis must be x,y or z" << endl;
    throw 0;
    }

  double percent_pos;
  int slicepos;
  std::string s( pos );
  size_t ipos = s.rfind("%");
  if(ipos == s.size() - 1)
    {
    const char *tok = strtok(pos, "%");
    percent_pos = atof( tok );
    slicepos = (int)round( ( percent_pos / 100.0 ) * (size[slicedir] -1)); 
    }
  else
    slicepos = atoi( pos );

  // Say what we are doing
  *c->verbose << "Extracting slice " << pos << " along " << axis << " axis " << endl;


  typedef itk::ExtractImageFilter< ImageType, SliceType > ExtractFilterType;
  typename ExtractFilterType::Pointer extractfilter = ExtractFilterType::New();
  typename ImageType::PointType origin = image->GetOrigin();
  size[slicedir] = 0;
  typename ImageType::IndexType start = inputRegion.GetIndex();
  start[slicedir] = slicepos;
  typename ImageType::RegionType desiredRegion;
  desiredRegion.SetSize(  size  );
  desiredRegion.SetIndex( start );
  extractfilter->SetExtractionRegion( desiredRegion );
  extractfilter->SetInput( image );
  extractfilter->Update();
  typename SliceType::Pointer slice = extractfilter->GetOutput();
  slice->SetSpacing(slicespacing);

  typedef itk::ImageFileWriter<SliceType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( slice );
  writer->SetFileName( filename );
  writer->Update();

}

// Invocations
template class ExtractSlice<double, 3>;
template class ExtractSlice<double, 2>;

