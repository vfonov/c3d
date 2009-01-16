#include "PadImage.h"
#include <itkConstantPadImageFilter.h>

template <class TPixel, unsigned int VDim>
void
PadImage<TPixel, VDim>
::operator() (IndexType padExtentLower, IndexType padExtentUpper, float padValue)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  typedef itk::ConstantPadImageFilter< ImageType, ImageType > ConstantPad;
  typename ConstantPad::Pointer padFilter = ConstantPad::New();

  // Pad first three dimensions only
  unsigned long lowerBound[VDim];
  unsigned long upperBound[VDim];

  for (int i = 0; i < 3; i++) {
    lowerBound[i] = padExtentLower[i];
    upperBound[i] = padExtentUpper[i];
  }

  for (unsigned int i = 3; i < VDim; i++) {
    lowerBound[i] = 0;
    upperBound[i] = 0;
  }

  padFilter->SetPadLowerBound(lowerBound);
  padFilter->SetPadUpperBound(upperBound);

  padFilter->SetConstant(static_cast<typename ImageType::PixelType>(padValue));

  padFilter->SetInput(img);

  padFilter->Update();

  ImagePointer output = padFilter->GetOutput();

  typename ImageType::IndexType originVox, paddedOriginVox;

  originVox.Fill(0);

  paddedOriginVox[0] = lowerBound[0];
  paddedOriginVox[1] = lowerBound[1];
  paddedOriginVox[2] = lowerBound[2];

  itk::Point<double, VDim> originMM, originMM_Padded; 

  
  
  img->TransformIndexToPhysicalPoint(originVox, originMM);
  output->TransformIndexToPhysicalPoint(paddedOriginVox, originMM_Padded);

  

  RealVector origin;
  for(size_t i = 0; i < 3; i++)
    {
      origin[i] = img->GetOrigin()[i] + (originMM[i] - originMM_Padded[i]);
    }

  padFilter->GetOutput()->SetOrigin(origin.data_block());

  cout << "INP_ORG: " << img->GetOrigin() << endl;
  cout << "OUT_ORG: " << output->GetOrigin() << endl;

  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(padFilter->GetOutput());

}

// Invocations
template class PadImage<double, 2>;
template class PadImage<double, 3>;
