#include "BiasFieldCorrection.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkN3MRIBiasFieldCorrectionImageFilter.h"

template <class TPixel, unsigned int VDim>
void
BiasFieldCorrection<TPixel, VDim>
::operator() ()
{
  // Check input availability
  if(c->m_ImageStack.size() < 1)
    throw ConvertException("No images on stack");

  int n3_shrink_factor=4;
  float  n3_spline_distance = 100;
  
  
  // Get image from stack
  ImagePointer mri = c->m_ImageStack.back();
  c->m_ImageStack.pop_back();

  // Report what the filter is doing
  *c->verbose << "N3 BiasFieldCorrection #" << c->m_ImageStack.size() << endl;
  
  
  // Set up a filter to shrink image by a factor
  typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkerType;
  typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput(mri);
  shrinker->SetShrinkFactors(n3_shrink_factor);
  shrinker->Update();

  // Compute mask using Otsu threshold
  typedef itk::OtsuThresholdImageFilter<ImageType, ImageType> ThresholderType;
  typename ThresholderType::Pointer otsu = ThresholderType::New();
  otsu->SetInput(mri);
  otsu->SetNumberOfHistogramBins( 200 );
  otsu->SetInsideValue( 0 );
  otsu->SetOutsideValue( 1 );
  otsu->Update();
  ImagePointer mask = otsu->GetOutput();
  *c->verbose << "   Otsu threshold: "<<otsu->GetThreshold()<<endl;
  

  // Shrink the mask
  typename ShrinkerType::Pointer maskshrinker = ShrinkerType::New();
  maskshrinker->SetInput( mask);
  maskshrinker->SetShrinkFactors(n3_shrink_factor);
  maskshrinker->Update();

  // Bias filter
  typedef itk::N3MRIBiasFieldCorrectionImageFilter<ImageType, ImageType, ImageType> CorrecterType;
  typename CorrecterType::Pointer correcter = CorrecterType::New();
  correcter->SetInput( shrinker->GetOutput() );
  correcter->SetMaskImage( maskshrinker->GetOutput() );

  typename CorrecterType::ArrayType numberOfControlPoints;

  typename ImageType::IndexType inputImageIndex =
    mri->GetLargestPossibleRegion().GetIndex();
  typename ImageType::SizeType inputImageSize =
    mri->GetLargestPossibleRegion().GetSize();

  typename ImageType::PointType newOrigin = mri->GetOrigin();

  unsigned long lowerBound[VDim];
  unsigned long upperBound[VDim];

  for( unsigned int d = 0; d < VDim; d++ )
    {
    float domain = static_cast<float>( mri->
      GetLargestPossibleRegion().GetSize()[d] - 1 ) * mri->GetSpacing()[d];
    unsigned int numberOfSpans = static_cast<unsigned int>(
      vcl_ceil( domain / n3_spline_distance ) );
    unsigned long extraPadding = static_cast<unsigned long>( ( numberOfSpans *
      n3_spline_distance - domain ) / mri->GetSpacing()[d] + 0.5 );
    lowerBound[d] = static_cast<unsigned long>( 0.5 * extraPadding );
    upperBound[d] = extraPadding - lowerBound[d];
    newOrigin[d] -= ( static_cast<float>( lowerBound[d] ) *
      mri->GetSpacing()[d] );

    numberOfControlPoints[d] = numberOfSpans + correcter->GetSplineOrder();
    }
  correcter->SetNumberOfControlPoints( numberOfControlPoints );
  
  
  *c->verbose << "  Shrink factor: " << n3_shrink_factor << endl;
  *c->verbose << "  Number of histogram bins: "<< correcter->GetNumberOfHistogramBins()<<endl;
  *c->verbose << "  Weiner Filter noise: "<< correcter->GetWeinerFilterNoise()<<endl;
  *c->verbose << "  Bias FWHM: "<< correcter->GetBiasFieldFullWidthAtHalfMaximum()<<endl;
  *c->verbose << "  Max Number of Iterations: "<< correcter->GetMaximumNumberOfIterations()<<endl;
  *c->verbose << "  Convergence Threshold: "<< correcter->GetConvergenceThreshold()<<endl;
  *c->verbose << "  Spline Order: "<< correcter->GetSplineOrder()<<endl;
  *c->verbose << "  Use Optimal Bias Field scaling: "<< correcter->GetUseOptimalBiasFieldScaling()<<endl;

  // Progress meter
  // typedef CommandIterationUpdate<CorrecterType> CommandType;
  // typename CommandType::Pointer observer = CommandType::New();
  // correcter->AddObserver( itk::IterationEvent(), observer );
  correcter->Update();

  /**
   * Reconstruct the bias field at full image resolution.  Divide
   * the original input image by the bias field to get the final
   * corrected image.
   */
  typedef itk::BSplineControlPointImageFilter<typename
    CorrecterType::BiasFieldControlPointLatticeType, typename
    CorrecterType::ScalarImageType> BSplinerType;

  typename BSplinerType::Pointer bspliner = BSplinerType::New();
  bspliner->SetInput( correcter->GetLogBiasFieldControlPointLattice() );
  bspliner->SetSplineOrder( correcter->GetSplineOrder() );
  bspliner->SetSize( mri->GetLargestPossibleRegion().GetSize() );
  bspliner->SetOrigin( mri->GetOrigin() );
  bspliner->SetDirection( mri->GetDirection() );
  bspliner->SetSpacing( mri->GetSpacing() );
  bspliner->Update();

  typename ImageType::Pointer logField = ImageType::New();
  logField->SetOrigin( bspliner->GetOutput()->GetOrigin() );
  logField->SetSpacing( bspliner->GetOutput()->GetSpacing() );
  logField->SetRegions( bspliner->GetOutput()->GetLargestPossibleRegion().GetSize() );
  logField->SetDirection( bspliner->GetOutput()->GetDirection() );
  logField->Allocate();

  itk::ImageRegionIterator<typename CorrecterType::ScalarImageType> ItB(
    bspliner->GetOutput(),
    bspliner->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItF( logField,
    logField->GetLargestPossibleRegion() );
  for( ItB.GoToBegin(), ItF.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItF )
    {
    ItF.Set( ItB.Get()[0] );
    }

  typedef itk::ExpImageFilter<ImageType, ImageType> ExpFilterType;
  typename ExpFilterType::Pointer expFilter = ExpFilterType::New();
  expFilter->SetInput( logField );
  expFilter->Update();

  typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> DividerType;
  typename DividerType::Pointer divider = DividerType::New();
  divider->SetInput1( mri );
  divider->SetInput2( expFilter->GetOutput() );
  divider->Update();

  // Update
  c->m_ImageStack.push_back(divider->GetOutput());
}

// Invocations
template class BiasFieldCorrection<double, 2>;
template class BiasFieldCorrection<double, 3>;
