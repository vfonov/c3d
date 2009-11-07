#include "LabelOverlapMeasures.h"
#include "itkLabelOverlapMeasuresImageFilter.h"

#include <iomanip>

template <class TPixel, unsigned int VDim>
void
LabelOverlapMeasures<TPixel, VDim>
::operator() ()
{
  // Get images from stack
  size_t n = c->m_ImageStack.size();
  if( n < 2 )
    throw ConvertException( "Label overlap measures require two image inputs" );
  ImagePointer target = c->m_ImageStack[n-1];
  ImagePointer source = c->m_ImageStack[n-2];

  // Create a short image for the labels
  typedef itk::Image<short, VDim> LabelImageType;
  typename LabelImageType::Pointer slabSource = LabelImageType::New();
  typename LabelImageType::Pointer slabTarget = LabelImageType::New();

  // Allocate the images
  slabSource->SetRegions( source->GetBufferedRegion() );
  slabSource->Allocate();

  slabTarget->SetRegions( target->GetBufferedRegion() );
  slabTarget->Allocate();

  // Round off doubles to create labels
  size_t nvS = slabSource->GetBufferedRegion().GetNumberOfPixels();
  for( size_t i = 0; i < nvS; i++ )
    {
    slabSource->GetBufferPointer()[i]
      = (short) ( source->GetBufferPointer()[i] + 0.5 );
    }

  size_t nvT = slabTarget->GetBufferedRegion().GetNumberOfPixels();
  for( size_t i = 0; i < nvT; i++ )
    {
    slabTarget->GetBufferPointer()[i]
      = (short) ( target->GetBufferPointer()[i] + 0.5 );
    }

  // Create the label statistics fltOverlap
  typedef itk::LabelOverlapMeasuresImageFilter<LabelImageType> OverlapFilter;
  typename OverlapFilter::Pointer fltOverlap = OverlapFilter::New();

  // Set the inputs
  fltOverlap->SetSourceImage( slabSource );
  fltOverlap->SetTargetImage( slabTarget );

  // Update the fltOverlap
  fltOverlap->Update();

  // Print the result

  std::cout << "                                  "
            << "************  All Labels *************" << std::endl;
  std::cout << std::setw( 17 ) << "   " << std::setw( 17 ) << "Total"
    << std::setw( 17 ) << "Union (jaccard)"
    << std::setw( 17 ) << "Mean (dice)" << std::setw( 17 ) << "False negative"
    << std::setw( 17 ) << "False positive" << std::endl;
  std::cout << std::setw( 17 ) << "Volumetric";
  std::cout << std::setw( 17 ) << fltOverlap->GetVolumeTotalOverlap();
  std::cout << std::setw( 17 ) << fltOverlap->GetVolumeUnionOverlap();
  std::cout << std::setw( 17 ) << fltOverlap->GetVolumeMeanOverlap();
  std::cout << std::setw( 17 ) << fltOverlap->GetVolumeFalseNegativeError();
  std::cout << std::setw( 17 ) << fltOverlap->GetVolumeFalsePositiveError();
  std::cout << std::endl;
  std::cout << std::setw( 17 ) << "Surface";
  std::cout << std::setw( 17 ) << fltOverlap->GetSurfaceTotalOverlap();
  std::cout << std::setw( 17 ) << fltOverlap->GetSurfaceUnionOverlap();
  std::cout << std::setw( 17 ) << fltOverlap->GetSurfaceMeanOverlap();
  std::cout << std::setw( 17 ) << fltOverlap->GetSurfaceFalseNegativeError();
  std::cout << std::setw( 17 ) << fltOverlap->GetSurfaceFalsePositiveError();
  std::cout << std::endl;
  std::cout << "Other related measures." << std::endl;
  std::cout << "  Volume similarity: "
    << fltOverlap->GetVolumeSimilarity() << std::endl;
  std::cout << std::endl;

  std::cout << "                               "
            << "************ Individual Labels *************" << std::endl;
  std::cout << "Volumetric overlap measures." << std::endl;
  std::cout << std::setw( 17 ) << "Label"
            << std::setw( 17 ) << "Total"
            << std::setw( 17 ) << "Union (jaccard)"
            << std::setw( 17 ) << "Mean (dice)"
            << std::setw( 17 ) << "False negative"
            << std::setw( 17 ) << "False positive" << std::endl;

  typename OverlapFilter::MapType labelMap = fltOverlap->GetLabelSetMeasures();
  typename OverlapFilter::MapType::const_iterator it;
  for( it = labelMap.begin(); it != labelMap.end(); ++it )
    {
    if( (*it).first == 0 )
      {
      continue;
      }

    short label = (*it).first;

    std::cout << std::setw( 17 ) << label;
    std::cout << std::setw( 17 ) << fltOverlap->GetVolumeTotalOverlap( label );
    std::cout << std::setw( 17 ) << fltOverlap->GetVolumeUnionOverlap( label );
    std::cout << std::setw( 17 ) << fltOverlap->GetVolumeMeanOverlap( label );
    std::cout << std::setw( 17 ) << fltOverlap->GetVolumeFalseNegativeError( label );
    std::cout << std::setw( 17 ) << fltOverlap->GetVolumeFalsePositiveError( label );
    std::cout << std::endl;
    }

  std::cout << "Surface overlap measures." << std::endl;
  std::cout << std::setw( 17 ) << "Label"
            << std::setw( 17 ) << "Total"
            << std::setw( 17 ) << "Union (jaccard)"
            << std::setw( 17 ) << "Mean (dice)"
            << std::setw( 17 ) << "False negative"
            << std::setw( 17 ) << "False positive" << std::endl;
  for( it = labelMap.begin(); it != labelMap.end(); ++it )
    {
    if( (*it).first == 0 )
      {
      continue;
      }
    short label = (*it).first;

    std::cout << std::setw( 17 ) << label;
    std::cout << std::setw( 17 ) << fltOverlap->GetSurfaceTotalOverlap( label );
    std::cout << std::setw( 17 ) << fltOverlap->GetSurfaceUnionOverlap( label );
    std::cout << std::setw( 17 ) << fltOverlap->GetSurfaceMeanOverlap( label );
    std::cout << std::setw( 17 ) << fltOverlap->GetSurfaceFalseNegativeError( label );
    std::cout << std::setw( 17 ) << fltOverlap->GetSurfaceFalsePositiveError( label );
    std::cout << std::endl;
    }

  std::cout << "Other related measures." << std::endl;
  std::cout << std::setw( 17 ) << "Label"
            << std::setw( 17 ) << "Volume sim."
            << std::setw( 17 ) << "Hausdorff"
            << std::setw( 17 ) << "Dir. Hausdorff"
            << std::setw( 17 ) << "Contour"
            << std::setw( 17 ) << "Dir. Contour" << std::endl;
  for( it = labelMap.begin(); it != labelMap.end(); ++it )
    {
    if( (*it).first == 0 )
      {
      continue;
      }
    short label = (*it).first;

    std::cout << std::setw( 17 ) << label;
    std::cout << std::setw( 17 ) << fltOverlap->GetVolumeSimilarity( label );
    std::cout << std::setw( 17 ) << fltOverlap->GetHausdorffDistance( label );
    std::cout << std::setw( 17 ) << fltOverlap->GetDirectedHausdorffDistance( label );
    std::cout << std::setw( 17 ) << fltOverlap->GetContourMeanDistance( label );
    std::cout << std::setw( 17 ) << fltOverlap->GetDirectedContourMeanDistance( label );
    std::cout << std::endl;
    }
}

// Invocations
template class LabelOverlapMeasures<double, 2>;
template class LabelOverlapMeasures<double, 3>;
