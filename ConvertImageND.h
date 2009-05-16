#ifndef __ConvertImageND_h_
#define __ConvertImageND_h_

#include "itkOrientedRASImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkInterpolateImageFunction.h"

#include <iostream>
#include <exception>
#include <cctype>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <cstdarg>

using namespace std;

template <class TPixel, unsigned int VDim> class ConvertAdapter;

class ConvertException : public std::exception
{
public:
  ConvertException(const char *fmt, ...)
    {
    char buffer[1024];
    va_list parg;
    va_start(parg, fmt);
    vsprintf(buffer, fmt, parg);
    va_end(parg);
    message=buffer;
    }

  virtual ~ConvertException() throw() {}

  virtual const char *what() const throw()
    { return message.c_str(); }

private:
  string message;
};

template<class TPixel, unsigned int VDim>
class ImageConverter
{
public:

  // Image typedef
  typedef itk::OrientedRASImage<TPixel, VDim> ImageType;
  typedef itk::Image<TPixel, VDim> UnorientedImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::RegionType RegionType;
  typedef vnl_vector_fixed<double, VDim> RealVector;
  typedef vnl_vector_fixed<int, VDim> IntegerVector;

  // Complex stuff
  typedef std::complex<TPixel> ComplexPixel;
  typedef itk::OrientedRASImage<ComplexPixel, VDim> ComplexImageType;

  // Iterators
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> ConstIterator;

  // Interpolator
  typedef itk::InterpolateImageFunction<ImageType, double> Interpolator;
  
  ImageConverter();
  int ProcessCommandLine(int argc, char *argv[]);

  friend class ConvertAdapter<TPixel, VDim>;

  // Copy image on stack
  void CopyImage();

  // Get bounding box of an image
  void GetBoundingBox(ImageType *image, RealVector &bb0, RealVector &bb1);

  // Print a matrix in a nice way
  void PrintMatrix(
    std::ostream &sout, vnl_matrix<double> mat, 
    const char *fmt = "%10.4f ", const char *prefix = "");

  // Set label set for split/merge operations
  typedef std::vector<double> LabelSet;
  LabelSet &GetSplitLabelSet() 
    { return m_SplitLabelSet; }

  Interpolator *GetInterpolator() const
    { return m_Interpolator; }

  void SetInterpolator(Interpolator *interp)
    { m_Interpolator = interp; }

private:

  // Internal functions
  void PrintCommandListing();
  int ProcessCommand(int argc, char *argv[]);

  // Read vectors, etc from command line
  SizeType ReadSizeVector(const char *vec);
  IndexType ReadIndexVector(const char *vec);
  RealVector ReadRealVector(const char *vec, bool is_point);
  RealVector ReadRealSize(const char *vec)
    {
    RealVector x = ReadRealVector(vec, false);
    for(size_t d = 0; d < VDim; d++)
      x[d] = fabs(x[d]);
    return x;
    }
  TPixel ReadIntensityValue(const char *vec);

  // Templated write function
  template<class TOutPixel> void TemplatedWriteImage(const char *file, double xRoundFactor);

  // Map of variable names
  typedef map<string, ImagePointer> ImageVariableMap;
  ImageVariableMap m_ImageVars;

  // Image interpolator for all interpolation commands
  itk::SmartPointer<Interpolator> m_Interpolator;

  // Implementation of the 'foreach' loop
  size_t ForEachLoop(int argc, char *argv[]);

  // Loop status flag
  bool m_InLoop;

  // Label set for split/merge
  LabelSet m_SplitLabelSet;

public:

  // Stack of images from the command line
  vector<ImagePointer> m_ImageStack;

  // Typeid of the image to be saved
  string m_TypeId;

  // Interpolation type
  string m_Interpolation;

  // Background intensity
  double m_Background;

  // Rounding factor (not used for float and double) is 0.5 unless -noround was specified
  double m_RoundFactor;

  // Whether SPM extensions are used
  bool m_FlagSPM;

  // Number of iterations for various algorithms
  size_t m_Iterations;

  // Root mean square error for anti-aliasing algorithm
  double m_AntiAliasRMS;

  // Level set algorithm parameters
  double m_LevSetCurvature, m_LevSetAdvection;

  // Verbose output stream
  std::ostringstream devnull;
  std::ostream *verbose;
};

#endif

