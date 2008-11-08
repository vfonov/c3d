#include "ConvertImageND.h"

#include "AddImages.h"
#include "AntiAliasImage.h"
#include "ApplyMetric.h"
#include "ClipImageIntensity.h"
#include "ComputeFFT.h"
#include "ComputeOverlaps.h"
#include "ConnectedComponents.h"
#include "ConvertAdapter.h"
#include "CopyTransform.h"
#include "CreateImage.h"
#include "ExtractRegion.h"
#include "GeneralLinearModel.h"
#include "ImageERF.h"
#include "ImageLaplacian.h"
#include "LabelStatistics.h"
#include "LevelSetSegmentation.h"
#include "MathematicalMorphology.h"
#include "MultiplyImages.h"
#include "PeronaMalik.h"
#include "PrintImageInfo.h"
#include "ReadImage.h"
#include "ReciprocalImage.h"
#include "ReplaceIntensities.h"
#include "ResampleImage.h"
#include "ResliceImage.h"
#include "SampleImage.h"
#include "ScaleShiftImage.h"
#include "SignedDistanceTransform.h"
#include "SmoothImage.h"
#include "SplitMultilabelImage.h"
#include "StapleAlgorithm.h"
#include "ThresholdImage.h"
#include "TrimImage.h"
#include "UpdateMetadataKey.h"
#include "Vote.h"
#include "WarpImage.h"
#include "WarpLabelImage.h"
#include "WriteImage.h"

#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"

#include <cstring>

// Helper function: read a double, throw exception if unreadable
double myatof(char *str)
{
  char *end = 0;
  double d = strtod(str, &end);
  if(*end != 0)
    throw "strtod conversion failed";
  return d;
};


template<class TPixel, unsigned int VDim>
ImageConverter<TPixel,VDim>
::ImageConverter()
  : verbose(&devnull)
{
  // Initialize to defaults
  m_TypeId = "float";
  m_Interpolation = "Linear";
  m_Background = 0.0;
  m_RoundFactor = 0.5;
  m_FlagSPM = false;
  m_AntiAliasRMS = 0.07;
  m_Iterations = 0;
  m_InLoop = false;

  m_LevSetCurvature = 0.2;
  m_LevSetAdvection = 0.0;

  SetInterpolator();
}

template<class TPixel, unsigned int VDim>
int
ImageConverter<TPixel, VDim>
::ProcessCommand(int argc, char *argv[])
{
  // Get the first command
  string cmd = argv[0];

  // cout << "COMMAND: " << cmd << endl;

  // Commands in alphabetical order
  if(cmd == "-add")
    {
    AddImages<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if(cmd == "-anisotropic-diffusion" || cmd == "-ad")
    {
    double cond = atof(argv[1]);
    int niter = atoi(argv[2]);
    PeronaMalik<TPixel, VDim> adapter(this);
    adapter(cond, (size_t) niter);
    return 2;
    }

  // Anti-alias a binary image, turning it into a smoother floating point image;
  // the argument is the iso-surface value
  // This command is affected by -iterations and -rms flags
  else if(cmd == "-antialias" || cmd == "-alias")
    {
    AntiAliasImage<TPixel, VDim> adapter(this);
    adapter(atof(argv[1]));
    return 1;
    }
   
  // Associate variable name with image currently at the top of
  // the stack
  else if(cmd == "-as" || cmd == "-set")
    {
    string var(argv[1]);
    if(m_ImageStack.size() == 0)
      throw ConvertException("No image to assign to variable %s", var.c_str());
    m_ImageVars[var] = m_ImageStack.back();
    return 1;
    }

  else if(cmd == "-background")
    { m_Background = atof(argv[1]); return 1; }

  // f(x) = (x == xBackground) ? 0 : 1
  else if(cmd == "-binarize")
    {
    ThresholdImage<TPixel, VDim> adapter(this);
    adapter(m_Background, m_Background, 0.0, 1.0);
    return 0;
    }

  else if(cmd == "-connected-components" || cmd == "-connected" || cmd == "-comp")
    {
    ConnectedComponents<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if(cmd == "-clear")
    {
    *verbose << "Clearing the stack" << endl;
    m_ImageStack.clear();
    return 0;
    }

  else if(cmd == "-clip")
    {
    double iMin = ReadIntensityValue(argv[1]);
    double iMax = ReadIntensityValue(argv[2]);
    ClipImageIntensity<TPixel, VDim>(this)(iMin, iMax);
    return 2;
    }

  else if(cmd == "-copy-transform" || cmd == "-ct")
    {
    CopyTransform<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  // Create a new image with given size and voxel size
  else if(cmd == "-create")
    {
    SizeType dims = ReadSizeVector(argv[1]);
    RealVector voxel = ReadRealVector(argv[2]);
    CreateImage<TPixel, VDim> adapter(this);
    adapter(dims, voxel);
    return 2;
    }

  else if(cmd == "-dilate")
    {
    MathematicalMorphology<TPixel,VDim> adapter(this);
    adapter(false, atof(argv[1]), ReadSizeVector(argv[2]));
    return 2;
    }

  else if(cmd == "-divide")
    {
    // Multiply by the reciprocal (for time being at least)
    ReciprocalImage<TPixel, VDim> reciprocal(this);
    reciprocal();

    MultiplyImages<TPixel, VDim> multiply(this);
    multiply();

    return 0;
    }

  else if(cmd == "-endfor")
    {
    // This command ends the for loop
    if(!m_InLoop)
      throw ConvertException("Out of place -endfor command");
    this->m_InLoop = false;
    return 0;
    }

  else if(cmd == "-erode")
    {
    MathematicalMorphology<TPixel,VDim> adapter(this);
    adapter(false, atof(argv[1]), ReadSizeVector(argv[2]));
    return 2;
    }

  else if(cmd == "-erf")
    {
    double thresh = atof(argv[1]);
    double scale = atof(argv[2]);
    ImageERF<TPixel, VDim> adapter(this);
    adapter(thresh, scale);
    return 2;
    }

  else if(cmd == "-fft")
    {
    ComputeFFT<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if(cmd == "-foreach")    
    {
    if(m_InLoop)
      throw ConvertException("Nested loops are not allowed");
    this->m_InLoop = true;
    return this->ForEachLoop(argc, argv);
    }

  else if (cmd == "-glm")
    {
    string mat(argv[1]);
    string con(argv[2]);
    GeneralLinearModel<TPixel, VDim> adapter(this);
    adapter(mat, con);
    return 2;
    }

  else if(cmd == "-info")
    { 
    PrintImageInfo<TPixel, VDim> adapter(this);
    adapter(false); 
    return 0; 
    }

  else if(cmd == "-info-full")
    { 
    PrintImageInfo<TPixel, VDim> adapter(this);
    adapter(true); 
    return 0; 
    }

  else if (cmd == "-insert" || cmd == "-ins")
    {
    string var(argv[1]);
    size_t pos = (size_t) atoi(argv[2]);
    typename ImageVariableMap::iterator img = m_ImageVars.find(var);

    // Check if the variable exists
    if(img == m_ImageVars.end())
      throw ConvertException("No image assigned to variable %s", var.c_str());

    // Check if the position is ok
    if(m_ImageStack.size() > pos)
      throw ConvertException("Can not insert at position %i in stack of size %i", pos, m_ImageStack.size());

    // Insert at the appropriate place
    typename vector<ImagePointer>::iterator it = m_ImageStack.end();
    for(size_t i = 0; i < pos; i++) --it;
    m_ImageStack.insert(it, img->second);
    
    return 2;
    }

  else if(cmd == "-interpolation" || cmd == "-interp" || cmd == "-int")
    { 
    m_Interpolation = argv[1];
    SetInterpolator();
    return 1; 
    }

  else if(cmd == "-iterations")
    { 
    m_Iterations = static_cast<size_t>(atoi(argv[1])); 
    return 1; 
    }

  else if(cmd == "-label-statistics" || cmd == "-lstat")
    {
    LabelStatistics<TPixel, VDim>(this)();
    return 0;
    }

  else if(cmd == "-laplacian" || cmd == "-laplace")
    {
    ImageLaplacian<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if(cmd == "-levelset")
    {
    int nIter = atoi(argv[1]);
    LevelSetSegmentation<TPixel, VDim> adapter(this);
    adapter(nIter);
    return 1;
    }

  else if(cmd == "-levelset-curvature")
    {
    m_LevSetCurvature = atof(argv[1]);
    return 1;
    }
  
  else if(cmd == "-levelset-advection")
    {
    m_LevSetAdvection = atof(argv[1]);
    return 1;
    }

  else if(cmd == "-merge")
    {
    Vote<TPixel, VDim> adapter(this);
    adapter(true);
    return 0;
    }

  else if(cmd == "-mi" || cmd == "-mutual-info")
    {
    ApplyMetric<TPixel, VDim> adapter(this);
    adapter("MI");
    return 0;
    }

  else if(cmd == "-multiply" || cmd == "-times")
    {
    MultiplyImages<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if(cmd == "-nmi" || cmd == "-normalized-mutual-info")
    {
    ApplyMetric<TPixel, VDim> adapter(this);
    adapter("NMI");
    return 0;
    }

  else if(cmd == "-noround")
    { m_RoundFactor = 0.0; return 0; }

  // Enable SPM extensions
  else if(cmd == "-nospm")
    { m_FlagSPM = false; return 0; }

  // Overwrite / Output command - save the image without checking if
  // it already exists.
  else if(cmd == "-o")
    {
    WriteImage<TPixel, VDim> adapter(this);
    adapter(argv[1], true);
    return 1;
    }

  else if(cmd == "-origin")
    {
    RealVector org = ReadRealVector(argv[1]);
    m_ImageStack.back()->SetOrigin(org.data_block());
    return 1;
    }

  else if(cmd == "-overlap")
    {
    double label = atof(argv[1]);
    ComputeOverlaps<TPixel, VDim> adapter(this);
    adapter(label);
    return 1;
    }

  else if(cmd == "-pixel")
    {
    // Get a pixel value - no interpolation
    typename RegionType::IndexType idx = ReadIndexVector(argv[1]);
    try 
      {
      double pix = m_ImageStack.back()->GetPixel(idx);
      cout << "Pixel " << idx << " has value " << pix << endl;
      }
    catch(...)
      {
      cerr << "Error: pixel " << idx << " can not be examined!" << endl;
      }
    return 1;
    }
  
  else if(cmd == "-pop")
    {
    *verbose << "Removing (popping) the last image from the stack" << endl;
    m_ImageStack.pop_back();
    return 0;
    }

  else if(cmd == "-popas")
    {
    string var(argv[1]);
    if(m_ImageStack.size() == 0)
      throw ConvertException("No image to assign to variable %s", var.c_str());
    m_ImageVars[var] = m_ImageStack.back();
    m_ImageStack.pop_back();
    return 1;
    }

  else if(cmd == "-probe")
    {
    // Get the probe point
    RealVector x = ReadRealVector(argv[1]);
    SampleImage<TPixel, VDim> adapter(this);
    adapter(x);
    return 1;
    }

  else if (cmd == "-push" || cmd == "-get")
    {
    string var(argv[1]);
    typename ImageVariableMap::iterator img = m_ImageVars.find(var);
    if(img == m_ImageVars.end())
      throw ConvertException("No image assigned to variable %s", var.c_str());
    m_ImageStack.push_back(img->second);
    return 1;
    }


  else if(cmd == "-reciprocal")
    {
    // Multiply by the reciprocal (for time being at least)
    ReciprocalImage<TPixel, VDim>(this)();

    return 0;
    }

  else if (cmd == "-region")
    {
    // Get the position and index for the region
    RegionType bbox;
    bbox.SetIndex(ReadIndexVector(argv[1]));
    bbox.SetSize(ReadSizeVector(argv[2]));

    *verbose << "Extracting Subregion in #" << m_ImageStack.size() << endl;
    ExtractRegion<TPixel, VDim> adapter(this);
    adapter(bbox);
    return 2;
    }

  else if (cmd == "-replace")
    {
    vector<double> vReplace;

    // Read a list of numbers from the command line
    for(int i = 1; i < argc; i++)
      {
      try 
        { vReplace.push_back(myatof(argv[i])); }
      catch(...)
        { break; }
      }

    // Make sure the number of replacement rules is even
    if(vReplace.size() % 2 == 1)
      {
      cerr << "The number of parameters to '-replace' must be even!" << endl;
      throw -1;
      }

    // Replace the intensities with values supplie
    ReplaceIntensities<TPixel, VDim> adapter(this);
    adapter(vReplace);

    return vReplace.size();
    }

  // Resample command. Retain the bounding box of the image
  // while changing voxel size
  else if(cmd == "-resample")
    {
    SizeType sz = ReadSizeVector(argv[1]);
    ResampleImage<TPixel, VDim> adapter(this);
    adapter(sz);
    return 1;
    }

  else if(cmd == "-resample-mm")
    {
    RealVector vox = ReadRealVector(argv[1]);
    SizeType sz = m_ImageStack.back()->GetBufferedRegion().GetSize();
    for(size_t i = 0; i < VDim; i++)
      {
      sz[i] = (size_t)(0.5 + sz[i] * m_ImageStack.back()->GetSpacing()[i]) / vox[i];
      }
    ResampleImage<TPixel, VDim> adapter(this);
    adapter(sz);
    return 1;
    }

  else if(cmd == "-reslice-itk")
    {
    string fn_tran(argv[1]);
    ResliceImage<TPixel, VDim> adapter(this);
    adapter("itk", fn_tran);
    return 1;
    }

  else if(cmd == "-reslice-matrix")
    {
    string fn_tran(argv[1]);
    ResliceImage<TPixel, VDim> adapter(this);
    adapter("matrix", fn_tran);
    return 1;
    }

  else if(cmd == "-reslice-identity")
    {
    ResliceImage<TPixel, VDim> adapter(this);
    adapter("identity", "");
    return 0;
    }

  else if(cmd == "-rms")
    { m_AntiAliasRMS = atof(argv[1]); return 1; }
 
  else if(cmd == "-round")
    { m_RoundFactor = 0.5; return 0; }

  else if(cmd == "-scale")
    {
    double factor = atof(argv[1]);
    ScaleShiftImage<TPixel, VDim> adapter(this);
    adapter(factor, 0.0);
    return 1;
    }

  else if(cmd == "-shift")
    {
    double x = atof(argv[1]);
    ScaleShiftImage<TPixel, VDim> adapter(this);
    adapter(1.0, x);
    return 1;
    }

  else if(cmd == "-signed-distance-transform" || cmd == "-sdt")
    {
    SignedDistanceTransform<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }
  
  // Gaussian smoothing command
  else if(cmd == "-smooth")
    {
    RealVector stdev = ReadRealVector(argv[1]);
    SmoothImage<TPixel, VDim> adapter(this);
    adapter(stdev);
    return 1;
    }

  else if(cmd == "-split")
    {
    SplitMultilabelImage<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if(cmd == "-staple")
    {
    // Perform the STAPLE algorithm on the data
    double ival = atof(argv[1]);
    StapleAlgorithm<TPixel, VDim> adapter(this);
    adapter(ival);
    return 1;
    }

  // Enable SPM extensions
  else if(cmd == "-spm")
    { m_FlagSPM = true; return 0; }

  // Stretch the intensity range
  else if(cmd == "-stretch")
    {
    double u1 = ReadIntensityValue(argv[1]);
    double u2 = ReadIntensityValue(argv[2]);
    double v1 = ReadIntensityValue(argv[3]);
    double v2 = ReadIntensityValue(argv[4]);
    double a = (v2 - v1) / (u2 - u1);
    double b = v1 - u1 * a;
    ScaleShiftImage<TPixel, VDim> adapter(this);
    adapter(a, b);
    return 4;
    }

  // Thresholding
  else if(cmd == "-threshold" || cmd == "-thresh")
    {
    double u1 = strcmp(argv[1],"-inf") == 0 ? -vnl_huge_val(0.0) : ReadIntensityValue(argv[1]);
    double u2 = strcmp(argv[2],"inf") == 0 ? vnl_huge_val(0.0) : ReadIntensityValue(argv[2]);
    double v1 = ReadIntensityValue(argv[3]);
    double v2 = ReadIntensityValue(argv[4]);
    ThresholdImage<TPixel, VDim> adapter(this);
    adapter(u1, u2, v1, v2);
    return 4;
    }

  // Trim the image (trim background values from the margin)
  else if(cmd == "-trim")
    {
    // Read the size of the wrap region
    RealVector margin = ReadRealVector(argv[1]);

    // Trim the image accordingly
    TrimImage<TPixel, VDim> adapter(this);
    adapter(margin);

    // Return the number of arguments consumed
    return 1;
    }

  // Output type specification
  else if(cmd == "-type")
    { 
    m_TypeId = argv[1]; 
    int (*pf)(int) = tolower;
    transform(m_TypeId.begin(), m_TypeId.end(), m_TypeId.begin(), pf);
    return 1; 
    }

  // Verbose mode
  else if(cmd == "-verbose")
    { verbose = &std::cout; return 0; }

  else if(cmd == "-vote")
    {
    Vote<TPixel, VDim> adapter(this);
    adapter(false);
    return 0;
    }

  else if(cmd == "-vote-label")
    {
    UpdateMetadataKey<TPixel, VDim> adapter(this);
    adapter("CND:VOTE_LABEL",argv[1]);
    return 1;
    }

  else if(cmd == "-voxel-sum")
    {
    // Simply print the sum of all voxels in the image
    double sum = 0;
    size_t n = m_ImageStack.back()->GetBufferedRegion().GetNumberOfPixels();
    for(size_t i = 0; i < n; i++)
      sum += m_ImageStack.back()->GetBufferPointer()[i];
    cout << "Voxel Sum: " << sum << endl;
    return 0;
    }

  // Warp image
  else if(cmd == "-warp")
    { 
    WarpImage<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  // Warp label image
  else if(cmd == "-warp-label" || cmd=="-warplabel" || cmd=="-wl")
    { 
    RealVector stdev = ReadRealVector(argv[1]);
    WarpLabelImage<TPixel, VDim> adapter(this);
    adapter(stdev);
    return 1;
    }

  // Unknown command
  else
    { cerr << "Unknown command " << cmd << endl; throw -1; }

  cerr << "Fell through!" << endl;
  throw -1;
}


template<class TPixel, unsigned int VDim>
int
ImageConverter<TPixel, VDim>
::ProcessCommandLine(int argc, char *argv[])
{
  // Check the number of arguments
  if(argc == 1)
    {
    cerr << "PICSL convert3d tool" << endl;
    cerr << "For documentation and usage examples, see" << endl;
    cerr << "    http://alliance.seas.upenn.edu/~pauly2/wiki/index.php" 
      << "?n=Main.Convert3DTool" << endl;
    return -1;
    }

  // Try processing command line
  try 
    {
    // The last parameter in the command line is the output file
    string fnOutput = argv[argc-1];

    // Command line processing
    for(int i = 1; i < argc; ++i)
      {
      string cmd = argv[i];
      if(cmd[0] == '-')
        {
        // A command has been supplied
        i += ProcessCommand(argc-(i+1), argv+i);
        }
      else
        {
        // An image file name has been provided. If this image is followed by commands
        // read it and push in the pipeline.
        if(i != argc-1)
          {
          ReadImage<TPixel, VDim> adapter(this);
          adapter(argv[i]);
          }
        else
          {
          // Write the image, but in safe mode
          WriteImage<TPixel, VDim> adapter(this);
          adapter(argv[i], false);
          }
        }
      }
    return 0;
    }
  catch (std::exception &exc)
    {
    cerr << "Exception caught of type " << typeid(exc).name() << endl;
    cerr << "What: " << exc.what() << endl;
    return -1;
    }
  catch (...) 
    { 
    cerr << "Unknown exception caught by convert3d" << endl;
    return -1; 
    }
}

bool str_ends_with(const std::string &s, const std::string &pattern)
{
  size_t ipos = s.rfind(pattern);
  return(ipos == s.size() - pattern.size());
}

template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel, VDim>::RealVector
ImageConverter<TPixel, VDim>
::ReadRealVector(const char *vec_in)
{
  // Create a copy of the input string
  char *vec = new char[strlen(vec_in)];
  strcpy(vec, vec_in);

  size_t i;

  // Output vector
  RealVector x, scale, offset;

  // Check the type of the vector
  if(str_ends_with(vec,"mm")) 
    {
    scale.fill(1.0);
    offset.fill(0.0);
    }
  else if(str_ends_with(vec,"vox"))
    for(i = 0; i < VDim; i++)
      {
      scale[i] = m_ImageStack.back()->GetSpacing()[i];
      offset[i] = m_ImageStack.back()->GetOrigin()[i];
      }
  else 
    throw ConvertException(
      "Invalid vector spec %s (does not end with 'mm' or 'vox')", vec_in);
  
  // Find all the 'x' in the string
  char *tok = strtok(vec, "xmvo");
  for(i = 0; i < VDim && tok != NULL; i++)
    {
    x[i] = atof(tok);
    tok = strtok(NULL, "xmvo");
    } 

  // Check if there is only one number
  if(i == 1)
    x.fill(x[0]);
  else if(i != VDim)
    throw ConvertException(
      "Invalid vector spec %s (must have 1 or %i components)", vec_in, VDim);

  // Scale the vector
  for(i = 0; i < VDim; i++)
    {
    x[i] *= scale[i];
    x[i] += offset[i];
    }

  delete vec;
  return x;
}

template<class TPixel, unsigned int VDim>
TPixel
ImageConverter<TPixel, VDim>
::ReadIntensityValue(const char *vec)
{
  // Return value
  TPixel rval = 0;

  // Check if the input is infinity first
  if(!strcmp(vec, "inf") || !strcmp(vec, "+inf") || !strcmp(vec, "Inf") || !strcmp(vec, "+Inf"))
    return vnl_huge_val(0.0);
  else if(!strcmp(vec, "-inf") || !strcmp(vec, "-Inf"))
    return -vnl_huge_val(0.0);

  // Read the double part
  char *endptr;

  // Read the floating point part
  TPixel val = (TPixel) strtod(vec, &endptr);

  // Check validity
  if(endptr == vec)
    throw ConvertException("Can't convert %s to an intensity spec", vec);

  // Check if there is a '%' specification
  if(*endptr == '%')
    {
    // Check valid quantile
    if(val < 0.0 || val > 100.0)
      throw ConvertException("Invalid quantile spec %s, must be between 0 and 100", vec);

    // Compute the specified quantile
    double qtile = 0.01 * val;
    if(m_ImageStack.size() == 0)
      throw ConvertException("Can't use intensity quantile spec with no image on stack");
    size_t n = m_ImageStack.back()->GetBufferedRegion().GetNumberOfPixels();
    TPixel *asort = new TPixel[n];
    std::copy(m_ImageStack.back()->GetBufferPointer(), m_ImageStack.back()->GetBufferPointer() + n, asort);
    std::sort(asort, asort+n);
    size_t k = (size_t) (qtile * n);
    val = asort[k];
    delete asort;
    }

  return val;
}

template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel, VDim>::SizeType
ImageConverter<TPixel, VDim>
::ReadSizeVector(const char *vec_in)
{
  // Create a copy of the input string
  char *vec = new char[strlen(vec_in)];
  strcpy(vec, vec_in);

  size_t i;

  typename ImageType::SizeType sz;

  // Check if the string ends with %
  if(str_ends_with(vec, "%"))
    {
    // Read floating point size
    RealVector factor;
    char *tok = strtok(vec, "x%");
    for(i = 0; i < VDim && tok != NULL; i++)
      {
      factor[i] = atof(tok);
      if(factor[i] <= 0)
        throw ConvertException("Non-positive percent size specification: %s", vec_in);
      tok = strtok(NULL, "x%");
      } 

    if(i == 1)
      factor.fill(factor[0]);

    // Get the size of the image in voxels
    for(size_t i = 0; i < VDim; i++)
      sz[i] = (unsigned long)(m_ImageStack.back()->GetBufferedRegion().GetSize(i) * 0.01 * factor[i] + 0.5);
    }
  else
    {
    // Find all the 'x' in the string
    char *tok = strtok(vec, "x");
    for(size_t i = 0; i < VDim; i++)
      {
      if(tok == NULL)
        throw ConvertException("Invalid size specification: %s", vec_in);
      int x = atoi(tok);
      if(x <= 0)
        throw ConvertException("Non-positive size specification: %s", vec_in);
      sz[i] = (unsigned long)(x);
      tok = strtok(NULL, "x");
      } 
    }
  
  delete vec;
  return sz;
}


template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel, VDim>::IndexType
ImageConverter<TPixel, VDim>
::ReadIndexVector(const char *vec_in)
{
  // Create a copy of the input string
  char *vec = new char[strlen(vec_in)];
  strcpy(vec, vec_in);

  size_t i;

  typename ImageType::IndexType idx;

  // Check if the string ends with %
  if(str_ends_with(vec, "%"))
    {
    // Read floating point size
    RealVector factor;
    char *tok = strtok(vec, "x%");
    for(i = 0; i < VDim && tok != NULL; i++)
      {
      factor[i] = atof(tok);
      tok = strtok(NULL, "x%");
      } 

    if(i == 1)
      factor.fill(factor[0]);

    // Get the size of the image in voxels
    for(size_t i = 0; i < VDim; i++)
      idx[i] = (long)(m_ImageStack.back()->GetBufferedRegion().GetSize(i) * 0.01 * factor[i] + 0.5);
    }
  else
    {
    // Find all the 'x' in the string
    char *tok = strtok(vec, "x");
    for(size_t i = 0; i < VDim; i++)
      {
      if(tok == NULL)
        throw ConvertException("Invalid index specification: %s", vec_in);
      int x = atoi(tok);
      idx[i] = (long)(x);
      tok = strtok(NULL, "x");
      } 
    }

  delete vec;
  return idx;
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::CopyImage()
{
  // Get the input image
  ImagePointer input = m_ImageStack.back();
  
  // Simply make a copy of the input image on the stack
  ImagePointer output = ImageType::New();
  output->SetRegions(input->GetBufferedRegion());
  output->SetSpacing(input->GetSpacing());
  output->SetOrigin(input->GetOrigin());
  output->SetDirection(input->GetDirection());
  output->SetMetaDataDictionary(input->GetMetaDataDictionary());
  output->Allocate();

  // Copy everything
  size_t n = input->GetBufferedRegion().GetNumberOfPixels();
  for(size_t i = 0; i < n; i++)
    output->GetBufferPointer()[i] = input->GetBufferPointer()[i];

  // Put on the end of the stack
  m_ImageStack.pop_back();
  m_ImageStack.push_back(output);
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::GetBoundingBox(ImageType *image, RealVector &bb0, RealVector &bb1)
{
  for(size_t i = 0; i < VDim; i++)
    {
    bb0[i] = image->GetOrigin()[i];
    bb1[i] = bb0[i] + image->GetSpacing()[i] * image->GetBufferedRegion().GetSize()[i];
    }
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::PrintMatrix(std::ostream &sout, vnl_matrix<double> mat, const char *fmt, const char *prefix)
{
  // Print each row and column of the matrix
  char buffer[256];
  for(size_t i = 0; i < mat.rows(); i++)
    {
    sout << prefix;
    for(size_t j = 0; j < mat.columns(); j++)
      {
      sprintf(buffer, fmt, mat(i,j));
      sout << buffer;
      }
    sout << endl;
    }
}

template<class TPixel, unsigned int VDim>
size_t
ImageConverter<TPixel, VDim>
::ForEachLoop(int argc, char *argv[])
{
  size_t narg;

  // Note: this is a rather lame implementation that uses recursion with
  // a state variable to repeat a range of command a bunch of times

  // Back up the current stack
  vector<ImagePointer> in_stack = m_ImageStack, out_stack;

  // Print out what's going on
  *verbose << "Repeating commands for all " << in_stack.size() << " images" << endl;

  // Loop over all images
  for(size_t i = 0; i < in_stack.size(); i++)
    {
    narg = 1;

    // Set up the image stack
    m_ImageStack.clear();
    m_ImageStack.push_back(in_stack[i]);

    // Set the in-loop flag
    m_InLoop = true;
    
    // When the -endfor is encountered, the InLoop flag will be switched
    while(m_InLoop)
      narg += 1 + this->ProcessCommand(argc-narg, argv+narg);

    // Place the result (if any) on the output stack
    if(m_ImageStack.size() > 1)
      throw ConvertException("Commands in the -foreach clause may not produce multiple outputs");
    else if(m_ImageStack.size() == 1)
      out_stack.push_back(m_ImageStack.back());
    }

  // Update the stack 
  m_ImageStack = out_stack;

  // Return the number of arguments to the next command
  return narg - 1;
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::SetInterpolator()
{
  // Get the interpolation type
  string terp = m_Interpolation.c_str();
  std::transform(terp.begin(), terp.end(), terp.begin(), tolower);

  // Create the interpolator depending on parameter
  if(terp == "nearestneighbor" || terp == "nn" || terp == "0")
    {
    typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double> NNInterpolatorType;
    m_Interpolator = NNInterpolatorType::New();
    }
  else if(terp == "linear" || terp == "1")
    {
    typedef itk::LinearInterpolateImageFunction<ImageType,double> LinearInterpolatorType;
    m_Interpolator = LinearInterpolatorType::New();
    }
  else if(terp == "cubic" || terp == "3")
    {
    typedef itk::BSplineInterpolateImageFunction<ImageType,double> CubicInterpolatorType;
    m_Interpolator = CubicInterpolatorType::New();
    }
  else if(terp == "sinc")
    {
    typedef itk::WindowedSincInterpolateImageFunction<ImageType, 4> SincInterpolatorType;
    m_Interpolator = SincInterpolatorType::New();
    }
  else
    { 
    throw ConvertException("Unknown interpolation type: %s", m_Interpolation.c_str());
    }
}


template class ImageConverter<double, 2>;
template class ImageConverter<double, 3>;
 
