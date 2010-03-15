#include "ConvertImageND.h"

#include "AddImages.h"
#include "AntiAliasImage.h"
#include "ApplyMetric.h"
#include "BiasFieldCorrection.h"
#include "BinaryImageCentroid.h"
#include "ClipImageIntensity.h"
#include "ComputeFFT.h"
#include "ComputeOverlaps.h"
#include "ConnectedComponents.h"
#include "ConvertAdapter.h"
#include "CopyTransform.h"
#include "CreateImage.h"
#include "CreateInterpolator.h"
#include "ExtractRegion.h"
#include "ExtractSlice.h"
#include "GeneralLinearModel.h"
#include "HistogramMatch.h"
#include "ImageERF.h"
#include "ImageLaplacian.h"
#include "LabelOverlapMeasures.h"
#include "LabelStatistics.h"
#include "LevelSetSegmentation.h"
#include "MathematicalMorphology.h"
#include "MixtureModel.h"
//#include "MRFVote.h"
#include "MultiplyImages.h"
#include "NormalizedCrossCorrelation.h"
#include "PadImage.h"
#include "PeronaMalik.h"
#include "PrintImageInfo.h"
#include "Rank.h"
#include "ReadImage.h"
#include "ReciprocalImage.h"
#include "ReplaceIntensities.h"
#include "ResampleImage.h"
#include "ResliceImage.h"
#include "SampleImage.h"
#include "ScaleShiftImage.h"
#include "SetSform.h"
#include "SignedDistanceTransform.h"
#include "SmoothImage.h"
#include "SplitMultilabelImage.h"
#include "StapleAlgorithm.h"
#include "ThresholdImage.h"
#include "TrimImage.h"
#include "UnaryMathOperation.h"
#include "UpdateMetadataKey.h"
#include "Vote.h"
#include "VoxelwiseRegression.h"
#include "WarpImage.h"
#include "WarpLabelImage.h"
#include "WriteImage.h"
#include "WeightedSum.h"

#include <cstring>
#include <algorithm>

// Support for regular expressions via KWSYS in ITK
#include <itksys/RegularExpression.hxx>

using namespace itksys;

extern const char *ImageConverter_VERSION_STRING;

// Helper function: read a double, throw exception if unreadable
double myatof(char *str)
{
  char *end = 0;
  double d = strtod(str, &end);
  if(*end != 0)
    throw "strtod conversion failed";
  return d;
};

std::string str_to_lower(const char *input)
{
  std::string s(input);
  std::transform(s.begin(), s.end(), s.begin(), (int(*)(int)) tolower);
  return s;
}


template<class TPixel, unsigned int VDim>
ImageConverter<TPixel,VDim>
::ImageConverter()
  : verbose(&devnull)
{
  // Initialize to defaults
  m_TypeId = "float";
  m_Background = 0.0;
  m_RoundFactor = 0.5;
  m_FlagSPM = false;
  m_MultiComponentSplit = false;
  m_AntiAliasRMS = 0.07;
  m_Iterations = 0;
  m_InLoop = false;
  m_PercentIntensityMode = PIM_QUANTILE;

  m_LevSetCurvature = 0.2;
  m_LevSetAdvection = 0.0;

  // Create an interpolator
  m_Interpolation = "linear";
  CreateInterpolator<TPixel, VDim>(this).CreateLinear();
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::PrintCommandListing()
{
  cout << "Command Listing: " << endl;
  cout << "    -add" << endl;
  cout << "    -anisotropic-diffusion, -ad" << endl;
  cout << "    -antialias, -alias" << endl;
  cout << "    -as, -set" << endl;
  cout << "    -background" << endl;
  cout << "    -biascorr" << endl;
  cout << "    -binarize" << endl;
  cout << "    -centroid" << endl;
  cout << "    -connected-components, -connected, -comp" << endl;
  cout << "    -clear" << endl;
  cout << "    -clip" << endl;
  cout << "    -copy-transform, -ct" << endl;
  cout << "    -create" << endl;
  cout << "    -dilate" << endl;
  cout << "    -divide" << endl;
  cout << "    -dup" << endl;
  cout << "    -endfor" << endl;
  cout << "    -erode" << endl;
  cout << "    -erf" << endl;
  cout << "    -exp" << endl;
  cout << "    -fft" << endl;
  cout << "    -foreach" << endl;
  cout << "    -glm" << endl;
  cout << "    -histmatch, -histogram-match" << endl;
  cout << "    -info" << endl;
  cout << "    -info-full" << endl;
  cout << "    -insert, -ins" << endl;
  cout << "    -interpolation, -interp, -int" << endl;
  cout << "    -iterations" << endl;
  cout << "    -label-overlap" << std::endl;
  cout << "    -label-statistics, -lstat" << endl;
  cout << "    -laplacian, -laplace" << endl;
  cout << "    -levelset" << endl;
  cout << "    -levelset-curvature" << endl;
  cout << "    -levelset-advection" << endl;
  cout << "    -ln, -log" << endl;
  cout << "    -log10" << endl;
  cout << "    -mcs, -multicomponent-split" << endl;
  cout << "    -mean" << endl;
  cout << "    -merge" << endl;
  cout << "    -mi, -mutual-info" << endl;
  cout << "    -mixture, -mixture-model" << endl;
  cout << "    -multiply, -times" << endl;
  cout << "    -ncc, -normalized-cross-correlation" << endl;
  cout << "    -nmi, -normalized-mutual-info" << endl;
  cout << "    -normpdf" << endl;
  cout << "    -noround" << endl;
  cout << "    -nospm" << endl;
  cout << "    -o" << endl;
  cout << "    -omc, -output-multicomponent" << endl;
  cout << "    -oo, -output-multiple" << endl;
  cout << "    -origin" << endl;
  cout << "    -overlap" << endl;
  cout << "    -pad" << endl;
  cout << "    -percent-intensity-mode, -pim" << endl;
  cout << "    -pixel" << endl;
  cout << "    -pop" << endl;
  cout << "    -popas" << endl;
  cout << "    -probe" << endl;
  cout << "    -push, -get" << endl;
  cout << "    -rank" << endl;
  cout << "    -reciprocal" << endl;
  cout << "    -region" << endl;
  cout << "    -replace" << endl;
  cout << "    -resample" << endl;
  cout << "    -resample-mm" << endl;
  cout << "    -reslice-itk" << endl;
  cout << "    -reslice-matrix" << endl;
  cout << "    -reslice-identity" << endl;
  cout << "    -rms" << endl;
  cout << "    -round" << endl;
  cout << "    -scale" << endl;
  cout << "    -set-sform" << endl;
  cout << "    -shift" << endl;
  cout << "    -signed-distance-transform, -sdt" << endl;
  cout << "    -slice" << endl;
  cout << "    -smooth" << endl;
  cout << "    -spacing" << endl;
  cout << "    -split" << endl;
  cout << "    -sqrt" << endl;
  cout << "    -staple" << endl;
  cout << "    -spm" << endl;
  cout << "    -stretch" << endl;
  cout << "    -threshold, -thresh" << endl;
  cout << "    -trim" << endl;
  cout << "    -trim-to-size" << endl;
  cout << "    -type" << endl;
  cout << "    -verbose" << endl;
  cout << "    -version" << endl;
  cout << "    -vote" << endl;
  cout << "    -vote-label" << endl;
  cout << "    -voxel-sum" << endl;
  cout << "    -voxel-integral, -voxel-int" << endl;
  cout << "    -voxelwise-regression, -voxreg" << endl;
  cout << "    -warp" << endl;
  cout << "    -weighted-sum, -wsum" << endl;
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
    {
    m_Background = atof(argv[1]);
    *verbose << "Background value set to " << m_Background << endl;
    return 1;
    }

  else if(cmd == "-biascorr")
    {
    BiasFieldCorrection<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  // f(x) = (x == xBackground) ? 0 : 1
  else if(cmd == "-binarize")
    {
    ThresholdImage<TPixel, VDim> adapter(this);
    adapter(m_Background, m_Background, 0.0, 1.0);
    return 0;
    }

  else if(cmd == "-centroid")
    {
    BinaryImageCentroid<TPixel, VDim> adapter(this);
    adapter();
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
    RealVector voxel = ReadRealSize(argv[2]);
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

  else if(cmd == "-dup" || cmd == "-duplicate")
    {
    m_ImageStack.push_back(m_ImageStack.back());
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

  else if(cmd == "-exp")
    {
    UnaryMathOperation<TPixel, VDim> adapter(this);
    adapter(&vcl_exp);
    return 0;
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

  else if (cmd == "-h")
    {
    PrintCommandListing();
    return 0;
    }


  else if(cmd == "-histmatch" || cmd == "-histogram-match")
    {
    size_t nmatch = atoi(argv[1]);
    HistogramMatch<TPixel, VDim> adapter(this);
    adapter(nmatch);
    return 1;
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
    // Adapter that creates interpolators
    CreateInterpolator<TPixel, VDim> adapter(this);

    // Interpret the interpolation type
    m_Interpolation = str_to_lower(argv[1]);

    if(m_Interpolation == "nearestneighbor" || m_Interpolation == "nn" || m_Interpolation == "0")
      {
      adapter.CreateNN();
      return 1;
      }
    else if(m_Interpolation == "linear" || m_Interpolation == "1")
      {
      adapter.CreateLinear();
      return 1;
      }
    else if(m_Interpolation == "cubic" || m_Interpolation == "3")
      {
      adapter.CreateCubic();
      return 1;
      }
    else if(m_Interpolation == "sinc")
      {
      adapter.CreateSinc();
      return 1;
      }
    else if(m_Interpolation == "gaussian")
      {
      RealVector sigma = ReadRealSize(argv[2]);
      adapter.CreateGaussian(sigma);
      return 2;
      }
    else
      {
      throw ConvertException("Unknown interpolation type: %s", m_Interpolation.c_str());
      }
    }

  else if(cmd == "-iterations")
    {
    m_Iterations = static_cast<size_t>(atoi(argv[1]));
    return 1;
    }

  else if(cmd == "-label-overlap")
    {
    LabelOverlapMeasures<TPixel, VDim>(this)();
    return 0;
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

  else if(cmd == "-ln" || cmd == "-log")
    {
    UnaryMathOperation<TPixel, VDim> adapter(this);
    adapter(&vcl_log);
    return 0;
    }

  else if(cmd == "-log10")
    {
    UnaryMathOperation<TPixel, VDim> adapter(this);
    adapter(&vcl_log10);
    return 0;
    }

  else if(cmd == "-mcs" || cmd == "-multicomponent-split")
    {
    m_MultiComponentSplit = true; return 0;
    }

  else if(cmd == "-mean")
    {
    size_t n = m_ImageStack.size();
    for(size_t i = 1; i < n; i++)
      {
      AddImages<TPixel,VDim> adapter(this);
      adapter();
      }
    ScaleShiftImage<TPixel, VDim> scaler(this);
    scaler(1.0 / n, 0.0);
    return 0;
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
    int nret = 0;
    string fn("none");
    if (argc == 2)
      {
      fn = argv[1];
      nret = 1;
      }
    adapter("MI", fn.c_str());
    return nret;
    }

  else if(cmd == "-mixture" || cmd == "-mixture-model")
    {
    int ncomp = atoi(argv[1]);
    if(ncomp == 0)
      throw ConvertException("Incorrect specification of mixture model initialization");

    std::vector<double> mu, sigma;
    for(int i = 0; i < ncomp; i++)
      {
      mu.push_back(ReadIntensityValue(argv[2 + 2 * i]));
      sigma.push_back(ReadIntensityValue(argv[3 + 2 * i]));
      }

    MixtureModel<TPixel, VDim> adapter(this);
    adapter(mu, sigma);

    return 1 + 2 * ncomp;
    }

  else if(cmd == "-mmi" || cmd == "-mattes-mutual-info")
    {
    ApplyMetric<TPixel, VDim> adapter(this);
    int nret = 0;
    string fn("none");
    if (argc == 2)
      {
      fn = argv[1];
      nret = 1;
      }
    adapter("MMI", fn.c_str());
    return nret;
    }
  else if(cmd == "-msq" || cmd == "-mean-square")
    {
    ApplyMetric<TPixel, VDim> adapter(this);
    int nret = 0;
    string fn("none");
    if (argc == 2)
      {
      fn = argv[1];
      nret = 1;
      }
    adapter("MSQ", fn.c_str());
    return nret;
    }
  else if(cmd == "-multiply" || cmd == "-times")
    {
    MultiplyImages<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if(cmd == "-ncc" || cmd == "-normalized-cross-correlation")
    {
    SizeType sz = ReadSizeVector(argv[1]);
    NormalizedCrossCorrelation<TPixel,VDim> adapter(this);
    adapter(sz);
    return 1;
    }

  else if(cmd == "-ncor" || cmd == "-normalized-correlation")
    {
    ApplyMetric<TPixel, VDim> adapter(this);
    int nret = 0;
    string fn("none");
    if (argc == 2)
      {
      fn = argv[1];
      nret = 1;
      }
    adapter("NCOR", fn.c_str());
    return nret;
    }
  else if(cmd == "-nmi" || cmd == "-normalized-mutual-info")
    {
    ApplyMetric<TPixel, VDim> adapter(this);
    int nret = 0;
    string fn("none");
    if (argc == 2)
      {
      fn = argv[1];
      nret = 1;
      }
    adapter("NMI", fn.c_str());
    return nret;
    }

  else if(cmd == "-nomcs" || cmd == "-no-multicomponent-split")
    {
    m_MultiComponentSplit = false; return 0;
    }

  else if(cmd == "-normpdf")
    {
    // Compute normal PDF of intensity values given sigma and mu
    double mu = atof(argv[1]);
    double s = atof(argv[2]);

    // Subtract mu
    ScaleShiftImage<TPixel, VDim> scale1(this);
    scale1(1.0, -mu);

    // Square
    m_ImageStack.push_back(m_ImageStack.back());
    MultiplyImages<TPixel, VDim> times(this);
    times();

    // Scale by -1/2s
    ScaleShiftImage<TPixel, VDim> scale2(this);
    scale2(-0.5 / s, 0.0);

    // Exponentiate
    UnaryMathOperation<TPixel, VDim> exp1(this);
    exp1(&vcl_exp);

    // Scale by factor
    ScaleShiftImage<TPixel, VDim> scale3(this);
	scale3(1.0 / sqrt(2 * vnl_math::pi * s * s), 0.0);
    }

  else if(cmd == "-noround")
    { m_RoundFactor = 0.0; return 0; }

  // Enable SPM extensions
  else if(cmd == "-nospm")
    { m_FlagSPM = false; return 0; }

  // Overwrite / Output command - save the image without checking if
  // it already exists.
  else if(cmd == "-o" || cmd == "-output")
    {
    WriteImage<TPixel, VDim> adapter(this);
    adapter(argv[1], true);
    return 1;
    }

  else if(cmd == "-omc" || cmd == "-output-multicomponent")
    {
    // Number of components (all by default)
    int nc, np;

    // A parameter can be optionally specified (how many components)
    RegularExpression re("[0-9]+");
    if(re.find(argv[1]))
      { nc = atoi(argv[1]); np=2; }
    else
      { nc = m_ImageStack.size(); np = 1; }

    // Create a writer
    WriteImage<TPixel, VDim> adapter(this);
    adapter.WriteMultiComponent(argv[np], nc);
    return np;
    }

  // Write mulptiple images
  else if(cmd == "-oo" || cmd == "-output-multiple")
    {
    // Check if the argument is a printf pattern
    char buffer[1024];
    sprintf(buffer, argv[1],0);
    if(strcmp(buffer, argv[1]))
      {
      // A pattern is specified. For each image on the stack, use pattern
      for(size_t i = 0; i < m_ImageStack.size(); i++)
        {
        sprintf(buffer, argv[1], i);
        WriteImage<TPixel, VDim> adapter(this);
        adapter(buffer, true, i);
        }
      return 1;
      }
    else
      {
      // Filenames are specified. Find out how many there are
      size_t nfiles = 0;
      for(int i = 1; i < argc; i++)
        {
        if(argv[i][0] != '-') nfiles++; else break;
        }


      if(nfiles == 0)
        throw ConvertException("No files specified to -oo command");

      if(nfiles > m_ImageStack.size())
        throw ConvertException("Too many files specified to -oo command");

      for(size_t j = 0; j < nfiles; j++)
        {
        WriteImage<TPixel, VDim> adapter(this);
        adapter(argv[j+1], true, m_ImageStack.size() - nfiles + j);
        }

      return nfiles;
      }
    }

  else if(cmd == "-orient")
    {
    // Read an RAS code
    RegularExpression re("[raslpi]{3}");
    if(re.find(str_to_lower(argv[1])))
      { cout << "You supplied a RAS code" << endl; }
    else
      { cout << "I am expecting a matrix" << endl; }
    return 1;
    }


  else if(cmd == "-origin")
    {
    RealVector org = ReadRealVector(argv[1], true);
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

  else if(cmd == "-pad")
    {
    // specify two sizes that give the padding in x,y,z
    // pad is the offset (in voxels) from the edge of the image to the
    // padExtent. Example: -pad 1x1x1vox 0x0x0vox pads on the left, posterior, inferior sides
    // by one voxel -pad 10x10x10% 10x10x10% adds 10% on all sides
    IndexType padExtentLower, padExtentUpper;

    padExtentLower = ReadIndexVector(argv[1]);
    padExtentUpper = ReadIndexVector(argv[2]);

    float padValue = atof(argv[3]);

    *verbose << "Padding image #" << m_ImageStack.size() << endl;

    PadImage<TPixel, VDim> adapter(this);
    adapter(padExtentLower, padExtentUpper, padValue);
    return 3;
    }

  else if(cmd == "-percent-intensity-mode" || cmd == "-pim")
    {
    // What does % mean when specifying intensities
    string pim = str_to_lower(argv[1]);
    if(pim == "quantile" || pim == "q")
      m_PercentIntensityMode = PIM_QUANTILE;
    else if(pim == "foregroundquantile" || pim == "fq")
      m_PercentIntensityMode = PIM_FGQUANTILE;
    else if(pim == "range" || pim == "r")
      m_PercentIntensityMode = PIM_RANGE;
    else
      throw ConvertException("Wrong -percent-intensity-mode spec %s. See help.", pim.c_str());
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
    RealVector x = ReadRealVector(argv[1], true);
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

  else if(cmd == "-rank")
    {
    Rank<TPixel,VDim> adapter(this);
    adapter();
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


  else if(cmd == "-reciprocal")
    {
    // Multiply by the reciprocal (for time being at least)
    ReciprocalImage<TPixel, VDim>(this)();

    return 0;
    }

  else if (cmd == "-set-sform")
    {
    string fn_tran( argv[1] );

    *verbose << "Setting sform of image #" << m_ImageStack.size() << endl;
    SetSform<TPixel, VDim> adapter(this);
    adapter( fn_tran );

    return 1;
    }

  else if (cmd == "-slice")
    {
    string axis( argv[1] );
    int pos = atoi(argv[2]);
    char * filename = argv[3];

    *verbose << "Extracting slice in #" << m_ImageStack.size() << endl;
    ExtractSlice<TPixel, VDim> adapter(this);
    adapter(axis, pos, filename);

    return 3;
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
    RealVector vox = ReadRealSize(argv[1]);
    SizeType sz = m_ImageStack.back()->GetBufferedRegion().GetSize();
    for(size_t i = 0; i < VDim; i++)
      {
      sz[i] = static_cast<size_t>((0.5 + sz[i] * m_ImageStack.back()->GetSpacing()[i]) / vox[i]);
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
    RealVector stdev = ReadRealSize(argv[1]);
    SmoothImage<TPixel, VDim> adapter(this);
    adapter(stdev);
    return 1;
    }

  else if(cmd == "-spacing")
    {
    RealVector org = ReadRealSize(argv[1]);
    m_ImageStack.back()->SetSpacing(org.data_block());
    return 1;
    }

  else if(cmd == "-split")
    {
    SplitMultilabelImage<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if(cmd == "-sqrt")
    {
    UnaryMathOperation<TPixel, VDim> adapter(this);
    adapter(&vcl_sqrt);
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
    RealVector margin = ReadRealSize(argv[1]);

    // Trim the image accordingly
    TrimImage<TPixel, VDim> adapter(this);
    adapter(margin, TrimImage<TPixel, VDim>::SPECIFY_MARGIN);

    // Return the number of arguments consumed
    return 1;
    }

  else if(cmd == "-trim-to-size")
    {
    // Read the size of the trim region
    RealVector size = ReadRealSize(argv[1]);

    // Trim the image accordingly
    TrimImage<TPixel, VDim> adapter(this);
    adapter(size, TrimImage<TPixel, VDim>::SPECIFY_FINALSIZE);

    // Return the number of arguments consumed
    return 1;
    }

  // Output type specification
  else if(cmd == "-type")
    {
    m_TypeId = str_to_lower(argv[1]);
    return 1;
    }

  // Verbose mode
  else if(cmd == "-verbose")
    { verbose = &std::cout; return 0; }

  else if(cmd == "-version")
    {
    cout << "Version " << ImageConverter_VERSION_STRING << endl;
    return 0;
    }

  else if(cmd == "-vote")
    {
    Vote<TPixel, VDim> adapter(this);
    adapter(false);
    return 0;
    }
//
//  else if(cmd == "-vote-mrf")
//    {
//    double beta = atof(argv[1]);
//    size_t iter = atoi(argv[2]);
//    MRFVote<TPixel, VDim> adapter(this);
//    adapter(beta, iter, false);
//    return 2;
//    }
//
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

  else if(cmd == "-voxel-integral" || cmd == "-voxel-int")
    {
    // Simply print the sum of all voxels in the image
    double sum = 0;
    size_t n = m_ImageStack.back()->GetBufferedRegion().GetNumberOfPixels();
    for(size_t i = 0; i < n; i++)
      sum += m_ImageStack.back()->GetBufferPointer()[i];
    double vol = 1.0;
    for(size_t d = 0; d < VDim; d++)
      vol *= m_ImageStack.back()->GetSpacing()[d];
    cout << "Voxel Integral: " << sum * vol << endl;
    return 0;
    }

  else if(cmd == "-voxelwise-regression" || cmd == "-voxreg")
    {
    // Get the order
    size_t order = atoi(argv[1]);
    VoxelwiseRegression<TPixel, VDim> adapter(this);
    adapter(order);
    return 1;
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
    RealVector stdev = ReadRealSize(argv[1]);
    WarpLabelImage<TPixel, VDim> adapter(this);
    adapter(stdev);
    return 1;
    }

  else if(cmd == "-weighted-sum" || cmd == "-wsum")
    {
    std::vector<double> weights;
    for(int i = 1; i < argc; i++)
      if(argv[i][0] != '-')
        weights.push_back(atof(argv[i]));
      else break;
    WeightedSum<TPixel,VDim> adapter(this);
    adapter(weights);
    return weights.size();
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
  // Disable multithreading
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);

  // Check the number of arguments
  if(argc == 1)
    {
    cerr << "PICSL convert3d tool" << endl;
    cerr << "For full documentation and usage examples, see" << endl;
    cerr << "    http://www.itksnap.org/pmwiki/pmwiki.php?n=Convert3D.Documentation" << endl;
    cerr << "For a brief summary of commands, call" << endl;
    cerr << "    " << argv[0] << " -h" << endl;
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
        i += ProcessCommand(argc-i, argv+i);
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

// How the specification is made
enum VecSpec { PHYSICAL, VOXELS, PERCENT, NONE };

template<unsigned int VDim>
void ReadVecSpec(const char *vec_in, vnl_vector_fixed<double,VDim> &vout, VecSpec &type)
{
  // Set up the regular expressions for numerical string parsing
  RegularExpression re1(
    "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)(mm|vox|%)?");
  RegularExpression re2(
    "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)x"
    "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)(mm|vox|%)?");
  RegularExpression re3(
    "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)x"
    "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)x"
    "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)(mm|vox|%)?");

  // Lowercase string
  string vec = str_to_lower(vec_in);
  string spec;

  // Check if it's a single-component specification
  if(VDim == 2 && re2.find(vec))
    {
    vout[0] = atof(re2.match(1).c_str());
    vout[1] = atof(re2.match(3).c_str());
    spec = re2.match(5);
    }
  else if(VDim == 3 && re3.find(vec))
    {
    vout[0] = atof(re3.match(1).c_str());
    vout[1] = atof(re3.match(3).c_str());
    vout[2] = atof(re3.match(5).c_str());
    spec = re3.match(7);
    }
  else if(re1.find(vec))
    {
    vout.fill(atof(re1.match(1).c_str()));
    spec = re1.match(3);
    }
  else throw ConvertException("Invalid vector specification %s", vec_in);

  // Get the type of spec. Luckily, all suffixes have the same length
  switch(spec.length()) {
  case 0: type = NONE; break;
  case 1: type = PERCENT; break;
  case 2: type = PHYSICAL; break;
  case 3: type = VOXELS; break;
  default: throw ConvertException("Internal error in VecSpec code");
  }
}

template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel, VDim>::RealVector
ImageConverter<TPixel, VDim>
::ReadRealVector(const char *vec_in, bool is_point)
{
  // Output vector
  RealVector x, scale, offset;
  VecSpec type;

  // Read the vector
  ReadVecSpec(vec_in, x, type);

  // Check the type of the vector
  if(type != VOXELS && type != PHYSICAL)
    throw ConvertException(
      "Invalid vector spec %s (must end with 'mm' or 'vox')", vec_in);

  // If the vector is in vox units, map it to physical units
  if(type == VOXELS)
    {
    // Get the matrix
    typename ImageType::TransformMatrixType MP =
      m_ImageStack.back()->GetVoxelSpaceToRASPhysicalSpaceMatrix();

    // Create the vector to multiply by
    vnl_vector_fixed<double, VDim+1> X, XP;
    for(size_t d = 0; d < VDim; d++)
      X[d] = x[d];
    X[VDim] = is_point ? 1.0 : 0.0;

    // Apply matrix
    XP = MP * X;
    for(size_t d = 0; d < VDim; d++)
      x[d] = XP[d];
    }

  return x;
}

template<class TPixel, unsigned int VDim>
TPixel
ImageConverter<TPixel, VDim>
::ReadIntensityValue(const char *vec)
{
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
    if(m_PercentIntensityMode == PIM_QUANTILE || m_PercentIntensityMode == PIM_FGQUANTILE)
      {
      // Check valid quantile
      if(val < 0.0 || val > 100.0)
        throw ConvertException("Invalid quantile spec %s, must be between 0 and 100", vec);

      // Compute the specified quantile
      double qtile = 0.01 * val;
      if(m_ImageStack.size() == 0)
        throw ConvertException("Can't use intensity quantile spec with no image on stack");

      // Allocate an array for sorting
      size_t n = m_ImageStack.back()->GetBufferedRegion().GetNumberOfPixels();
      TPixel *asort = new TPixel[n], *p = asort, *q = m_ImageStack.back()->GetBufferPointer();

      // Copy intensity values that are legit, ignore nans
      for(size_t i = 0; i < n; i++, q++)
        {
        // We don't include nans and if FGQUANTILE, background values
        if(!vnl_math_isnan(*q))
          if(m_PercentIntensityMode == PIM_QUANTILE || *q != m_Background)
            {*p = *q; ++p;}
        }

      // Get the size of the sort array
      size_t np = p - asort;
      if(np == 0)
        {
        if(m_PercentIntensityMode == PIM_QUANTILE)
          throw ConvertException("Quantile could not be computed because the image has only NANs");
        else
          throw ConvertException("Foreground quantile could not be computed because the image has only background");
        }

      // Sort the acceptable values
      std::sort(asort, p);

      // Get the quantile
      size_t k = (size_t) (qtile * np);
      val = asort[k];
      delete asort;

      if(m_PercentIntensityMode == PIM_QUANTILE)
        *verbose << "Quantile " << qtile << " maps to " << val << endl;
      else
        *verbose << "Foreground quantile " << qtile <<  " (over "
        << np << " voxels) maps to " << val << endl;
      }
    else
      {
      // A range specification. We need the min and max of the image
      double imin = numeric_limits<double>::max();
      double imax = - numeric_limits<double>::max();
      size_t n = m_ImageStack.back()->GetBufferedRegion().GetNumberOfPixels();
      TPixel *q = m_ImageStack.back()->GetBufferPointer();
      for(size_t i = 0; i < n; i++, q++)
        {
        if(*q < imin) imin = *q;
        if(*q > imax) imax = *q;
        }

      double rspec = val * 0.01;
      val = imin + rspec * (imax - imin);
      *verbose << "Intensity range spec " << rspec << " maps to " << val << endl;
      }
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
  size_t narg = 0;

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

template class ImageConverter<double, 2>;
template class ImageConverter<double, 3>;

