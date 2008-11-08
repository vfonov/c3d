#include "itkOrientedRASImage.h"
#include "itkImageFileReader.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_det.h"
#include "vnl/vnl_inverse.h"
#include <iostream>

using namespace std;

// Typedefs
typedef itk::OrientedRASImage<double, 3> ImageType;
typedef vnl_matrix_fixed<double, 4, 4> MatrixType;
typedef std::vector<MatrixType> MatrixStack;

int usage()
{
  cout << "RAS Affine Transform Tool" << endl;
  cout << "Usage: " << endl;
  cout << "  c3d_affine_tool [transform_files | options] " << endl;
  cout << "Options: " << endl;
  cout << "  -ref image    Set reference (fixed) image" << endl;
  cout << "  -src image    Set source (moving) image" << endl;
  cout << "  -o matfile    Write output matrix" << endl;
  cout << "  -fsl2ras      Convert FSL to RAS" << endl;
  cout << "  -mult         Multiply matrices" << endl;
  cout << "  -inv          Invert matrix" << endl;
  return -1;
}

void ras_read(MatrixStack &vmat, const char *fname)
{
  MatrixType mat;

  ifstream fin(fname);
  for(size_t i = 0; i < 4; i++)
    for(size_t j = 0; j < 4; j++)
      if(fin.good())
        {
        fin >> mat[i][j];
        }
      else
        {
        throw "Unable to read matrix";
        }
  fin.close();

  vmat.push_back(mat);
}

void ras_write(MatrixStack &vmat, const char *fname)
{
  MatrixType mat = vmat.back();

  ofstream fout(fname);
  for(size_t i = 0; i < 4; i++)
    for(size_t j = 0; j < 4; j++)
      fout << mat[i][j] << (j < 3 ? " " : "\n");

  fout.close();

  vmat.pop_back();
}

void fsl_to_ras(MatrixStack &vmat, ImageType *ref, ImageType *mov)
{
  MatrixType m_fsl, m_spcref, m_spcmov, m_swpref, m_swpmov, m_ref, m_mov, m_out;
  m_fsl = vmat.back();

  // Set the ref/mov matrices
  m_ref = ref->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix();
  m_mov = mov->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix();

  // Set the swap matrices
  m_swpref.set_identity();
  if(vnl_det(m_ref) > 0)
    {
    m_swpref(0,0) = -1.0;
    m_swpref(0,3) = (ref->GetBufferedRegion().GetSize(0) - 1) * ref->GetSpacing()[0];
    }

  m_swpmov.set_identity();
  if(vnl_det(m_mov) > 0)
    {
    m_swpmov(0,0) = -1.0;
    m_swpmov(0,3) = (mov->GetBufferedRegion().GetSize(0) - 1) * mov->GetSpacing()[0];
    }

  // Set the spacing matrices
  m_spcref.set_identity();
  m_spcmov.set_identity();
  for(size_t i = 0; i < 3; i++)
    {
    m_spcref(i,i) = ref->GetSpacing()[i];
    m_spcmov(i,i) = mov->GetSpacing()[i];
    }

  // Compute the output matrix
  m_out = 
    m_mov * vnl_inverse(m_spcmov) * m_swpmov *
    vnl_inverse(m_fsl) * 
    m_swpref * m_spcref * vnl_inverse(m_ref);

  // Put it on the stack
  vmat.pop_back();
  vmat.push_back(m_out);
}

void ras_inv(MatrixStack &vmat)
{
  MatrixType m = vmat.back();
  MatrixType minv = vnl_inverse(m);
  vmat.pop_back();
  vmat.push_back(minv);
}

void ras_mult(MatrixStack &vmat)
{
  MatrixType A = vmat[vmat.size() - 2];
  MatrixType B = vmat[vmat.size() - 1];
  MatrixType AB = A * B;
  vmat.pop_back();
  vmat.pop_back();
  vmat.push_back(AB);
}

int main(int argc, char *argv[])
{
  // Show usage
  if(argc < 2) return usage();

  // Set up the images that might be loaded
  ImageType::Pointer ref = NULL, src = NULL;

  // Set up the matrix stack
  vector<MatrixType> vmat;

  // Parse the command line
  for(int iarg = 1; iarg < argc; iarg++)
    {
    string arg = argv[iarg];
    if(arg == "-ref")
      {
      // Read the reference image
      typedef itk::ImageFileReader<ImageType> ReaderType;
      ReaderType::Pointer read_ref = ReaderType::New();
      read_ref->SetFileName(argv[++iarg]);
      read_ref->Update();
      ref = read_ref->GetOutput();
      }
    else if(arg == "-src" || arg == "-mov")
      {
      typedef itk::ImageFileReader<ImageType> ReaderType;
      // Read the moving image
      ReaderType::Pointer read_src = ReaderType::New();
      read_src->SetFileName(argv[++iarg]);
      read_src->Update();
      src = read_src->GetOutput();
      }
    else if(arg == "-fsl2ras")
      {
      // Convert FSL to RAS 
      fsl_to_ras(vmat, ref, src);
      }
    else if(arg == "-mult")
      {
      ras_mult(vmat);
      }
    else if(arg == "-inv")
      {
      ras_inv(vmat);
      }
    else if(arg == "-o")
      {
      ras_write(vmat, argv[++iarg]);
      }
    else if(arg[0] != '-')
      {
      ras_read(vmat, arg.c_str());
      }
    else
      {
      cerr << "Unknown option " << arg << endl;
      return usage();
      }
    }
  
}
