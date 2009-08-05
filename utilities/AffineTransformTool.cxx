#include "itkOrientedRASImage.h"
#include "itkImageFileReader.h"
#include <itkAffineTransform.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkTransformFactory.h>
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_det.h"
#include "vnl/vnl_inverse.h"
#include "vnl/algo/vnl_real_eigensystem.h"
#include "vnl/algo/vnl_qr.h"
#include <iostream>

using namespace std;

// Typedefs
typedef itk::OrientedRASImage<double, 3> ImageType;
typedef vnl_matrix_fixed<double, 4, 4> MatrixType;
typedef std::vector<MatrixType> MatrixStack;

int usage()
{
  cout << 
    "RAS Affine Transform Tool"
    "Usage: \n"
    "  c3d_affine_tool [transform_files | options] \n"
    "Options: \n"
    "  -fsl2ras      Convert FSL to RAS\n"
    "  -ref image    Set reference (fixed) image - only for -fsl2ras\n"
    "  -src image    Set source (moving) image - only for -fsl2ras\n"
    "  -sform image  Read matrix from NIfTI sform\n"
    "  -o matfile    Write output matrix\n"
    "  -det          Print the determinant\n"
    "  -inv          Invert matrix\n"
    "  -mult         Multiply matrices\n"
    "  -sqrt         Matrix square root (i.e., Q s.t. A = Q * Q)\n"
    "  -itk file     Import ITK transform\n"
    "  -oitk file    Export ITK transform\n"
    "  -info         Print matrix\n"
    ;
  return -1;
}

void itk_read(MatrixStack &vmat, const char *fname)
{
  typedef itk::MatrixOffsetTransformBase<double, 3, 3> MOTBType;
  typedef itk::AffineTransform<double, 3> AffTran;
  itk::TransformFactory<MOTBType>::RegisterTransform();
  itk::TransformFactory<AffTran>::RegisterTransform();

  itk::TransformFileReader::Pointer fltReader = itk::TransformFileReader::New();
  fltReader->SetFileName(fname);
  fltReader->Update();

  itk::TransformBase *base = fltReader->GetTransformList()->front();
  typedef itk::MatrixOffsetTransformBase<double, 3, 3> MOTBType;
  MOTBType *motb = dynamic_cast<MOTBType *>(base);

  MatrixType mat;
  mat.set_identity();
  if(motb)
    {
    for(size_t r = 0; r < 3; r++)
      {
      for(size_t c = 0; c < 3; c++)
        {
        mat(r,c) = motb->GetMatrix()(r,c);
        }
      mat(r,3) = motb->GetOffset()[r];
      }
    mat(2,0) *= -1; mat(2,1) *= -1; 
    mat(0,2) *= -1; mat(1,2) *= -1;
    mat(0,3) *= -1; mat(1,3) *= -1;
    vmat.push_back(mat);
    }
  else
    {
    cerr << "Unable to read ITK transform file" << endl;
    exit(-1);
    }
}

void itk_write(MatrixStack &vmat, const char *fname)
{
  // Get the current matrix
  MatrixType mat = vmat.back();
  
  // Flip the entries that must be flipped
  mat(2,0) *= -1; mat(2,1) *= -1; 
  mat(0,2) *= -1; mat(1,2) *= -1;
  mat(0,3) *= -1; mat(1,3) *= -1;

  // Create an ITK affine transform
  typedef itk::MatrixOffsetTransformBase<double, 3> AffTran;
  AffTran::Pointer atran = AffTran::New();

  // Populate its matrix
  AffTran::MatrixType amat = atran->GetMatrix();
  AffTran::OffsetType aoff = atran->GetOffset();

  for(size_t r = 0; r < 3; r++)
    {
    for(size_t c = 0; c < 3; c++)
      {
      amat(r,c) = mat(r,c);
      }
    aoff[r] = mat(r,3);
    }

  atran->SetMatrix(amat);
  atran->SetOffset(aoff);

  // Write the transform
  itk::TransformFileWriter::Pointer wrt = itk::TransformFileWriter::New();
  wrt->SetInput(atran);
  wrt->SetFileName(fname);
  wrt->Update();
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

void ras_sqrt(MatrixStack &vmat)
{
  MatrixType m = vmat.back();
  
  // Extract the rotation matrix
  vnl_matrix<double> A = m.extract(3,3);
   
  // Compute the eigensystem of the rotation matrix
  vnl_real_eigensystem eig(A);
  vnl_matrix<vcl_complex<double> > D = eig.D;
  vnl_matrix<vcl_complex<double> > V = eig.V;

  // Compute the square root of the rotation matrix
  D.apply(vcl_sqrt);
  vnl_matrix<vcl_complex<double> > QC = V * D * V.transpose();

  // Take the real part of QC
  vnl_matrix_fixed<double, 3, 3> Q, QI;
  for(size_t i = 0; i < 3; i++)
    for(size_t j = 0; j < 3; j++)
      Q(i,j) = vcl_real(QC(i,j));

  // Compute Q + I
  QI.set_identity();
  QI+=Q;

  // Compute the offset
  vnl_vector<double> b = m.get_column(3).extract(3);
  vnl_qr<double> qr(QI);
  vnl_vector<double> p = qr.solve(b);

  // Create the output matrix
  MatrixType Z;
  Z.set_identity();
  Z.update(Q, 0, 0);
  
  for(size_t i = 0; i < 3; i++)
    Z(i,3) = p(i);

  vmat.pop_back();
  vmat.push_back(Z);
}

void ras_det(MatrixStack &vmat)
{
  MatrixType m = vmat.back();
  double det = vnl_det(m);
  cout << "Det: " << det << endl;
}

void ras_sform(MatrixStack &vmat, const char *fname)
{
  typedef itk::OrientedRASImage<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fname);
  try 
    {
    reader->Update();
    }
  catch(std::exception &exc) 
    {
    cerr << "Can not read image " << fname << endl;
    cerr << "Reason: " << exc.what() << endl;
    }

  ImageType::Pointer img = reader->GetOutput();
 
  // Set up the directions
  MatrixType m = img->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix();
  vmat.push_back(m);
}

void ras_print(MatrixStack &vmat)
{
  MatrixType m = vmat.back();
  for(size_t i = 0; i < 4; i++)
    printf("%12.5f   %12.5f   %12.5f   %12.5f\n", m(i,0), m(i,1), m(i,2), m(i,3));
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
    else if(arg == "-det")
      {
      ras_det(vmat);
      }
    else if(arg == "-inv")
      {
      ras_inv(vmat);
      }
    else if(arg == "-sqrt")
      {
      ras_sqrt(vmat);
      }
    else if(arg == "-info")
      {
      ras_print(vmat);
      }
    else if(arg == "-sform")
      {
      ras_sform(vmat, argv[++iarg]);
      }
    else if(arg == "-itk")
      {
      itk_read(vmat, argv[++iarg]);
      }
    else if(arg == "-o")
      {
      ras_write(vmat, argv[++iarg]);
      }
    else if(arg == "-oitk")
      {
      itk_write(vmat, argv[++iarg]);
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
