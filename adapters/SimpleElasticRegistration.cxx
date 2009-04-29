#include "SimpleElasticRegistration.h"

template <class TPixel, unsigned int VDim>
void
SimpleElasticRegistration<TPixel, VDim>
::operator() ()
{
  // Get reference and moving images from stack
  this->imov = c->m_ImageStack.back(); c->m_ImageStack.pop_back();
  this->iref = c->m_ImageStack.back(); c->m_ImageStack.pop_back();
  
  // Create initial deformation field
  this->v = VectorImageType::New();
  this->v->SetRegions(iref->GetBufferedRegion());
  this->v->Allocate();

  // Initialize field to zero
  VectorType vzero; vzero.Fill(0.0f);
  this->v->FillBuffer(vzero);

  // Create the sampling function
  giref = GaussInterpolator::New();
  double sigma[] = {1.0, 1.0, 1.0};
  giref->SetParameters(sigma, 4.0);
  giref->SetInputImage(iref);

  gimov = GaussInterpolator::New();
  gimov->SetParameters(sigma, 4.0);
  gimov->SetInputImage(imov);
  
  // Compute the objective function
  *c->verbose << "Performing elastic registration" << endl;
  *c->verbose << "  Initial Objective " << ComputeObjective() << endl;

  // Do some processing ...
  // ImagePointer result = ...;
  
  // Put result on stack
  // c->m_ImageStack.pop_back();
  // c->m_ImageStack.push_back(result);
}

template <class TPixel, unsigned int VDim>
double
SimpleElasticRegistration<TPixel, VDim>
::ComputeObjective()
{
  // Compute the image objective
  double f_img = 0.0;
  Iterator itref(iref, iref->GetBufferedRegion());
  VectorIterator itv(v, v->GetBufferedRegion());
  for(; !itref.IsAtEnd(); ++itref, ++itv)
    {
    // Get the index in reference space
    itk::Point<double, VDim> p;
    itk::ContinuousIndex<double, VDim> q, qr;
    VectorType vi = itv.Get();
    IndexType idx_ref = itref.GetIndex();
    iref->TransformIndexToPhysicalPoint(idx_ref, p);
    iref->TransformPhysicalPointToContinuousIndex(p, qr);
    for(size_t d = 0; d < VDim; d++)
      p[d] += vi[d];
    imov->TransformPhysicalPointToContinuousIndex(p, q);
    
    // Get the two image values
    if(gimov->IsInsideBuffer(q))
      {
      double I_ref = giref->EvaluateAtContinuousIndex(qr);
      double I_mov = gimov->EvaluateAtContinuousIndex(q);
      f_img += (I_ref - I_mov) * (I_ref - I_mov);
      }
    } 

  return f_img / iref->GetBufferedRegion().GetNumberOfPixels();
}

template <class TPixel, unsigned int VDim>
double
SimpleElasticRegistration<TPixel, VDim>
::ComputeGradient()
{
  // Compute the gradient of the image objective
  return 0.0;



}


// Invocations
template class SimpleElasticRegistration<double, 2>;
template class SimpleElasticRegistration<double, 3>;
