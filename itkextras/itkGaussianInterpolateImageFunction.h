/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaussianInterpolateImageFunction.h,v $
  Language:  C++
  Date:      $Date: 2009/04/29 16:55:34 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaussianInterpolateImageFunction_h
#define __itkGaussianInterpolateImageFunction_h

#include "itkInterpolateImageFunction.h"

namespace itk
{

/** \class GaussianInterpolateImageFunction
 * \brief Gaussianly interpolate an image at specified positions.
 *
 * GaussianInterpolateImageFunction linearly interpolates image intensity at
 * a non-integer pixel position. This class is templated
 * over the input image type and the coordinate representation type 
 * (e.g. float or double).
 *
 * This function works for N-dimensional images.
 *
 * \ingroup ImageFunctions ImageInterpolators 
 */
template <class TInputImage, class TCoordRep = double>
class ITK_EXPORT GaussianInterpolateImageFunction : 
  public InterpolateImageFunction<TInputImage,TCoordRep> 
{
public:
  /** Standard class typedefs. */
  typedef GaussianInterpolateImageFunction Self;
  typedef InterpolateImageFunction<TInputImage,TCoordRep> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(GaussianInterpolateImageFunction, InterpolateImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** OutputType typedef support. */
  typedef typename Superclass::OutputType OutputType;

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;

  /** RealType typedef support. */
  typedef typename Superclass::RealType RealType;

  /** Dimension underlying input image. */
  itkStaticConstMacro(VDim, unsigned int,Superclass::ImageDimension);

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Order specification (for derivatives) */
  typedef Index<VDim> OrderType;

  /** Set input */
  virtual void SetInputImage(const TInputImage *img)
    {
    // Call parent method
    Superclass::SetInputImage(img);

    if(img == NULL) return;

    // Set the bounding box
    for(size_t d = 0; d < VDim; d++)
      {
      bb_start[d] = 0;
      bb_end[d] = img->GetBufferedRegion().GetSize()[d] - 1;
      nt[d] = (int)(bb_end[d] - bb_start[d] + 0.5);
      dx[d] = new double[nt[d]];
      sf[d] = 1.0 / (sqrt(2.0) * sigma[d]);
      cut[d] = sigma[d] * alpha;
      }
    }

  void SetParameters(double *sigma, double alpha)
    {
    for(size_t d = 0; d < VDim; d++)
      this->sigma[d] = sigma[d];
    this->alpha = alpha;
    }


  void SetDerivativeOrder(OrderType order)
    {
    this->order = order;
    }

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the linearly interpolated image intensity at a 
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex( 
    const ContinuousIndexType & index ) const
    {
      // The bound variables for x, y, z
      int i0[VDim], i1[VDim];

      // Compute the ERF difference arrays
      for(size_t d = 0; d < VDim; d++)
        compute_erf_array(dx[d], i0[d], i1[d], bb_start[d], nt[d], cut[d], index[d], sf[d], order[d]);

      // Get a pointer to the output value
      double sum_me = 0.0, sum_m = 0.0;

      // Loop over the voxels in the region identified
      ImageRegion<VDim> region;
      for(size_t d = 0; d < VDim; d++)
        {
        region.SetIndex(d, i0[d]);
        region.SetSize(d, i1[d] - i0[d]);
        }
      
      for(
        ImageRegionConstIteratorWithIndex<InputImageType> it(this->GetInputImage(), region);
        !it.IsAtEnd(); ++it)
        {
        double w = 1.0;
        for(size_t d = 0; d < VDim; d++)
          {
          w *= dx[d][it.GetIndex()[d]];
          }
        sum_me += w * it.Get();
        sum_m += w;
        }

      return sum_me / sum_m;
    }

protected:
  GaussianInterpolateImageFunction()
    {
    for(size_t i = 0; i < VDim; i++)
      order[i] = 0.0;
    }
  ~GaussianInterpolateImageFunction(){};
  void PrintSelf(std::ostream& os, Indent indent) const
    { this->Superclass::PrintSelf(os,indent); }

private:
  GaussianInterpolateImageFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long  m_Neighbors;  


  double bb_start[VDim], bb_end[VDim], *dx[VDim], sf[VDim], cut[VDim];
  int nt[VDim], stride[VDim];
  double sigma[VDim], alpha;

  OrderType order;

  double erf_deriv(double t, int order) const
    {
    if(order == 0)
      return erf(t);
    if(order == 1)
      return 1.128379167095513 * exp(- t * t);
    if(order == 2)
      return -2.256758334191025 * exp(- t * t);
    return 0;
    }

  void compute_erf_array (
    double *dx_erf,         // The output array of erf(p+i+1) - erf(p+i)
    int &k0, int &k1,       // The range of integration 0 <= k0 < k1 <= n
    double b,               // Lower bound of the bounding box
    int n,                  // Size of the bounding box in steps
    double cut,             // The distance at which to cut off
    double p,               // the value p
    double sfac,            // scaling factor 1 / (Sqrt[2] sigma)
    int order = 0) const      
      {
      // Determine the range of voxels along the line where to evaluate erf
      k0 = (int) floor(p - b - cut);
      k1 = (int) ceil(p - b + cut);
      if(k0 < 0) k0 = 0;
      if(k1 > n) k1 = n;

      // Start at the first voxel
      double t = (b - p + k0) * sfac;
      double e_last = erf_deriv(t, order);
      for(int i = k0; i < k1; i++)
        {
        t += sfac;
        double e_now = erf_deriv(t, order);
        dx_erf[i] = e_now - e_last;
        e_last = e_now;
        }
      }

};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_GaussianInterpolateImageFunction(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT GaussianInterpolateImageFunction< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef GaussianInterpolateImageFunction< ITK_TEMPLATE_2 x > \
    GaussianInterpolateImageFunction##y; } \
}

#endif
