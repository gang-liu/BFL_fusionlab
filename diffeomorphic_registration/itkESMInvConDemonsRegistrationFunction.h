/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkESMInvConDemonsRegistrationFunction.h,v $
  Language:  C++
  Date:      $Date: 2008-07-04 22:54:43 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkESMInvConDemonsRegistrationFunction_h
#define __itkESMInvConDemonsRegistrationFunction_h

#include "itkPDEDeformableRegistrationFunction.h"
#include "itkPoint.h"
#include "itkCovariantVector.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkWarpImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageDuplicator.h"

namespace itk {

/**
 * \class ESMInvConDemonsRegistrationFunction
 *
 * \brief Fast implementation of the symmetric demons registration force
 *
 * This class provides a substantially faster implementation of the inverse consistent
 * symmetric demons registration force. Speed is improved by keeping
 * a deformed copy of the moving image for gradient evaluation.
 *
 * Symmetric forces simply means using the mean of the gradient
 * of the fixed image and the gradient of the warped moving
 * image.
 *
 * Note that this class also enables the use of fixed, mapped moving
 * and warped moving images forces by using a call to SetUseGradientType
 *
 * The moving image should not be saturated. We indeed use
 * Numeritypedef typename itk::NearestNeighborInterpolateImageFunction<MovingImageType,CoordRepType> DefaultInterpolatorType;cTraits<MovingPixelType>::Max() as a special value.
 *
 * \author Tom Vercauteren, INRIA & Mauna Kea Technologies
 *
 * This implementation was taken from the Insight Journal paper:
 * http://hdl.handle.net/1926/510
 *
 * \sa SymmetricForcesDemonsRegistrationFunction
 * \sa SymmetricForcesDemonsRegistrationFilter
 * \sa DemonsRegistrationFilter
 * \sa DemonsRegistrationFunction
 * \ingroup FiniteDifferenceFunctions
 *
 */
template<class TFixedImage, class TMovingImage, class TDeformationField>
class ITK_EXPORT ESMInvConDemonsRegistrationFunction : 
    public PDEDeformableRegistrationFunction< TFixedImage,
                                              TMovingImage, TDeformationField>
{
public:
  /** Standard class typedefs. */
  typedef ESMInvConDemonsRegistrationFunction               Self;
  typedef PDEDeformableRegistrationFunction< 
    TFixedImage, TMovingImage, TDeformationField >    Superclass;
  typedef SmartPointer<Self>                          Pointer;
  typedef SmartPointer<const Self>                    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( ESMInvConDemonsRegistrationFunction, 
                PDEDeformableRegistrationFunction );

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType      MovingImageType;
  typedef typename Superclass::MovingImagePointer   MovingImagePointer;
  typedef typename MovingImageType::PixelType       MovingPixelType;

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType       FixedImageType;
  typedef typename Superclass::FixedImagePointer    FixedImagePointer;
  typedef typename FixedImageType::PixelType        FixedPixelType;
  typedef typename FixedImageType::IndexType        IndexType;
  typedef typename FixedImageType::SizeType         SizeType;
  typedef typename FixedImageType::SpacingType      SpacingType;
  
  /** Deformation field type. */
  typedef typename Superclass::DeformationFieldType    DeformationFieldType;
  typedef typename Superclass::DeformationFieldTypePointer   
  DeformationFieldTypePointer;

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** Inherit some enums from the superclass. */
  typedef typename Superclass::PixelType            PixelType;
  typedef typename Superclass::RadiusType           RadiusType;
  typedef typename Superclass::NeighborhoodType     NeighborhoodType;
  typedef typename Superclass::FloatOffsetType      FloatOffsetType;
  typedef typename Superclass::TimeStepType         TimeStepType;
  //FixedImagePointer         m_sdm_movingImage;
 // FixedImagePointer         m_sdm_fixedImage;
  /** Interpolator type. */
  typedef double                                    CoordRepType;
  typedef InterpolateImageFunction<
    MovingImageType,CoordRepType>                   InterpolatorType;
  typedef typename InterpolatorType::Pointer        InterpolatorPointer;
  typedef typename InterpolatorType::PointType      PointType;
  //typedef LinearInterpolateImageFunction<
   // MovingImageType,CoordRepType>                   DefaultInterpolatorType;
  typedef typename itk::NearestNeighborInterpolateImageFunction<MovingImageType,CoordRepType> DefaultInterpolatorType;

  /** Warper type */
  typedef WarpImageFilter<
    MovingImageType,
    MovingImageType,DeformationFieldType>           WarperType;
  typedef typename WarperType::Pointer              WarperPointer;

  /** Covariant vector type. */
  typedef CovariantVector<double,itkGetStaticConstMacro(ImageDimension)> CovariantVectorType;

  /** Fixed image gradient calculator type. */
  typedef CentralDifferenceImageFunction<FixedImageType> GradientCalculatorType;
  typedef typename GradientCalculatorType::Pointer   GradientCalculatorPointer;

  /** Moving image gradient (unwarped) calculator type. */
  typedef CentralDifferenceImageFunction<MovingImageType,CoordRepType>
      MovingImageGradientCalculatorType;
  typedef typename MovingImageGradientCalculatorType::Pointer
      MovingImageGradientCalculatorPointer;
  
    /** Set the deformation field image. */
  void SetInvDeformationField(  DeformationFieldTypePointer ptr )
    { m_invDeformationField = ptr; }

  /** Get the deformation field. This function should have been
   *  declared const. It is not for backward compatibility reasons. */
  DeformationFieldType * GetInvDeformationField(void)
    { return m_invDeformationField; }
  
  
    /* Set/Get JacobianDetField */
    
  void SetJacobianDetImage( const FixedImageType * ptr){
      m_JacobianDetImage = ptr;
  }
  
 const FixedImageType * GetJacobianDetImage(void) const{
      return m_JacobianDetImage;
  }
  
   /* Set/Get ForwardWeight  */
    
  void SetFwWeightImage( const FixedImageType * ptr){
      m_FwWeightImage = ptr;
  }
  
   const FixedImageType * GetFwWeightImage(void) const{
      return m_FwWeightImage;
  }
   /* Set/Get FixedSDMImage  */
  void SetFixedSDMImage( const FixedImageType * ptr){
      m_sdm_fixedImage = ptr;
  }
   const FixedImageType * GetFixedSDMImage(void) const{
      return m_sdm_fixedImage;
  }
   /* Set/Get MovingSDMImage  */
   void SetMovingSDMImage( const FixedImageType * ptr){
      m_sdm_movingImage = ptr;
  }
   const FixedImageType * GetMovingSDMImage(void) const{
      return m_sdm_movingImage;
  }
   /* Set/Get orignalFixedSDMImage  */
  void SetorignalFixedSDMImage( const FixedImageType * ptr){
      m_sdm_orignalfixedImage = ptr;
  }
   const FixedImageType * GetorignalFixedSDMImage(void) const{
      return m_sdm_orignalfixedImage;
  }
   /* Set/Get orignalMovingSDMImage  */
   void SetorignalMovingSDMImage( const FixedImageType * ptr){
      m_sdm_orignalmovingImage = ptr;
  }
   const FixedImageType * GetorignalMovingSDMImage(void) const{
      return m_sdm_orignalmovingImage;
  }
   
    /* Set/Get orignalFixedImage  */
  //void SetorignalFixedImage( const FixedImageType * ptr){
  //    m_orignalfixedImage = ptr;
 // }
  // const FixedImageType * GetorignalFixedImage(void) const{
  //    return m_orignalfixedImage;
  ///}
  /* Set/Get orignalMovingImage  */
 //  void SetorignalMovingImage( const FixedImageType * ptr){
  //    m_orignalmovingImage = ptr;
 // }
 //  const FixedImageType * GetorignalMovingImage(void) const{
 //     return m_orignalmovingImage;
  //}
   
  /** Set the moving image interpolator. */
  void SetMovingImageInterpolator( InterpolatorType * ptr )
  { m_MovingImageInterpolator = ptr; m_MovingImageWarper->SetInterpolator( ptr ); }
   
  /** Get the moving image interpolator. */
  InterpolatorType * GetMovingImageInterpolator(void)
  { return m_MovingImageInterpolator; }
  
  /** Set the moving image interpolator. */
  void SetFixedImageInterpolator( InterpolatorType * ptr )
  { m_FixedImageInterpolator = ptr; m_FixedImageWarper->SetInterpolator( ptr ); }
   
  /** Get the moving image interpolator. */
  InterpolatorType * GetFixedImageInterpolator(void)
  { return m_FixedImageInterpolator; }

  /** This class uses a constant timestep of 1. */
  virtual TimeStepType ComputeGlobalTimeStep(void * itkNotUsed(GlobalData)) const
  { return m_TimeStep; }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  virtual void *GetGlobalDataPointer() const
    {
    GlobalDataStruct *global = new GlobalDataStruct();
    global->m_SumOfSquaredDifference  = 0.0;
    global->m_NumberOfPixelsProcessed = 0L;
    global->m_SumOfSquaredChange      = 0;
    return global;
    }

  /** Release memory for global data structure. */
  virtual void ReleaseGlobalDataPointer( void *GlobalData ) const;
  void SignedDistanceMap_fixedImage();
  void SignedDistanceMap_movingImage();
  void SignedDistanceMap_orignalfixedImage();
  void SignedDistanceMap_orignalmovingImage();

  /** Set the object's state before each iteration. */
  virtual void InitializeIteration();

  /** This method is called by a finite difference solver image filter at
   * each pixel that does not lie on a data set boundary */
  virtual PixelType  ComputeUpdate(const NeighborhoodType &neighborhood,
                                   void *globalData,
                                   const FloatOffsetType &offset = FloatOffsetType(0.0));

  /** Get the metric value. The metric value is the mean square difference 
   * in intensity between the fixed image and transforming moving image 
   * computed over the the overlapping region between the two images. */
  virtual double GetMetric() const
    { return m_Metric; }

  /** Get the rms change in deformation field. */
  virtual const double &GetRMSChange() const
    { return m_RMSChange; }

  /** Set/Get the threshold below which the absolute difference of
   * intensity yields a match. When the intensities match between a
   * moving and fixed image pixel, the update vector (for that
   * iteration) will be the zero vector. Default is 0.001. */
  virtual void SetIntensityDifferenceThreshold(double);
  virtual double GetIntensityDifferenceThreshold() const;

  /** Set/Get the maximum update step length. In Thirion this is 0.5.
   *  Setting it to 0 implies no restriction (beware of numerical
   *  instability in this case. */
  virtual void SetMaximumUpdateStepLength(double sm)
    {
    this->m_MaximumUpdateStepLength =sm;
    }
  virtual double GetMaximumUpdateStepLength() const 
    {
    return this->m_MaximumUpdateStepLength;
    }

  /** Type of available image forces */
  enum GradientType {
     Symmetric=0,
     Fixed,
     WarpedMoving,
     MappedMoving
  };
  
  void SetUseJacobian(bool flag)
{
    m_UseJacobian = flag;
  }
  
  bool GetUseJacobian () const {
      
      return m_UseJacobian;
  }
  
   void SetUseFwWeight(bool flag)
  {
    m_UseFwWeight = flag;
  }
  
  bool GetUseFwWeight () const {
      
      return m_UseFwWeight;
  }
  
  // Set/Get RegWeight
  
  void SetRegWeight( double w){
      m_RegWeight = w;
  }
  
  double GetRegWeight(void){
      return m_RegWeight;
  }

  /** Set/Get the type of used image forces */
  virtual void SetUseGradientType( GradientType gtype )
    { m_UseGradientType = gtype; }
  virtual GradientType GetUseGradientType() const
    { return m_UseGradientType; }

protected:
  ESMInvConDemonsRegistrationFunction();
  ~ESMInvConDemonsRegistrationFunction() {}
  void PrintSelf(std::ostream& os, Indent indent) const;
  

  
  /** FixedImage image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<FixedImageType> FixedImageNeighborhoodIteratorType;

  /** A global data type for this class of equation. Used to store
   * iterators for the fixed image. */
  struct GlobalDataStruct
    {
    double          m_SumOfSquaredDifference;
    unsigned long   m_NumberOfPixelsProcessed;
    double          m_SumOfSquaredChange;
    };

private:
  ESMInvConDemonsRegistrationFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  /** Cache fixed image information. */
  SpacingType                     m_FixedImageSpacing;
  PointType                       m_FixedImageOrigin;
  double                          m_Normalizer;
  //void SignedDistanceMap_fixedImage(){}
  //void SignedDistanceMap_movingImage(){} 
  /** The deformation field. */
  DeformationFieldTypePointer     m_invDeformationField;
  
  /** Function to compute derivatives of the fixed image. */
  GradientCalculatorPointer       m_FixedImageGradientCalculator;

  /** Function to compute derivatives of the moving image (unwarped). */
  MovingImageGradientCalculatorPointer m_MappedMovingImageGradientCalculator;

  GradientType                    m_UseGradientType;

  /** Function to interpolate the moving image. */
  InterpolatorPointer             m_MovingImageInterpolator;

  /** Filter to warp moving image for fast gradient computation. */
  WarperPointer                   m_MovingImageWarper;
  
  /** Function to interpolate the moving image. */
  InterpolatorPointer             m_FixedImageInterpolator;

  /** Filter to warp moving image for fast gradient computation. */
  WarperPointer                   m_FixedImageWarper;

  /** The global timestep. */
  TimeStepType                    m_TimeStep;

  /** Threshold below which the denominator term is considered zero. */
  double                          m_DenominatorThreshold;

  /** Threshold below which two intensity value are assumed to match. */
  double                          m_IntensityDifferenceThreshold;

  /** Maximum update step length in pixels (default is 0.5 as in Thirion). */
  double                          m_MaximumUpdateStepLength;
  
  double                          m_RegWeight;
    /* Whether to use Jacobian weighting or not*/
  
  bool                          m_UseJacobian;
  bool                          m_UseFwWeight;
    
  FixedImagePointer         m_JacobianDetImage;
  FixedImagePointer         m_FwWeightImage;
  FixedImagePointer         m_sdm_movingImage;
  FixedImagePointer         m_sdm_fixedImage;
  FixedImagePointer         m_sdm_orignalmovingImage;
  FixedImagePointer         m_sdm_orignalfixedImage;
  //FixedImagePointer         m_orignalmovingImage;
  //FixedImagePointer         m_orignalfixedImage;
  /** The metric value is the mean square difference in intensity between
   * the fixed image and transforming moving image computed over the 
   * the overlapping region between the two images. */
  mutable double                  m_Metric;
  mutable double                  m_SumOfSquaredDifference;
  mutable unsigned long           m_NumberOfPixelsProcessed;
  mutable double                  m_RMSChange;
  mutable double                  m_SumOfSquaredChange;

  /** Mutex lock to protect modification to metric. */
  mutable SimpleFastMutexLock     m_MetricCalculationLock;

};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkESMInvConDemonsRegistrationFunction.txx"
#endif

#endif
