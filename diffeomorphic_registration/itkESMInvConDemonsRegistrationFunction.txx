/*=========================================================================
 
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkESMInvConDemonsRegistrationFunction.txx,v $
  Language:  C++
  Date:      $Date: 2008-07-11 19:02:04 $
  Version:   $Revision: 1.4 $
 
  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.
 
=========================================================================*/

#ifndef __itkESMInvConDemonsRegistrationFunction_txx
#define __itkESMInvConDemonsRegistrationFunction_txx

#include "itkESMInvConDemonsRegistrationFunction.h"
#include "itkExceptionObject.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/vnl_inverse.h"
#include "vnl/vnl_trace.h"
#include "vnl/algo/vnl_determinant.h"

#include "itkWarpImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"

namespace itk {
    
/**
 * Default constructor
 */
    template <class TFixedImage, class TMovingImage, class TDeformationField>
    ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::ESMInvConDemonsRegistrationFunction()
{
    
    RadiusType r;
    unsigned int j;
    for( j = 0; j < ImageDimension; j++ )
    {
        r[j] = 0;
    }
    this->SetRadius(r);
    
    m_TimeStep = 1.0;
    m_DenominatorThreshold = 1e-9;
    m_IntensityDifferenceThreshold = 0.001;
    m_MaximumUpdateStepLength = 0.5;
    m_RegWeight = 1000;
    
    this->SetMovingImage(NULL);
    this->SetFixedImage(NULL);
    this->SetInvDeformationField(NULL);
    this->SetJacobianDetImage(NULL);
    this->SetFwWeightImage(NULL);
    
    m_FixedImageSpacing.Fill( 1.0 );
    m_FixedImageOrigin.Fill( 0.0 );
    m_Normalizer = 0.0;
    m_FixedImageGradientCalculator = GradientCalculatorType::New();
    m_MappedMovingImageGradientCalculator = MovingImageGradientCalculatorType::New();
    
    this->m_UseGradientType = Symmetric;
    
    typename DefaultInterpolatorType::Pointer interp =
    DefaultInterpolatorType::New();
    
    m_MovingImageInterpolator = static_cast<InterpolatorType*>(
    interp.GetPointer() );
    
    m_MovingImageWarper = WarperType::New();
    m_MovingImageWarper->SetInterpolator( m_MovingImageInterpolator );
    m_MovingImageWarper->SetEdgePaddingValue( NumericTraits<MovingPixelType>::max() );
    
    typename DefaultInterpolatorType::Pointer interp_fixed =
    DefaultInterpolatorType::New();
    
    m_FixedImageInterpolator = static_cast<InterpolatorType*>(
    interp_fixed.GetPointer() );
    
    m_FixedImageWarper = WarperType::New();
    m_FixedImageWarper->SetInterpolator( m_FixedImageInterpolator );
    m_FixedImageWarper->SetEdgePaddingValue( NumericTraits<FixedPixelType>::max() );
    
    m_Metric = NumericTraits<double>::max();
    m_SumOfSquaredDifference = 0.0;
    m_NumberOfPixelsProcessed = 0L;
    m_RMSChange = NumericTraits<double>::max();
    m_SumOfSquaredChange = 0.0;
    m_UseJacobian = false;
    m_UseFwWeight = false;
    }
    
    
/*
 * Standard "PrintSelf" method.
 */
    template <class TFixedImage, class TMovingImage, class TDeformationField>
    void
    ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);
    
    os << indent << "UseGradientType: ";
    os << m_UseGradientType << std::endl;
    os << indent << "MaximumUpdateStepLength: ";
    os << m_MaximumUpdateStepLength << std::endl;
    
    os << indent << "MovingImageIterpolator: ";
    os << m_MovingImageInterpolator.GetPointer() << std::endl;
    os << indent << "FixedImageGradientCalculator: ";
    os << m_FixedImageGradientCalculator.GetPointer() << std::endl;
    os << indent << "MappedMovingImageGradientCalculator: ";
    os << m_MappedMovingImageGradientCalculator.GetPointer() << std::endl;
    os << indent << "DenominatorThreshold: ";
    os << m_DenominatorThreshold << std::endl;
    os << indent << "IntensityDifferenceThreshold: ";
    os << m_IntensityDifferenceThreshold << std::endl;
    
    os << indent << "Metric: ";
    os << m_Metric << std::endl;
    os << indent << "SumOfSquaredDifference: ";
    os << m_SumOfSquaredDifference << std::endl;
    os << indent << "NumberOfPixelsProcessed: ";
    os << m_NumberOfPixelsProcessed << std::endl;
    os << indent << "RMSChange: ";
    os << m_RMSChange << std::endl;
    os << indent << "SumOfSquaredChange: ";
    os << m_SumOfSquaredChange << std::endl;
    
    }
    
/**
 *
 */
    template <class TFixedImage, class TMovingImage, class TDeformationField>
    void
    ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::SetIntensityDifferenceThreshold(double threshold)
{
    m_IntensityDifferenceThreshold = threshold;
    }
    
/**
 *
 */
    template <class TFixedImage, class TMovingImage, class TDeformationField>
    double
    ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::GetIntensityDifferenceThreshold() const
{
    return m_IntensityDifferenceThreshold;
    }
    
/**
 * Set the function state values before each iteration
 */
    template <class TFixedImage, class TMovingImage, class TDeformationField>
    void
    ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::InitializeIteration()
{
    if( !this->GetMovingImage() || !this->GetFixedImage()
    || !m_MovingImageInterpolator )
    {
        itkExceptionMacro(
        << "MovingImage, FixedImage and/or Interpolator not set" );
    }
    
    // cache fixed image information
    m_FixedImageSpacing = this->GetFixedImage()->GetSpacing();
    m_FixedImageOrigin  = this->GetFixedImage()->GetOrigin();
    
    // compute the normalizer
    if( m_MaximumUpdateStepLength > 0.0 )
    {
        m_Normalizer = 0.0;
        for( unsigned int k = 0; k < ImageDimension; k++ )
        {
            m_Normalizer += m_FixedImageSpacing[k] * m_FixedImageSpacing[k];
        }
        m_Normalizer *= m_MaximumUpdateStepLength * m_MaximumUpdateStepLength /
        static_cast<double>( ImageDimension );
    }
    else
    {
        // set it to minus one to denote a special case
        // ( unrestricted update length )
        m_Normalizer = -1.0;
    }
    
    /*if ( this->m_UseJacobian && !this->GetJacobianDetImage()){
      itkExceptionMacro(
        << "Jacobian Image needs to be set ! ");  
    }*/
    // setup gradient calculator
    m_FixedImageGradientCalculator->SetInputImage( this->GetFixedImage() );
    m_MappedMovingImageGradientCalculator->SetInputImage( this->GetMovingImage() );
    

//this->SetorignalMovingImage( this->GetMovingImage() );
//this->SetorignalFixedImage( this->GetFixedImage() );


//typedef float      PixelType;
//typedef itk::Image< PixelType, ImageDimension >    ImageType;
//typedef itk::ImageFileWriter< ImageType >  WriterType;
//typename WriterType::Pointer writer3 = WriterType::New();
//writer3->SetFileName( "fixed_beforewarp.nii.gz" );
//writer3->SetInput( this->GetFixedImage() );
//writer3->Update();


    // Compute warped moving image
    m_MovingImageWarper->SetOutputSpacing( this->GetFixedImage()->GetSpacing() );
    m_MovingImageWarper->SetOutputOrigin( this->GetFixedImage()->GetOrigin() );
    m_MovingImageWarper->SetOutputDirection( this->GetFixedImage()->GetDirection() );
    m_MovingImageWarper->SetInput( this->GetMovingImage() );
    m_MovingImageWarper->SetDeformationField( this->GetDeformationField() );
    m_MovingImageWarper->GetOutput()->SetRequestedRegion( this->GetDeformationField()->GetRequestedRegion() );
    m_MovingImageWarper->Update();
 

  // typename WriterType::Pointer writer4 = WriterType::New();
//writer4->SetFileName( "fixed_afterwarp1.nii.gz" );
//writer4->SetInput( this->GetFixedImage() );
//writer4->Update();

    
    // Compute warped moving image
    m_FixedImageWarper->SetOutputSpacing( this->GetFixedImage()->GetSpacing() );
    m_FixedImageWarper->SetOutputOrigin( this->GetFixedImage()->GetOrigin() );
    m_FixedImageWarper->SetOutputDirection( this->GetFixedImage()->GetDirection() );
    m_FixedImageWarper->SetInput( this->GetFixedImage() );
    m_FixedImageWarper->SetDeformationField( this->GetInvDeformationField() );
    m_FixedImageWarper->GetOutput()->SetRequestedRegion( this->GetInvDeformationField()->GetRequestedRegion() );
    m_FixedImageWarper->Update();
    
  //typename WriterType::Pointer writer5 = WriterType::New();
//writer5->SetFileName( "fixed_afterwarp2.nii.gz" );
//writer5->SetInput( this->GetFixedImage() );
//writer5->Update();
 // typename WriterType::Pointer writer6 = WriterType::New();
//writer6->SetFileName( "fixed_warped.nii.gz" );
//writer6->SetInput( m_FixedImageWarper->GetOutput() );
//writer6->Update();




//std::cout << "54321" << std::endl;
SignedDistanceMap_fixedImage();
SignedDistanceMap_movingImage();
SignedDistanceMap_orignalfixedImage();
SignedDistanceMap_orignalmovingImage();

//std::cout << "2222222222222222222" << std::endl;
//  typedef float      PixelType;
 // typedef itk::Image< PixelType, ImageDimension >    ImageType;
  //typedef itk::ImageDuplicator< ImageType > DuplicatorType;
  //typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  //duplicator->SetInputImage(m_sdm_fixedImage.GetPointer());
  //duplicator->Update();



//std::cout << "12345" << std::endl;
    
//IndexType First = duplicator->GetOutput()->GetLargestPossibleRegion().GetIndex();
//float aa=duplicator->GetOutput()->GetPixel( First );
//std::cout << "aa=" << aa << std::endl;





    // setup moving image interpolator for further access
    m_MovingImageInterpolator->SetInputImage( this->GetMovingImage() );
    m_FixedImageInterpolator->SetInputImage( this->GetFixedImage() );
    
    // initialize metric computation variables
    m_SumOfSquaredDifference  = 0.0;
    m_NumberOfPixelsProcessed = 0L;
    m_SumOfSquaredChange      = 0.0;
    }
    
    
/**
 * Compute update at a non boundary neighbourhood
 */
    template <class TFixedImage, class TMovingImage, class TDeformationField>
    typename ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::PixelType
    ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::ComputeUpdate(const NeighborhoodType &it, void * gd,
    const FloatOffsetType& itkNotUsed(offset))
{
 
    GlobalDataStruct *globalData = (GlobalDataStruct *)gd;
    PixelType update;
    IndexType FirstIndex = this->GetFixedImage()->GetLargestPossibleRegion().GetIndex();
    IndexType LastIndex = this->GetFixedImage()->GetLargestPossibleRegion().GetIndex() +
    this->GetFixedImage()->GetLargestPossibleRegion().GetSize();
    
    const IndexType index = it.GetIndex();
    //std::cout << "222" << std::endl;
    // Get fixed image related information
    // Note: no need to check if the index is within
    // fixed image buffer. This is done by the external filter.
    const double fixedValue = static_cast<double>(this->GetorignalFixedSDMImage()->GetPixel( index ) );
    double jacobian = 1.0f;
    double fw_weight = 1.0f;
    if (m_UseJacobian)
    {
        jacobian = static_cast<double>(this->GetJacobianDetImage()->GetPixel(index));
    }
    
    if (m_UseFwWeight)
    {
        fw_weight = static_cast<double>(this->GetFwWeightImage()->GetPixel(index));
    }
 //std::cout << "333" << std::endl;
    //MovingPixelType warpedFixedPixValue = m_FixedImageWarper->GetOutput()->GetPixel( index );
    //MovingPixelType warpedFixedPixValue= m_sdm_fixedImage->GetPixel( index );
 MovingPixelType warpedFixedPixValue= this->GetFixedSDMImage()->GetPixel( index );
        //std::cout << "444" << std::endl;
    if( warpedFixedPixValue == NumericTraits <FixedPixelType>::max() )
    {
        update.Fill( 0.0 );
        return update;
    }
    double warpedFixedValue = static_cast<double>( warpedFixedPixValue );
    //std::cout << "3" << std::endl;
    const double movingValue = static_cast<double>(this->GetorignalMovingSDMImage()->GetPixel( index ) );
    // Get moving image related information
    // check if the point was mapped outside of the moving image using
    // the "special value" NumericTraits<MovingPixelType>::max()
    
    //MovingPixelType warpedMovingPixValue = m_MovingImageWarper->GetOutput()->GetPixel( index );
    //MovingPixelType warpedMovingPixValue = m_sdm_movingImage->GetPixel( index );
MovingPixelType warpedMovingPixValue = this->GetMovingSDMImage()->GetPixel( index );
    if( warpedMovingPixValue == NumericTraits <MovingPixelType>::max() )
    {
        update.Fill( 0.0 );
        return update;
    }
    
    
    double warpedMovingValue = static_cast<double>( warpedMovingPixValue );
    //const double movingValue = static_cast<double>( movingPixValue );
    
    CovariantVectorType usedGradientTimes2;
    
    // we don't use a CentralDifferenceImageFunction here to be able to
    // check for NumericTraits<MovingPixelType>::max()
    CovariantVectorType warpedMovingGradient;
    IndexType tmpIndex = index;
    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
    {
        // bounds checking
        if( FirstIndex[dim]==LastIndex[dim]  || index[dim]<FirstIndex[dim] || index[dim]>=LastIndex[dim] )
        {
            warpedMovingGradient[dim] = 0.0;
            continue;
        }
        else if ( index[dim] == FirstIndex[dim] )
        {
            // compute derivative
            tmpIndex[dim] += 1;
            //warpedMovingPixValue = m_MovingImageWarper->GetOutput()->GetPixel( tmpIndex );
            //warpedMovingPixValue = m_sdm_movingImage->GetPixel( tmpIndex );
warpedMovingPixValue = this->GetMovingSDMImage()->GetPixel( tmpIndex );
            if( warpedMovingPixValue == NumericTraits <MovingPixelType>::max() )
            {
                // weird crunched border case
                warpedMovingGradient[dim] = 0.0;
            }
            else
            {
                // forward difference
                warpedMovingGradient[dim] = static_cast<double>( warpedMovingPixValue ) - warpedMovingValue;
                warpedMovingGradient[dim] /= m_FixedImageSpacing[dim];
            }
            tmpIndex[dim] -= 1;
            continue;
        }
        else if ( index[dim] == (LastIndex[dim]-1) )
        {
            // compute derivative
            tmpIndex[dim] -= 1;
            //warpedMovingPixValue = m_MovingImageWarper->GetOutput()->GetPixel( tmpIndex );
            //warpedMovingPixValue = m_sdm_movingImage->GetPixel( tmpIndex );
warpedMovingPixValue = this->GetMovingSDMImage()->GetPixel( tmpIndex );
            if( warpedMovingPixValue == NumericTraits<MovingPixelType>::max() )
            {
                // weird crunched border case
                warpedMovingGradient[dim] = 0.0;
            }
            else
            {
                // backward difference
                warpedMovingGradient[dim] = warpedMovingValue - static_cast<double>( warpedMovingPixValue );
                warpedMovingGradient[dim] /= m_FixedImageSpacing[dim];
            }
            tmpIndex[dim] += 1;
            continue;
        }
        
       // std::cout << "4" << std::endl;
        // compute derivative
        tmpIndex[dim] += 1;
        //warpedMovingPixValue = m_MovingImageWarper->GetOutput()->GetPixel( tmpIndex );
        //warpedMovingPixValue = m_sdm_movingImage->GetPixel( tmpIndex );
warpedMovingPixValue = this->GetMovingSDMImage()->GetPixel( tmpIndex );
        if ( warpedMovingPixValue == NumericTraits<MovingPixelType>::max() )
        {
            // backward differencem_MovingImageWarper->GetOutput()->GetPixel( 
            warpedMovingGradient[dim] = warpedMovingValue;
    
            tmpIndex[dim] -= 2;
            //warpedMovingPixValue = m_MovingImageWarper->GetOutput()->GetPixel( tmpIndex );
            //warpedMovingPixValue = m_sdm_movingImage->GetPixel( tmpIndex );
warpedMovingPixValue = this->GetMovingSDMImage()->GetPixel( tmpIndex );
            if( warpedMovingPixValue == NumericTraits<MovingPixelType>::max() )
            {
            // weird crunched border case
                warpedMovingGradient[dim] = 0.0;
            }
            else
            {
                // backward difference
                warpedMovingGradient[dim] -= static_cast<double>(
                this->GetMovingSDMImage()->GetPixel( tmpIndex ) );
        
                warpedMovingGradient[dim] /= m_FixedImageSpacing[dim];
            }
        }
        else
        {
            warpedMovingGradient[dim] = static_cast<double>( warpedMovingPixValue );
            
            tmpIndex[dim] -= 2;
            //warpedMovingPixValue = m_MovingImageWarper->GetOutput()->GetPixel( tmpIndex );
            //warpedMovingPixValue = m_sdm_movingImage->GetPixel( tmpIndex );
warpedMovingPixValue = this->GetMovingSDMImage()->GetPixel( tmpIndex );
            if ( warpedMovingPixValue == NumericTraits<MovingPixelType>::max() )
            {
                // forward difference
                warpedMovingGradient[dim] -= warpedMovingValue;
                warpedMovingGradient[dim] /= m_FixedImageSpacing[dim];
            }
            else
            {
                // normal case, central difference
                warpedMovingGradient[dim] -= static_cast<double>( warpedMovingPixValue );
                warpedMovingGradient[dim] *= 0.5 / m_FixedImageSpacing[dim];
            }
        }
        tmpIndex[dim] += 1;
    }
    //std::cout << "5" << std::endl;
    CovariantVectorType warpedFixedGradient;
    IndexType tmpIndex2 = index;
    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
    {
        // bounds checking
        if( FirstIndex[dim]==LastIndex[dim] || index[dim]<FirstIndex[dim] || index[dim]>=LastIndex[dim] )
        {
            warpedFixedGradient[dim] = 0.0;
            continue;
        }
        else if ( index[dim] == FirstIndex[dim] )
        {
            // compute derivative
            tmpIndex2[dim] += 1;
            //warpedFixedPixValue = m_FixedImageWarper->GetOutput()->GetPixel( tmpIndex2 );
            //warpedFixedPixValue = m_sdm_fixedImage->GetPixel( tmpIndex2 );
warpedFixedPixValue = this->GetFixedSDMImage()->GetPixel( tmpIndex2 );
            if( warpedFixedPixValue == NumericTraits <FixedPixelType>::max() )
            {
                // weird crunched border case
                warpedFixedGradient[dim] = 0.0;
            }
            else
            {
                // forward difference
                warpedFixedGradient[dim] = static_cast<double>( warpedFixedPixValue ) - warpedFixedValue;
                warpedFixedGradient[dim] /= m_FixedImageSpacing[dim];
            }
            tmpIndex2[dim] -= 1;
            continue;
        }
        else if ( index[dim] == (LastIndex[dim]-1) )
        {
            // compute derivative
            tmpIndex2[dim] -= 1;
            //warpedFixedPixValue = m_FixedImageWarper->GetOutput()->GetPixel( tmpIndex2 );
            //warpedFixedPixValue = m_sdm_fixedImage->GetPixel( tmpIndex2 );
warpedFixedPixValue = this->GetFixedSDMImage()->GetPixel( tmpIndex2 );
            if( warpedFixedPixValue == NumericTraits<FixedPixelType>::max() )
            {
                // weird crunched border case
                warpedFixedGradient[dim] = 0.0;
            }
            else
            {
                // backward difference
                warpedFixedGradient[dim] = warpedFixedValue - static_cast<double>( warpedFixedPixValue );
                warpedFixedGradient[dim] /= m_FixedImageSpacing[dim];
            }
            tmpIndex2[dim] += 1;
            continue;
        }
        
        
        // compute derivative
        tmpIndex2[dim] += 1;
        //warpedFixedPixValue = m_FixedImageWarper->GetOutput()->GetPixel( tmpIndex2 );
       // warpedFixedPixValue = m_sdm_fixedImage->GetPixel( tmpIndex2 );
warpedFixedPixValue = this->GetFixedSDMImage()->GetPixel( tmpIndex2 );
        if ( warpedFixedPixValue == NumericTraits<FixedPixelType>::max() )
        {
            // backward difference
            warpedFixedGradient[dim] = warpedFixedValue;
    
            tmpIndex2[dim] -= 2;
            //warpedFixedPixValue = m_FixedImageWarper->GetOutput()->GetPixel( tmpIndex2 );
            //warpedFixedPixValue = m_sdm_fixedImage->GetPixel( tmpIndex2 );
warpedFixedPixValue = this->GetFixedSDMImage()->GetPixel( tmpIndex2 );
            if( warpedFixedPixValue == NumericTraits<FixedPixelType>::max() )
            {
                // weird crunched border case
                 warpedFixedGradient[dim] = 0.0;
            }
            else
            {
                // backward difference
                warpedFixedGradient[dim] -= static_cast<double>(
                this->GetFixedSDMImage()->GetPixel( tmpIndex2 ) );
        
                warpedFixedGradient[dim] /= m_FixedImageSpacing[dim];
            }
        }
        else
        {
            warpedFixedGradient[dim] = static_cast<double>( warpedFixedPixValue );
            
            tmpIndex2[dim] -= 2;
            //warpedFixedPixValue = m_FixedImageWarper->GetOutput()->GetPixel( tmpIndex2 );
           // warpedFixedPixValue = m_sdm_fixedImage->GetPixel( tmpIndex2 );
warpedFixedPixValue = this->GetFixedSDMImage()->GetPixel( tmpIndex2 );
            if ( warpedFixedPixValue == NumericTraits<FixedPixelType>::max() )
            {
                // forward difference
                warpedFixedGradient[dim] -= warpedFixedValue;
                warpedFixedGradient[dim] /= m_FixedImageSpacing[dim];
            }
            else
            {
                // normal case, central difference
                warpedFixedGradient[dim] -= static_cast<double>( warpedFixedPixValue );
                warpedFixedGradient[dim] *= 0.5 / m_FixedImageSpacing[dim];
            }
        }
        tmpIndex2[dim] += 1;
    }
    
    vnl_vector<double> warpedFixedGradient_vec(ImageDimension);
    vnl_vector<double> warpedMovingGradient_vec(ImageDimension);
    
    vnl_matrix<double> Hess(ImageDimension, ImageDimension);
    vnl_matrix<double> invHess(ImageDimension, ImageDimension);
    
    
    const double speedValue = fixedValue - warpedMovingValue;
    const double speedValue2 = warpedFixedValue - movingValue;
    
    //std::cout << "Reg Weight (inside) " << this->m_RegWeight << std::endl;
    
    vnl_diag_matrix<double> Reg(ImageDimension, this->m_RegWeight);
    warpedFixedGradient_vec = warpedFixedGradient.GetVnlVector();
    warpedMovingGradient_vec = warpedMovingGradient.GetVnlVector();
    
     
    //std::cout << "Jacobian " << jacobian << "Fw weight: " << fw_weight << std::endl;
    
    Hess = 0.25*fw_weight*outer_product(warpedMovingGradient_vec, warpedMovingGradient_vec) + 
        jacobian * 0.25*outer_product(warpedFixedGradient_vec, warpedFixedGradient_vec);
    
    //Hess = 0.25*outer_product(warpedMovingGradient_vec, warpedMovingGradient_vec) + 
    //    0.25*outer_product(warpedFixedGradient_vec, warpedFixedGradient_vec);
    
    if (jacobian == 0)
    {
       Hess = 4.0*Hess;
    }
    
    invHess = vnl_inverse(Hess + Reg);
    vnl_matrix<double> temp_mat(ImageDimension,1);



          temp_mat.set_column(0,0.25*fw_weight*(fixedValue - warpedMovingValue)*warpedMovingGradient_vec + 
                jacobian * 0.25*(warpedFixedValue - movingValue)*warpedFixedGradient_vec);
    

    if (jacobian == 0)
    {
        temp_mat = 4.0*temp_mat;
    }
    
    vnl_matrix<double> update_mat(ImageDimension,1);
    
    update_mat = invHess * temp_mat;
    
           
  /**
   * Compute Update.
   * We avoid the mismatch in units between the two terms.
   * and avoid large step using a normalization term.
   */
        
    if ( vnl_math_abs(speedValue) < m_IntensityDifferenceThreshold || vnl_math_abs(speedValue2) < m_IntensityDifferenceThreshold)
    {
        update.Fill( 0.0 );
    }
    else
    {
        for( unsigned int dim = 0; dim < ImageDimension; dim++ )
        {
            update[dim] = update_mat[dim][0];
            
        }
    }
    
        
    // WARNING!! We compute the global data without taking into account the current update step.
    // There are several reasons for that: If an exponential, a smoothing or any other operation
    // is applied on the update field, we cannot compute the newMappedCenterPoint here; and even
    // if we could, this would be an often unnecessary time-consuming task.
    if ( globalData )
    {
        globalData->m_SumOfSquaredDifference += vnl_math_sqr( speedValue );
        globalData->m_NumberOfPixelsProcessed += 1;
        globalData->m_SumOfSquaredChange += update.GetSquaredNorm();
    }
    
    return update;
}
 

  template <class TFixedImage, class TMovingImage, class TDeformationField>
void ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::SignedDistanceMap_movingImage()
{

float eps = 0.0001;
typedef float      PixelType;
typedef unsigned int maximum,minimum;
//const   unsigned int        Dimension = 3;
typedef itk::Image< PixelType, ImageDimension >    ImageType;
//typedef itk::Image<float, ImageDimension>          FloatImageType;

typedef  itk::SignedMaurerDistanceMapImageFilter< ImageType, ImageType  > SignedMaurerDistanceMapImageFilterType;

  typedef itk::ImageDuplicator< ImageType > DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(m_MovingImageWarper->GetOutput());
  duplicator->Update();
  typename ImageType::Pointer image_label = duplicator->GetOutput();

  typename ImageType::Pointer m_movingImage = ImageType::New();
  typename ImageType::Pointer image2 = duplicator->GetOutput();

   typedef itk::CastImageFilter< ImageType, ImageType > CastFilterType;
   typename CastFilterType::Pointer castFilter = CastFilterType::New();
   castFilter->SetInput(image2);
   castFilter->Update();
   m_movingImage = castFilter->GetOutput();

itk::ImageRegionIterator<ImageType> imageIterator_1(m_movingImage,m_movingImage->GetLargestPossibleRegion());
while(!imageIterator_1.IsAtEnd())
  {
     imageIterator_1.Set(0);
       ++imageIterator_1;
  }


//get the array of labels
unsigned int array_label[100]={0};
unsigned int i=0,j,pixelValue_label;
itk::ImageRegionIterator<ImageType> imageIterator_2(m_MovingImageWarper->GetOutput(),m_MovingImageWarper->GetOutput()->GetLargestPossibleRegion());
while(!imageIterator_2.IsAtEnd())
  {
        pixelValue_label = imageIterator_2.Get();
	  unsigned int k=0;
	  if (pixelValue_label !=0)
  	    {
    		for(j=1;j<100;j++)
     		   {
       	              if( pixelValue_label == array_label[j-1])
         	         {
               		    k=k+1;
         		 }
     		    }
   		if(k==0 )
        	 {
                    array_label[i]=pixelValue_label;
                   // std::cout << array_label[i]<< std::endl;
                    i=i+1;
         	  }

            }
   ++imageIterator_2;
  }


// calculate sdt
unsigned int m;
for(m=0;array_label[m]>0;m++)
{
//get the single label region
itk::ImageRegionIterator<ImageType> imageIterator_3(m_MovingImageWarper->GetOutput(),m_MovingImageWarper->GetOutput()->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_4(image_label,image_label->GetLargestPossibleRegion());
while(!imageIterator_3.IsAtEnd())
  {    
      unsigned int pixelValue = imageIterator_3.Get();
      if (pixelValue == array_label[m]) 
        {
        imageIterator_4.Set(0);
        }
      else
        {
	  imageIterator_4.Set(1);
	}
    ++imageIterator_3;
    ++imageIterator_4;
  }
 // sdt the single label region
typename SignedMaurerDistanceMapImageFilterType::Pointer distanceMapImageFilter = SignedMaurerDistanceMapImageFilterType::New();
distanceMapImageFilter->SetInput(image_label);

//get the minimum and maximum in the single label sdm

typename ImageType::Pointer image_sdt = ImageType::New();
image_sdt = distanceMapImageFilter->GetOutput();// because the array_label have a lot zeros label
image_sdt->Update();

float minimum=0, maximum, pixelValue_am;
bool flag1=true;
itk::ImageRegionIterator<ImageType> imageIterator_5(image_sdt,image_sdt->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_6(image_label,image_label->GetLargestPossibleRegion());
while(!imageIterator_5.IsAtEnd())
  {
        pixelValue_am = imageIterator_5.Get();
           if (imageIterator_6.Get()==0)
           {
             if(flag1)
               {
                      flag1=false;
              //        minimum = pixelValue_am;
                      maximum = pixelValue_am;
                }
             if(pixelValue_am < minimum)
                {
              //    minimum = pixelValue_am;
                }
             else if(pixelValue_am > maximum)
                {
                  maximum = pixelValue_am;
                }
           }
    ++imageIterator_5;
    ++imageIterator_6;
    }

// norm the sdt single region and send it to totalsignedmap
float pixelValue_sdt, pixelValue_sdt_norm,diff;
diff=(maximum-minimum)+eps;  // may be zeros, should be think again;
//std::cout<< "moving the diff is: "<<diff<<std::endl;
//std::cout<< "the minmum is: "<<minimum<<std::endl;
bool flag2; // for debug only
itk::ImageRegionIterator<ImageType> imageIterator_7(image_sdt,image_sdt->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_8(image_label,image_label->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_9(m_movingImage,m_movingImage->GetLargestPossibleRegion());
while(!imageIterator_7.IsAtEnd())
  {
        pixelValue_sdt = imageIterator_7.Get();
        if (imageIterator_8.Get()==0)
	{
	      pixelValue_sdt_norm = (pixelValue_sdt-minimum)/diff;
          pixelValue_sdt_norm = pixelValue_sdt_norm*255;
	      imageIterator_9.Set(pixelValue_sdt_norm);
	}
    ++imageIterator_7;
    ++imageIterator_8;
    ++imageIterator_9;
  }

}

this->SetMovingSDMImage(m_movingImage);
} 
    

 template <class TFixedImage, class TMovingImage, class TDeformationField>

  void  ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::SignedDistanceMap_fixedImage()
{

float eps = 0.0001;
typedef float      PixelType;
typedef unsigned int maximum,minimum;
//const   unsigned int        Dimension = 3;
typedef itk::Image< PixelType, ImageDimension >    ImageType;
//typedef itk::Image<float, ImageDimension>          FloatImageType;

typedef  itk::SignedMaurerDistanceMapImageFilter< ImageType, ImageType  > SignedMaurerDistanceMapImageFilterType;

  typedef itk::ImageDuplicator< ImageType > DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(m_FixedImageWarper->GetOutput());
  duplicator->Update();
  typename ImageType::Pointer image_label = duplicator->GetOutput();


  typename ImageType::Pointer m_fixedImage = ImageType::New();
  typename ImageType::Pointer image2 = duplicator->GetOutput();
  
//typedef itk::ImageFileWriter< ImageType >  WriterType;
//typename WriterType::Pointer writer = WriterType::New();
//writer->SetFileName( "test_warpedfixedimage.nii.gz" );
//writer->SetInput( m_FixedImageWarper->GetOutput() );
//writer->Update();



//image_label = image;

   typedef itk::CastImageFilter< ImageType, ImageType > CastFilterType;
   typename CastFilterType::Pointer castFilter = CastFilterType::New();
   castFilter->SetInput(image2);
   castFilter->Update();
   m_fixedImage = castFilter->GetOutput();

itk::ImageRegionIterator<ImageType> imageIterator_1(m_fixedImage,m_fixedImage->GetLargestPossibleRegion());
while(!imageIterator_1.IsAtEnd())
  {
     imageIterator_1.Set(0);
       ++imageIterator_1;
  }


//get the array of labels
unsigned int array_label[100]={0};
unsigned int i=0,j,pixelValue_label;
itk::ImageRegionIterator<ImageType> imageIterator_2(m_FixedImageWarper->GetOutput(),m_FixedImageWarper->GetOutput()->GetLargestPossibleRegion());
while(!imageIterator_2.IsAtEnd())
  {
        pixelValue_label = imageIterator_2.Get();
	  unsigned int k=0;
	  if (pixelValue_label !=0)
  	    {
    		for(j=1;j<100;j++)
     		   {
       	              if( pixelValue_label == array_label[j-1])
         	         {
               		    k=k+1;
         		 }
     		    }
   		if(k==0 )
        	 {
                    array_label[i]=pixelValue_label;
                   // std::cout << array_label[i]<< std::endl;
                    i=i+1;
         	  }

      
      }
   ++imageIterator_2;
  }


// calculate sdt
unsigned int m;
for(m=0;array_label[m]>0;m++)
{
//get the single label region
itk::ImageRegionIterator<ImageType> imageIterator_3(m_FixedImageWarper->GetOutput(),m_FixedImageWarper->GetOutput()->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_4(image_label,image_label->GetLargestPossibleRegion());
while(!imageIterator_3.IsAtEnd())
  {    
      unsigned int pixelValue = imageIterator_3.Get();
      if (pixelValue == array_label[m]) 
        {
        imageIterator_4.Set(0);
        }
      else
        {
	  imageIterator_4.Set(1);
	}
    ++imageIterator_3;
    ++imageIterator_4;
  }
 // sdt the single label region
typename SignedMaurerDistanceMapImageFilterType::Pointer distanceMapImageFilter = SignedMaurerDistanceMapImageFilterType::New();
distanceMapImageFilter->SetInput(image_label);

//get the minimum and maximum in the single label sdm

typename ImageType::Pointer image_sdt = ImageType::New();
image_sdt = distanceMapImageFilter->GetOutput();// because the array_label have a lot zeros label
image_sdt->Update();

float minimum=0, maximum, pixelValue_am;
bool flag1=true;
itk::ImageRegionIterator<ImageType> imageIterator_5(image_sdt,image_sdt->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_6(image_label,image_label->GetLargestPossibleRegion());
while(!imageIterator_5.IsAtEnd())
  {
        pixelValue_am = imageIterator_5.Get();
           if (imageIterator_6.Get()==0)
           {
             if(flag1)
               {
                      flag1=false;
              //        minimum = pixelValue_am;
                      maximum = pixelValue_am;
                }
             if(pixelValue_am < minimum)
                {
              //    minimum = pixelValue_am;
                }
             else if(pixelValue_am > maximum)
                {
                  maximum = pixelValue_am;
                }
           }
    ++imageIterator_5;
    ++imageIterator_6;
    }

// norm the sdt single region and send it to totalsignedmap
float pixelValue_sdt, pixelValue_sdt_norm,diff;
diff=(maximum-minimum)+eps;  // may be zeros, should be think again;
//std::cout<< "fixed the diff is: "<<diff<<std::endl;
//std::cout<< "the minmum is: "<<minimum<<std::endl;
bool flag2; // for debug only
itk::ImageRegionIterator<ImageType> imageIterator_7(image_sdt,image_sdt->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_8(image_label,image_label->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_9(m_fixedImage,m_fixedImage->GetLargestPossibleRegion());
while(!imageIterator_7.IsAtEnd())
  {
        pixelValue_sdt = imageIterator_7.Get();
        if (imageIterator_8.Get()==0)
	{
	      pixelValue_sdt_norm = (pixelValue_sdt-minimum)/diff;
          pixelValue_sdt_norm = pixelValue_sdt_norm*255;
	      imageIterator_9.Set(pixelValue_sdt_norm);
	}
    ++imageIterator_7;
    ++imageIterator_8;
    ++imageIterator_9;
  }

}
typedef itk::ImageFileWriter< ImageType >  WriterType;
typename WriterType::Pointer writer1 = WriterType::New();
writer1->SetFileName( "test_warpedfixedsdmimage.nii.gz" );
writer1->SetInput( m_fixedImage.GetPointer() );
writer1->Update();
this->SetFixedSDMImage(m_fixedImage);

}


template <class TFixedImage, class TMovingImage, class TDeformationField>
void ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::SignedDistanceMap_orignalmovingImage()
{

float eps = 0.0001;
typedef float      PixelType;
typedef unsigned int maximum,minimum;
typedef itk::Image< PixelType, ImageDimension >    ImageType;

typedef  itk::SignedMaurerDistanceMapImageFilter< ImageType, ImageType  > SignedMaurerDistanceMapImageFilterType;

  typedef itk::ImageDuplicator< ImageType > DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(this->GetMovingImage());
  duplicator->Update();

  typename ImageType::Pointer image_label2 = duplicator->GetOutput();

  typename DuplicatorType::Pointer duplicator1 = DuplicatorType::New();
  duplicator1->SetInputImage(this->GetMovingImage());
  duplicator1->Update();
  typename ImageType::Pointer image_initial = duplicator1->GetOutput();

  typename ImageType::Pointer m_orignalmovingImage = ImageType::New();
  typename ImageType::Pointer image22 = duplicator->GetOutput();

   typedef itk::CastImageFilter< ImageType, ImageType > CastFilterType;
   typename CastFilterType::Pointer castFilter = CastFilterType::New();
   castFilter->SetInput(image22);
   castFilter->Update();
   m_orignalmovingImage = castFilter->GetOutput();

itk::ImageRegionIterator<ImageType> imageIterator_1(m_orignalmovingImage,m_orignalmovingImage->GetLargestPossibleRegion());
while(!imageIterator_1.IsAtEnd())
  {
     imageIterator_1.Set(0);
       ++imageIterator_1;
  }


//get the array of labels
unsigned int array_label[100]={0};
unsigned int i=0,j,pixelValue_label;
itk::ImageRegionIterator<ImageType> imageIterator_2(image_initial,image_initial->GetLargestPossibleRegion());
while(!imageIterator_2.IsAtEnd())
  {
        pixelValue_label = imageIterator_2.Get();
	  unsigned int k=0;
	  if (pixelValue_label !=0)
  	    {
    		for(j=1;j<100;j++)
     		   {
       	              if( pixelValue_label == array_label[j-1])
         	         {
               		    k=k+1;
         		 }
     		    }
   		if(k==0 )
        	 {
                    array_label[i]=pixelValue_label;
                   // std::cout << array_label[i]<< std::endl;
                    i=i+1;
         	  }

            }
   ++imageIterator_2;
  }


// calculate sdt
unsigned int m;
for(m=0;array_label[m]>0;m++)
{
//get the single label region
itk::ImageRegionIterator<ImageType> imageIterator_3(image_initial,image_initial->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_4(image_label2,image_label2->GetLargestPossibleRegion());
while(!imageIterator_3.IsAtEnd())
  {    
      unsigned int pixelValue = imageIterator_3.Get();
      if (pixelValue == array_label[m]) 
        {
        imageIterator_4.Set(0);
        }
      else
        {
	  imageIterator_4.Set(1);
	}
    ++imageIterator_3;
    ++imageIterator_4;
  }
 // sdt the single label region
typename SignedMaurerDistanceMapImageFilterType::Pointer distanceMapImageFilter = SignedMaurerDistanceMapImageFilterType::New();
distanceMapImageFilter->SetInput(image_label2);

//get the minimum and maximum in the single label sdm

typename ImageType::Pointer image_sdt2 = ImageType::New();
image_sdt2 = distanceMapImageFilter->GetOutput();// because the array_label have a lot zeros label
image_sdt2->Update();

float minimum=0, maximum, pixelValue_am;
bool flag1=true;
itk::ImageRegionIterator<ImageType> imageIterator_5(image_sdt2,image_sdt2->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_6(image_label2,image_label2->GetLargestPossibleRegion());
while(!imageIterator_5.IsAtEnd())
  {
        pixelValue_am = imageIterator_5.Get();
           if (imageIterator_6.Get()==0)
           {
             if(flag1)
               {
                      flag1=false;
              //        minimum = pixelValue_am;
                      maximum = pixelValue_am;
                }
             if(pixelValue_am < minimum)
                {
              //    minimum = pixelValue_am;
                }
             else if(pixelValue_am > maximum)
                {
                  maximum = pixelValue_am;
                }
           }
    ++imageIterator_5;
    ++imageIterator_6;
    }

// norm the sdt single region and send it to totalsignedmap
float pixelValue_sdt, pixelValue_sdt_norm,diff;
diff=(maximum-minimum)+eps;  // may be zeros, should be think again;
//std::cout<< "moving the diff is: "<<diff<<std::endl;
//std::cout<< "the minmum is: "<<minimum<<std::endl;
bool flag2; // for debug only
itk::ImageRegionIterator<ImageType> imageIterator_7(image_sdt2,image_sdt2->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_8(image_label2,image_label2->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_9(m_orignalmovingImage,m_orignalmovingImage->GetLargestPossibleRegion());
while(!imageIterator_7.IsAtEnd())
  {
        pixelValue_sdt = imageIterator_7.Get();
        if (imageIterator_8.Get()==0)
	{
	      pixelValue_sdt_norm = (pixelValue_sdt-minimum)/diff;
          pixelValue_sdt_norm = pixelValue_sdt_norm*255;
	      imageIterator_9.Set(pixelValue_sdt_norm);
	}
    ++imageIterator_7;
    ++imageIterator_8;
    ++imageIterator_9;
  }

}

this->SetorignalMovingSDMImage(m_orignalmovingImage);
} 
    

template <class TFixedImage, class TMovingImage, class TDeformationField>
void ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::SignedDistanceMap_orignalfixedImage()
{

float eps = 0.0001;
typedef float      PixelType;
typedef unsigned int maximum,minimum;
typedef itk::Image< PixelType, ImageDimension >    ImageType;
typedef  itk::SignedMaurerDistanceMapImageFilter< ImageType, ImageType  > SignedMaurerDistanceMapImageFilterType;

  typedef itk::ImageDuplicator< ImageType > DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(this->GetFixedImage());
  duplicator->Update();

//typedef itk::ImageFileWriter< ImageType >  WriterType;
//typename WriterType::Pointer writer1 = WriterType::New();
//writer1->SetFileName( "test_fixedimage.nii.gz" );
//writer1->SetInput( this->GetFixedImage() );
//writer1->Update();

  
  typename ImageType::Pointer image_label = duplicator->GetOutput();

  typename DuplicatorType::Pointer duplicator2 = DuplicatorType::New();
  duplicator2->SetInputImage(this->GetFixedImage());
  duplicator2->Update();
typename ImageType::Pointer image_initial = duplicator2->GetOutput();
  typename ImageType::Pointer m_orignalfixedImage = ImageType::New();
  typename ImageType::Pointer image2 = duplicator->GetOutput();

   typedef itk::CastImageFilter< ImageType, ImageType > CastFilterType;
   typename CastFilterType::Pointer castFilter = CastFilterType::New();
   castFilter->SetInput(image2);
   castFilter->Update();
   m_orignalfixedImage = castFilter->GetOutput();

itk::ImageRegionIterator<ImageType> imageIterator_1(m_orignalfixedImage,m_orignalfixedImage->GetLargestPossibleRegion());
while(!imageIterator_1.IsAtEnd())
  {
     imageIterator_1.Set(0);
       ++imageIterator_1;
  }


//get the array of labels
unsigned int array_label[100]={0};
unsigned int i=0,j,pixelValue_label;
itk::ImageRegionIterator<ImageType> imageIterator_2(image_initial,image_initial->GetLargestPossibleRegion());
while(!imageIterator_2.IsAtEnd())
  {
        pixelValue_label = imageIterator_2.Get();
	  unsigned int k=0;
	  if (pixelValue_label !=0)
  	    {
    		for(j=1;j<100;j++)
     		   {
       	              if( pixelValue_label == array_label[j-1])
         	         {
               		    k=k+1;
         		 }
     		    }
   		if(k==0 )
        	 {
                    array_label[i]=pixelValue_label;
                   // std::cout << array_label[i]<< std::endl;
                    i=i+1;
         	  }

            }
   ++imageIterator_2;
  }


// calculate sdt
unsigned int m;
for(m=0;array_label[m]>0;m++)
{
//get the single label region
itk::ImageRegionIterator<ImageType> imageIterator_3(image_initial,image_initial->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_4(image_label,image_label->GetLargestPossibleRegion());
while(!imageIterator_3.IsAtEnd())
  {    
      unsigned int pixelValue = imageIterator_3.Get();
      if (pixelValue == array_label[m]) 
        {
        imageIterator_4.Set(0);
        }
      else
        {
	  imageIterator_4.Set(1);
	}
    ++imageIterator_3;
    ++imageIterator_4;
  }
 // sdt the single label region
typename SignedMaurerDistanceMapImageFilterType::Pointer distanceMapImageFilter = SignedMaurerDistanceMapImageFilterType::New();
distanceMapImageFilter->SetInput(image_label);

//get the minimum and maximum in the single label sdm

typename ImageType::Pointer image_sdt = ImageType::New();
image_sdt = distanceMapImageFilter->GetOutput();// because the array_label have a lot zeros label
image_sdt->Update();

float minimum=0, maximum, pixelValue_am;
bool flag1=true;
itk::ImageRegionIterator<ImageType> imageIterator_5(image_sdt,image_sdt->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_6(image_label,image_label->GetLargestPossibleRegion());
while(!imageIterator_5.IsAtEnd())
  {
        pixelValue_am = imageIterator_5.Get();
           if (imageIterator_6.Get()==0)
           {
             if(flag1)
               {
                      flag1=false;
              //        minimum = pixelValue_am;
                      maximum = pixelValue_am;
                }
             if(pixelValue_am < minimum)
                {
              //    minimum = pixelValue_am;
                }
             else if(pixelValue_am > maximum)
                {
                  maximum = pixelValue_am;
                }
           }
    ++imageIterator_5;
    ++imageIterator_6;
    }

// norm the sdt single region and send it to totalsignedmap
float pixelValue_sdt, pixelValue_sdt_norm,diff;
diff=(maximum-minimum)+eps;  // may be zeros, should be think again;
//std::cout<< "moving the diff is: "<<diff<<std::endl;
//std::cout<< "the minmum is: "<<minimum<<std::endl;
bool flag2; // for debug only
itk::ImageRegionIterator<ImageType> imageIterator_7(image_sdt,image_sdt->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_8(image_label,image_label->GetLargestPossibleRegion());
itk::ImageRegionIterator<ImageType> imageIterator_9(m_orignalfixedImage,m_orignalfixedImage->GetLargestPossibleRegion());
while(!imageIterator_7.IsAtEnd())
  {
        pixelValue_sdt = imageIterator_7.Get();
        if (imageIterator_8.Get()==0)
	{
	      pixelValue_sdt_norm = (pixelValue_sdt-minimum)/diff;
          pixelValue_sdt_norm = pixelValue_sdt_norm*255;
	      imageIterator_9.Set(pixelValue_sdt_norm);
	}
    ++imageIterator_7;
    ++imageIterator_8;
    ++imageIterator_9;
  }

}
//typedef itk::ImageFileWriter< ImageType >  WriterType;
//typename WriterType::Pointer writer = WriterType::New();
//writer->SetFileName( "test_fixedsdmimage.nii.gz" );
//writer->SetInput( m_orignalfixedImage.GetPointer() );
//writer->Update();


this->SetorignalFixedSDMImage(m_orignalfixedImage);
} 
    

/**
 * Update the metric and release the per-thread-global data.
 */
    template <class TFixedImage, class TMovingImage, class TDeformationField>
    void
    ESMInvConDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
    ::ReleaseGlobalDataPointer( void *gd ) const
{
    GlobalDataStruct * globalData = (GlobalDataStruct *) gd;
    
    m_MetricCalculationLock.Lock();
    m_SumOfSquaredDifference += globalData->m_SumOfSquaredDifference;
    m_NumberOfPixelsProcessed += globalData->m_NumberOfPixelsProcessed;
    m_SumOfSquaredChange += globalData->m_SumOfSquaredChange;
    if( m_NumberOfPixelsProcessed )
    {
        m_Metric = m_SumOfSquaredDifference /
        static_cast<double>( m_NumberOfPixelsProcessed );
        m_RMSChange = vcl_sqrt( m_SumOfSquaredChange /
        static_cast<double>( m_NumberOfPixelsProcessed ) );
    }
    m_MetricCalculationLock.Unlock();
    
    delete globalData;
    }
    
} // end namespace itk

#endif
