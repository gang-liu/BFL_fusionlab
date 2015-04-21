#include "itkESMInvConDemonsRegistrationFunction.h"
//#include "fastimagewriting.h"

#include <itkNeighborhoodAlgorithm.h>

//#include <boost/timer.hpp>

#include <mex.h>

template <class MatlabPixelType, unsigned int Dimension>
void invcondemonsforces(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{  //typedef unsigned int Dimension=3;
   typedef float PixelType;
   typedef itk::Image< PixelType, Dimension >           ImageType;
   
   typedef float                                        VectorComponentType;
   typedef itk::Vector<VectorComponentType, Dimension>  VectorPixelType;
   typedef itk::Image<VectorPixelType, Dimension>       DeformationFieldType;

   typedef itk::ESMInvConDemonsRegistrationFunction
      <ImageType,ImageType,DeformationFieldType>        DemonsRegistrationFunctionType;

   
   //boost::timer timer;


   // Allocate images and deformation field
   typename ImageType::Pointer fixedimage
      = ImageType::New();
   typename ImageType::Pointer movingimage
      = ImageType::New();
   typename ImageType::Pointer jacobianimage
      = ImageType::New();  
   typename DeformationFieldType::Pointer field
      = DeformationFieldType::New();
      
   typename DeformationFieldType::Pointer inv_field
      = DeformationFieldType::New();
      
   typename DeformationFieldType::Pointer update
      = DeformationFieldType::New();
   
   typename DeformationFieldType::SpacingType spacing;
   spacing.Fill( 1.0 );
   
   typename DeformationFieldType::PointType origin;
   origin.Fill( 0.0 );
   
   typename DeformationFieldType::RegionType     region;
   typename DeformationFieldType::SizeType       size;
   typename DeformationFieldType::IndexType      start;

   unsigned int numPix(1u);
   const MatlabPixelType * fixinptr =  static_cast<const MatlabPixelType *>(mxGetData(prhs[0]));
   const MatlabPixelType * movinptr =  static_cast<const MatlabPixelType *>(mxGetData(prhs[1]));
   const MatlabPixelType * fieldinptrs[Dimension];
   const MatlabPixelType * inv_fieldinptrs[Dimension];
   
   
   mwSize matlabdims[Dimension];
   for (unsigned int d=0; d<Dimension; d++)
   {
      matlabdims[d]= mxGetDimensions(prhs[0])[d];
      size[d] = matlabdims[d];
      start[d] = 0;
      numPix *= size[d];

      fieldinptrs[d] = static_cast<const MatlabPixelType *>(mxGetData(prhs[2+d]));
      inv_fieldinptrs[d] = static_cast<const MatlabPixelType *>(mxGetData(prhs[2+Dimension+d]));
   }
   const MatlabPixelType * jacobianptr = static_cast<const MatlabPixelType *> (mxGetData(prhs[2*Dimension + 2]));
   
   const unsigned int UseJacFlag = static_cast<unsigned int>( mxGetPr(prhs[2*Dimension+3])[0] );
   const double RegWeight = static_cast<double>( mxGetPr(prhs[2*Dimension+4])[0] );
   
   //std::cout << "RegWeight : " << RegWeight << std::endl;
   
   region.SetSize( size );
   region.SetIndex( start );

   fixedimage->SetOrigin( origin );
   fixedimage->SetSpacing( spacing );
   fixedimage->SetRegions( region );
   fixedimage->Allocate();

   movingimage->SetOrigin( origin );
   movingimage->SetSpacing( spacing );
   movingimage->SetRegions( region );
   movingimage->Allocate();
   
   if (UseJacFlag > 0)
   {
    jacobianimage->SetOrigin( origin );
    jacobianimage->SetSpacing( spacing );
    jacobianimage->SetRegions( region );
    jacobianimage->Allocate();
   }
   
   field->SetOrigin( origin );
   field->SetSpacing( spacing );
   field->SetRegions( region );
   field->Allocate();
   
   inv_field->SetOrigin( origin );
   inv_field->SetSpacing( spacing );
   inv_field->SetRegions( region );
   inv_field->Allocate();

   update->SetOrigin( origin );
   update->SetSpacing( spacing );
   update->SetRegions( region );
   update->Allocate();


   //mexPrintf("done Allocate(); %f sec\n");
   //timer.restart();

   
   PixelType * fixptr = fixedimage->GetBufferPointer();
   const PixelType * const fixbuff_end = fixptr + numPix;
   PixelType * movptr = movingimage->GetBufferPointer();
   PixelType * jacptr = NULL;
   if (UseJacFlag > 0)
   {
       jacptr = jacobianimage->GetBufferPointer();
   }
   VectorPixelType * fieldptr = field->GetBufferPointer();
   VectorPixelType * inv_fieldptr = inv_field->GetBufferPointer();
   
   while ( fixptr != fixbuff_end )
   {
      *fixptr++ = *fixinptr++;
      *movptr++ = *movinptr++;
      if (UseJacFlag > 0)
      {
        *jacptr++ = *jacobianptr++; 
      }
      for (unsigned int d=0; d<Dimension; d++)
      {
         (*fieldptr)[d] = *(fieldinptrs[d])++;
      }
      ++fieldptr;
      
      for (unsigned int d=0; d<Dimension; d++)
      {
         (*inv_fieldptr)[d] = *(inv_fieldinptrs[d])++;
      }
      ++inv_fieldptr;
   }

   // Create demons function
   typename DemonsRegistrationFunctionType::Pointer drfp
      = DemonsRegistrationFunctionType::New();

   //mexPrintf("step size: %f\n",mxGetPr(prhs[2*Dimension+2])[0]);
   //drfp->SetMaximumUpdateStepLength( mxGetPr(prhs[2*Dimension+2])[0] );

   typename DemonsRegistrationFunctionType::GradientType gtype = DemonsRegistrationFunctionType::Symmetric;
   drfp->SetUseGradientType( gtype );

   drfp->SetDeformationField( field );
   drfp->SetInvDeformationField( inv_field);
   drfp->SetFixedImage( fixedimage );
   drfp->SetMovingImage( movingimage );
   
   drfp->SetRegWeight(RegWeight);
   
   if (UseJacFlag > 0)
   {
       drfp->SetUseJacobian(true);
       drfp->SetJacobianDetImage(jacobianimage);
   }
   else
   {
       drfp->SetUseJacobian(false);
   }
   
   drfp->InitializeIteration();

   //mexPrintf("done demons function init %f sec\n", timer.elapsed());
   //timer.restart();

   const itk::Size<Dimension> radius = drfp->GetRadius();

   // Break the input into a series of regions.  The first region is free
   // of boundary conditions, the rest with boundary conditions.  We operate
   // on the output region because input has been copied to output.
   typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
      <DeformationFieldType> FaceCalculatorType;
   typedef typename FaceCalculatorType::FaceListType FaceListType;
   typedef typename DemonsRegistrationFunctionType::NeighborhoodType
      NeighborhoodIteratorType;
   typedef itk::ImageRegionIterator<DeformationFieldType> UpdateIteratorType;

   FaceCalculatorType faceCalculator;
   
   FaceListType faceList = faceCalculator(field, region, radius);
   typename FaceListType::iterator fIt = faceList.begin();
   
   // Ask the function object for a pointer to a data structure it
   // will use to manage any global values it needs.  We'll pass this
   // back to the function object at each calculation and then
   // again so that the function object can use it to determine a
   // time step for this iteration.
   void * globalData = drfp->GetGlobalDataPointer();
   
   // Process the non-boundary region.
   NeighborhoodIteratorType nD(radius, field, *fIt);
   UpdateIteratorType       nU(update,  *fIt);
   nD.GoToBegin();
   while( !nD.IsAtEnd() )
   {
      nU.Value() = drfp->ComputeUpdate(nD, globalData);
      ++nD;
      ++nU;
   }

   // Process each of the boundary faces.
   
   NeighborhoodIteratorType bD;
   UpdateIteratorType   bU;
   for (++fIt; fIt != faceList.end(); ++fIt)
   {
      bD = NeighborhoodIteratorType(radius, field, *fIt);
      bU = UpdateIteratorType(update, *fIt);
      
      bD.GoToBegin();
      bU.GoToBegin();
      while ( !bD.IsAtEnd() )
      {
         bU.Value() = drfp->ComputeUpdate(bD, globalData);
         ++bD;
         ++bU;
      }
   }

   // Ask the finite difference function to compute the time step for
   // this iteration.  We give it the global data pointer to use, then
   // ask it to free the global data memory.
   //timeStep = df->ComputeGlobalTimeStep(globalData);
   drfp->ReleaseGlobalDataPointer(globalData);


   //mexPrintf("done actual computations %f sec\n", timer.elapsed());
   //timer.restart();


   // Allocate outputs
   const mxClassID classID = mxGetClassID(prhs[0]);
   MatlabPixelType * outptrs[Dimension];
   for (unsigned int d=0; d<Dimension; d++)
   {
      plhs[d] = mxCreateNumericArray(
         Dimension, matlabdims, classID, mxREAL);

      outptrs[d] = static_cast<MatlabPixelType *>(mxGetData(plhs[d]));
   }

   
   //mexPrintf("done allocate outputs %f sec\n", timer.elapsed());
   //timer.restart();
   

   // put result into outputs
   const VectorPixelType * upptr = update->GetBufferPointer();
   const VectorPixelType * const upbuff_end = upptr + numPix;
   
   while ( upptr != upbuff_end )
   {
      for (unsigned int d=0; d<Dimension; d++)
      {
         *(outptrs[d])++ = (*upptr)[d];
      }
      ++upptr;
   }

   //mexPrintf("done outputs copy %f sec\n", timer.elapsed());
}


void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   if (nrhs!=9 and nrhs!=11)
   {
      mexErrMsgTxt("9 or 11 inputs required.");
   }

   const int dim=(nrhs-5)/2;
   //mexPrintf("Dimension of images: %i\n",dim);
   const mxClassID classID = mxGetClassID(prhs[0]);

   /* The first inputs must be noncomplex double matrices.*/
   for (int n=0; n<2*dim+2; n++)
   {
      if ( mxGetClassID(prhs[n])!=classID || mxIsComplex(prhs[n]) )
      {
         mexErrMsgTxt("Input must be a noncomplex floating point.");
      }

      if ( mxGetNumberOfDimensions(prhs[n]) != dim )
      {
         mexErrMsgTxt("The dimension of the inputs must agree with the number of inputs.");
      }

      for (int dd=0; dd<dim; dd++)
      {
         if ( mxGetDimensions(prhs[n])[dd] != mxGetDimensions(prhs[0])[dd] )
         {
            mexErrMsgTxt("Inputs must have the same size.");
         }
      }
   }
    
   /* Check last two sinput */
   for (int n=2*dim+3; n<2*dim+5; n++)
   {
      if ( !mxIsDouble(prhs[n]) || mxIsComplex(prhs[n]) )
      {
         mexErrMsgTxt("Two last inputs must be double.");
      }
      if ( mxGetM(prhs[n]) != 1 || mxGetN(prhs[n]) != 1 )
      {
         mexErrMsgTxt("The dimension of the last input must be one.");
      }
   }
   
  
   if (static_cast<unsigned int>( mxGetPr(prhs[2*dim+3])[0] )  > 0)
   {
       int ii = 2*dim + 2;
       
      if ( mxGetClassID(prhs[ii])!=classID || mxIsComplex(prhs[ii]) )
      {
         mexErrMsgTxt("Input must be a noncomplex floating point.");
      }

      if ( mxGetNumberOfDimensions(prhs[ii]) != dim )
      {
         mexErrMsgTxt("The dimension of the inputs must agree with the number of inputs.");
      }

      for (int dd=0; dd<dim; dd++)
      {
         if ( mxGetDimensions(prhs[ii])[dd] != mxGetDimensions(prhs[0])[dd] )
         {
            mexErrMsgTxt("Inputs must have the same size.");
         }
      }
       
   }
   
   if (nlhs != dim)
   {
      mexErrMsgTxt("Number of outputs must agree with the number of inputs.");
   }
   //unsigned int dim=3;
   switch ( dim )
   {
   case 2:
      switch ( classID )
      {
         case mxSINGLE_CLASS:    
            invcondemonsforces<float,2>(nlhs, plhs, nrhs, prhs);
            break;
         case mxDOUBLE_CLASS:
            invcondemonsforces<double,2>(nlhs, plhs, nrhs, prhs);
            break;
         default:
            mexErrMsgTxt("Pixel type unsupported.");
      }
      break;
  case 3:
      switch ( classID )
      {
         case mxSINGLE_CLASS:    
            invcondemonsforces<float,3>(nlhs, plhs, nrhs, prhs);
            break;
         case mxDOUBLE_CLASS:
            invcondemonsforces<double,3>(nlhs, plhs, nrhs, prhs);
            break;
         default:
            mexErrMsgTxt("Pixel type unsupported.");
      }
      break;
   default:
     mexErrMsgTxt("Dimension unsupported.");
   }

   return;
}


