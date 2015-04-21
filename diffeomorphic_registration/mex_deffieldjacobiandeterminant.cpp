//#include "fastcopy.h"

#include <itkDisplacementFieldJacobianDeterminantFilter.h>

//#include <boost/timer.hpp>

#include <dlfcn.h>

#include <mex.h>

template <class MatlabPixelType, unsigned int Dimension>
void deffieldjacobiandeterminant(int nlhs,
                                 mxArray *plhs[],
                                 int nrhs,
                                 const mxArray *prhs[])
{
   typedef float                                        VectorComponentType;
   typedef itk::Vector<VectorComponentType, Dimension>  VectorPixelType;
   typedef itk::Image<VectorPixelType, Dimension>       DeformationFieldType;

   typedef float PixelType;
   typedef itk::Image< PixelType, Dimension >           ImageType;

   typedef itk::DisplacementFieldJacobianDeterminantFilter
      <DeformationFieldType, PixelType>                 JacobianDeterminantFilterType;

   //boost::timer timer;

   // Allocate deformation field
   typename DeformationFieldType::Pointer field =
      DeformationFieldType::New();
   
   typename DeformationFieldType::SpacingType spacing;
   spacing.Fill( 1.0 );
   
   typename DeformationFieldType::PointType origin;
   origin.Fill( 0.0 );
   
   typename DeformationFieldType::RegionType     region;
   typename DeformationFieldType::SizeType       size;
   typename DeformationFieldType::IndexType      start;

   unsigned int numPix(1u);
   const MatlabPixelType * inptrs[Dimension];
   mwSize matlabdims[Dimension];
   for (unsigned int d=0; d<Dimension; d++)
   {
      matlabdims[d]= mxGetDimensions(prhs[0])[d];
      size[d] = matlabdims[d];
      start[d] = 0;
      numPix *= size[d];

      inptrs[d] = static_cast<const MatlabPixelType *>(mxGetData(prhs[d]));
   }
   
   region.SetSize( size );
   region.SetIndex( start );
   
   field->SetOrigin( origin );
   field->SetSpacing( spacing );
   field->SetRegions( region );
   field->Allocate();

   //mexPrintf("done field->Allocate(); %f sec\n", timer.elapsed());
   //timer.restart();

   VectorPixelType * ptr = field->GetBufferPointer();
   const VectorPixelType * const buff_end = ptr + numPix;
   
   while ( ptr != buff_end )
   {
      for (unsigned int d=0; d<Dimension; d++)
      {
         (*ptr)[d] = *(inptrs[d])++;
      }
      ++ptr;
   }

   //mexPrintf("done inputs copy %f sec\n", timer.elapsed());
   //timer.restart();
   
   // Compute jacobian determinant
   typename JacobianDeterminantFilterType::Pointer
      jacdetfilter = JacobianDeterminantFilterType::New();

   jacdetfilter->SetInput( field );
   jacdetfilter->SetUseImageSpacing( false );

   jacdetfilter->UpdateLargestPossibleRegion();

   //mexPrintf("done jacobianFilter->UpdateLargestPossibleRegion; %f sec\n", timer.elapsed());
   //timer.restart();

   const mxClassID classID = mxGetClassID(prhs[0]);
   plhs[0] = mxCreateNumericArray(
      Dimension, matlabdims, classID, mxREAL);
   MatlabPixelType * outptr = static_cast<MatlabPixelType *>(mxGetData(plhs[0]));

   //mexPrintf("done allocate output %f sec\n", timer.elapsed());
   //timer.restart();
   

   // put result into output
   const PixelType * jacdetptr = jacdetfilter->GetOutput()->GetBufferPointer();
   //mz::copy(jacdetptr, jacdetptr + numPix, outptr);
   std::copy(jacdetptr, jacdetptr + numPix, outptr);

   //mexPrintf("done output copy %f sec\n", timer.elapsed());
}


void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
   dlopen("libITKAlgorithms.so", RTLD_LAZY|RTLD_GLOBAL);
   dlopen("libITKBasicFilters.so", RTLD_LAZY|RTLD_GLOBAL);
   dlopen("libITKCommon.so", RTLD_LAZY|RTLD_GLOBAL);
   dlopen("libITKIO.so", RTLD_LAZY|RTLD_GLOBAL);
   dlopen("libITKNumerics.so", RTLD_LAZY|RTLD_GLOBAL);
   dlopen("libMosaicing.so", RTLD_LAZY|RTLD_GLOBAL);
   dlopen("libCalibration.so", RTLD_LAZY|RTLD_GLOBAL);
   dlopen("libDataIO.so", RTLD_LAZY|RTLD_GLOBAL);
   dlopen("libCommon.so", RTLD_LAZY|RTLD_GLOBAL);
   dlopen("libExporter.so", RTLD_LAZY|RTLD_GLOBAL);
   dlopen("libItkProcessing.so", RTLD_LAZY|RTLD_GLOBAL);
   dlopen("libProcessing.so", RTLD_LAZY|RTLD_GLOBAL);
   dlopen("libRobustEstimation.so", RTLD_LAZY|RTLD_GLOBAL);
   
   /* Check for proper number of arguments. */
   if (nrhs<2 or nrhs>3)
   {
      mexErrMsgTxt("Two or three inputs required.");
   }

   const int dim=nrhs;

   const mxClassID classID = mxGetClassID(prhs[0]);

   /* The inputs must be noncomplex double matrices.*/
   for (int d=0; d<dim; d++)
   {
      if ( mxGetClassID(prhs[d])!=classID || mxIsComplex(prhs[d]) )
      {
         mexErrMsgTxt("Input must be a noncomplex floating point.");
      }

      if ( mxGetNumberOfDimensions(prhs[d]) != dim )
      {
         mexErrMsgTxt("The dimension of the inputs must match the number of inputs.");
      }

      for (int dd=0; dd<dim; dd++)
      {
         if ( mxGetDimensions(prhs[d])[dd] != mxGetDimensions(prhs[0])[dd] )
         {
            mexErrMsgTxt("Inputs must have the same size.");
         }
      }
   }

   if (nlhs != 1)
   {
      mexErrMsgTxt("Number of outputs must be one.");
   }

   switch ( dim )
   {
   case 2:
      switch ( classID )
      {
         case mxSINGLE_CLASS:    
            deffieldjacobiandeterminant<float,2>(nlhs, plhs, nrhs, prhs);
            break;
         case mxDOUBLE_CLASS:
            deffieldjacobiandeterminant<double,2>(nlhs, plhs, nrhs, prhs);
            break;
         default:
            mexErrMsgTxt("Pixel type unsupported.");
      }
      break;
   case 3:
      switch ( classID )
      {
         case mxSINGLE_CLASS:    
            deffieldjacobiandeterminant<float,3>(nlhs, plhs, nrhs, prhs);
            break;
         case mxDOUBLE_CLASS:
            deffieldjacobiandeterminant<double,3>(nlhs, plhs, nrhs, prhs);
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


