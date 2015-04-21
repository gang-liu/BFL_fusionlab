//#include "fastcopy.h"

#include <itkWarpImageFilter.h>

//#include <boost/timer.hpp>

#include <mex.h>

template <class MatlabPixelType, unsigned int Dimension>
void warpimage(int nlhs,
               mxArray *plhs[],
               int nrhs,
               const mxArray *prhs[])
{
   typedef float PixelType;
   typedef itk::Image< PixelType, Dimension >           ImageType;
   
   typedef float                                        VectorComponentType;
   typedef itk::Vector<VectorComponentType, Dimension>  VectorPixelType;
   typedef itk::Image<VectorPixelType, Dimension>       DeformationFieldType;

   typedef itk::WarpImageFilter
      <ImageType, ImageType, DeformationFieldType>      WarperType;

   
   //boost::timer timer;


   // Allocate images and deformation field
   typename ImageType::Pointer image
      = ImageType::New();
   typename DeformationFieldType::Pointer field
      = DeformationFieldType::New();
   
   typename DeformationFieldType::SpacingType spacing;
   spacing.Fill( 1.0 );
   
   typename DeformationFieldType::PointType origin;
   origin.Fill( 0.0 );
   
   typename DeformationFieldType::RegionType     region;
   typename DeformationFieldType::SizeType       size;
   typename DeformationFieldType::IndexType      start;

   unsigned int numPix(1u);
   const MatlabPixelType * iminptr =  static_cast<const MatlabPixelType *>(mxGetData(prhs[0]));
   const MatlabPixelType * fieldinptrs[Dimension];
   mwSize matlabdims[Dimension];
   for (unsigned int d=0; d<Dimension; d++)
   {
      matlabdims[d]= mxGetDimensions(prhs[0])[d];
      size[d] = matlabdims[d];
      start[d] = 0;
      numPix *= size[d];

      fieldinptrs[d] = static_cast<const MatlabPixelType *>(mxGetData(prhs[1+d]));
   }
   
   region.SetSize( size );
   region.SetIndex( start );

   image->SetOrigin( origin );
   image->SetSpacing( spacing );
   image->SetRegions( region );
   image->Allocate();
   
   field->SetOrigin( origin );
   field->SetSpacing( spacing );
   field->SetRegions( region );
   field->Allocate();


   //mexPrintf("done Allocate(); %f sec\n", timer.elapsed());
   //timer.restart();

   
   PixelType * imptr = image->GetBufferPointer();
   const PixelType * const buff_end = imptr + numPix;
   VectorPixelType * fieldptr = field->GetBufferPointer();
   
   while ( imptr != buff_end )
   {
      *imptr++ = *iminptr++;
         
      for (unsigned int d=0; d<Dimension; d++)
      {
         (*fieldptr)[d] = *(fieldinptrs[d])++;
      }
      ++fieldptr;
   }

   //mz::writeRaw<ImageType>(fixedimage,"fixedimage.mha");
   //mz::writeRaw<ImageType>(movingimage,"movingimage.mha");
   //mexPrintf("done inputs copy %f sec\n", timer.elapsed());
   //timer.restart();
   
   
   // Warp the image
   typename WarperType::Pointer warper = WarperType::New();
   warper->SetInput( image );
   warper->SetOutputSpacing( spacing );
   warper->SetOutputOrigin( origin );
   warper->SetDeformationField( field );
   if ( std::numeric_limits<PixelType>::has_quiet_NaN )
   {
      warper->SetEdgePaddingValue( std::numeric_limits<PixelType>::quiet_NaN() );
   }

   warper->UpdateLargestPossibleRegion();

   //mexPrintf("done warper->UpdateLargestPossibleRegion(); %f sec\n", timer.elapsed());
   //timer.restart();


   // Allocate output
   const mxClassID classID = mxGetClassID(prhs[0]);
   plhs[0] = mxCreateNumericArray(
      Dimension, matlabdims, classID, mxREAL);
   MatlabPixelType * outptr = static_cast<MatlabPixelType *>(mxGetData(plhs[0]));
   
   //mexPrintf("done allocate output %f sec\n", timer.elapsed());
   //timer.restart();
   

   // put result into output
   const PixelType * warpptr = warper->GetOutput()->GetBufferPointer();
   //mz::copy(warpptr, warpptr + numPix, outptr);
   std::copy(warpptr, warpptr + numPix, outptr);

   //mexPrintf("done outputs copy %f sec\n", timer.elapsed());
}


void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   if (nrhs!=3 and nrhs!=4)
   {
      mexErrMsgTxt("3 or 4 inputs required.");
   }

   const int dim=nrhs-1;

   const mxClassID classID = mxGetClassID(prhs[0]);

   /* The inputs must be noncomplex double matrices.*/
   for (int n=0; n<nrhs; n++)
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
            warpimage<float,2>(nlhs, plhs, nrhs, prhs);
            break;
         case mxDOUBLE_CLASS:
            warpimage<double,2>(nlhs, plhs, nrhs, prhs);
            break;
         default:
            mexErrMsgTxt("Pixel type unsupported.");
      }
      break;
   case 3:
      switch ( classID )
      {
         case mxSINGLE_CLASS:    
            warpimage<float,3>(nlhs, plhs, nrhs, prhs);
            break;
         case mxDOUBLE_CLASS:
            warpimage<double,3>(nlhs, plhs, nrhs, prhs);
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


