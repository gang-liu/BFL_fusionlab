//#include "fastcopy.h"
#include "itkImageFileReader.h"
#include "itkCommand.h"
#include "itkImage.h"

#include <itkWarpImageFilter.h>

//#include <boost/timer.hpp>

#include <mex.h>

template <class MatlabPixelType, unsigned int Dimension>
void readwarpfile(int nlhs,
               mxArray *plhs[],
               int nrhs,
               const mxArray *prhs[])
{
   typedef int PixelType;
   typedef itk::Image< PixelType, Dimension >           ImageType;
 //  mexPrintf("the dimentsion is ", Dimension);   
 //  typedef                                           PixelType;
   typedef itk::Vector<PixelType, Dimension>  VectorPixelType;
   typedef itk::Image<VectorPixelType, Dimension>       DeformationFieldType;

   typedef itk::WarpImageFilter
      <ImageType, ImageType, DeformationFieldType>      WarperType;
   mexPrintf("input number is: %d", nrhs);
   mexPrintf("output parameternumber is: %d", nlhs);
 // mexPrintf("the dimentsion is :%d", Dimension); 
   //boost::timer timer;
   typename DeformationFieldType::SizeType       size;
   typename DeformationFieldType::IndexType      start;
   unsigned int numPix(1u);
   char * inputDfFilename; 
   mwSize matlabdims[Dimension];

// here, we should read image parameter
   for (unsigned int d=0; d<Dimension; d++)
   {
      matlabdims[d]= mxGetDimensions(prhs[1])[d];
      size[d] = matlabdims[d];
      start[d] = 0;
      numPix *= size[d];

//    inptrs[d] = static_cast<const MatlabPixelType *>(mxGetData(prhs[d]));
   }
   inputDfFilename = mxArrayToString(prhs[0]);
   mexPrintf("output file name is %s", inputDfFilename);
   typename itk::ImageFileReader<DeformationFieldType>::Pointer df_reader = itk::ImageFileReader<DeformationFieldType>::New();
   df_reader->SetFileName(inputDfFilename);
   df_reader->Update();

   // Allocate images and deformation field
   typename ImageType::Pointer image
      = ImageType::New();
   typename DeformationFieldType::Pointer field = DeformationFieldType::New();
//  typename DeformationFieldType::Pointer field = 0; 
   typename DeformationFieldType::SpacingType spacing;
   spacing.Fill( 1.0 );
   
   typename DeformationFieldType::PointType origin;
   origin.Fill( 0.0 );
   
   typename DeformationFieldType::RegionType     region;
   field = df_reader->GetOutput();


   region = field->GetLargestPossibleRegion();

// initialize output point
   const mxClassID classID = mxGetClassID(prhs[1]);
   MatlabPixelType * fieldoutptrs[Dimension];
   for (unsigned int d=0; d<Dimension; d++)
   {
      plhs[d] = mxCreateNumericArray(
         Dimension, matlabdims, classID, mxREAL);

      fieldoutptrs[d] = static_cast<MatlabPixelType *>(mxGetData(plhs[d]));
   }


// read field
   VectorPixelType * fieldptr = field->GetBufferPointer();
   const VectorPixelType * const buff_end = fieldptr + numPix;
 

     while ( fieldptr != buff_end )
   {
      for (unsigned int d=0; d<Dimension; d++)
      {
         *(fieldoutptrs[d])++ = (*fieldptr)[d];
      }
      ++fieldptr;
   }

}


void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   if (nrhs!=2)
   {
      mexErrMsgTxt("2 input required.");
   }

   const int dim=nlhs;

   const mxClassID classID = mxGetClassID(prhs[1]);

   /* The inputs must be noncomplex double matrices.*/

   for (int n=1; n<nrhs-1; n++)
   {
      if ( mxGetClassID(prhs[n])!=classID || mxIsComplex(prhs[n]) )
      {
         mexErrMsgTxt("Input must be a noncomplex floating point.");
      }

//      if ( mxGetNumberOfDimensions(prhs[n]) != dim )
    //  {
  //       mexErrMsgTxt("The dimension of the  must agree with the number of inputs.");
    //  }

      for (int dd=0; dd<dim; dd++)
      {
         if ( mxGetDimensions(prhs[n])[dd] != mxGetDimensions(prhs[0])[dd] )
         {
            mexErrMsgTxt("Inputs must have the same size.");
         }
      }
   }

   if (nlhs != 3)
   {
      mexErrMsgTxt("Number of outputs must be 3.");
   }

   switch ( dim )
   {
   case 2:
      switch ( classID )
      {
         case mxSINGLE_CLASS:    
            readwarpfile<float,2>(nlhs, plhs, nrhs, prhs);
            break;
         case mxDOUBLE_CLASS:
            readwarpfile<double,2>(nlhs, plhs, nrhs, prhs);
            break;
         default:
            mexErrMsgTxt("Pixel type unsupported.");
      }
      break;
   case 3:
      switch ( classID )
      {
      //   std:cout<<"classID: "<<classID<<std::endl;
         case mxSINGLE_CLASS:    
            readwarpfile<float,3>(nlhs, plhs, nrhs, prhs);
            break;
         case mxDOUBLE_CLASS:
            readwarpfile<double,3>(nlhs, plhs, nrhs, prhs);
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


