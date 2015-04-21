//#include "fastcopy.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itkImage.h"

#include <itkWarpImageFilter.h>

//#include <boost/timer.hpp>

#include <mex.h>

template <class MatlabPixelType, unsigned int Dimension>
void writewarpfile(int nlhs,
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

 // mexPrintf("the dimentsion is :%d", Dimension); 
   //boost::timer timer;


   // Allocate images and deformation fieldnnnn
   typename ImageType::Pointer image
      = ImageType::New();
   typename DeformationFieldType::Pointer field = DeformationFieldType::New();
//  typename DeformationFieldType::Pointer field = 0; 
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
      
      fieldinptrs[d] = static_cast<const MatlabPixelType *>(mxGetData(prhs[d]));
   //   mexPrintf("the dimentsion is :%d", d);

   }
//   mexPrintf("er of pixel is:  %d \n", numPix);
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
//   const PixelType * const buff_begin = imptr;
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
   // here the deformation field is not like logdomain deformation field. why?
   //mz::writeRaw<ImageType>(fixedimage,"fixedimage.mha");
   //mz::writeRaw<ImageType>(movingimage,"movingimage.mha");
   //mexPrintf("done inputs copy %f sec\n", timer.elapsed());
   //timer.restart();
    char * outputDfFilename;
  
    outputDfFilename = mxArrayToString(prhs[3]); 
    mexPrintf("the output file name is: %s", outputDfFilename);
   //write image
   //
   //typename itk::ImageFileWriter<ImageType>::Pointer writer = itk::ImageFileWriter<ImageType>::New();
    //writer->SetFileName(outputDfFilename);
      //           writer->SetInput(image);
        //        std::cout<<"outputfile name: "<<outputDfFilename<<std::endl;
          //       writer->SetUseCompression(true);
      //           writer->Update();
  
   //write deformation field 
     
                typename itk::ImageFileWriter<DeformationFieldType>::Pointer df_writer = itk::ImageFileWriter<DeformationFieldType>::New();
                df_writer->SetFileName(outputDfFilename);
                df_writer->SetInput(field);
            
                df_writer->SetUseCompression(true);
           //     df_writer->UpdateLargestPossibleRegion();
               df_writer->Update();
/*               try
              {
                 df_writer->Update();
                        }
    catch( itk::ExceptionObject& err )
      {
      std::cout << "Unexpected error." << std::endl;
      std::cout << err << std::endl;
      exit( EXIT_FAILURE );
      }
*/
                mxFree(outputDfFilename);
}


void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   if (nrhs!=4 and nrhs!=3)
   {
      mexErrMsgTxt("4 or 3 inputs required.");
   }

   const int dim=nrhs-1;

   const mxClassID classID = mxGetClassID(prhs[0]);

   /* The inputs must be noncomplex double matrices.*/

   for (int n=0; n<nrhs-1; n++)
   {
      if ( mxGetClassID(prhs[n])!=classID || mxIsComplex(prhs[n]) )
      {
         mexErrMsgTxt("Input must be a noncomplex floating point.");
      }

      if ( mxGetNumberOfDimensions(prhs[n]) != dim )
      {
         mexErrMsgTxt("The dimension of the  must agree with the number of inputs.");
      }

      for (int dd=0; dd<dim; dd++)
      {
         if ( mxGetDimensions(prhs[n])[dd] != mxGetDimensions(prhs[0])[dd] )
         {
            mexErrMsgTxt("Inputs must have the same size.");
         }
      }
   }

   if (nlhs != 0)
   {
      mexErrMsgTxt("Number of outputs must be zero.");
   }

   switch ( dim )
   {
   case 2:
      switch ( classID )
      {
         case mxSINGLE_CLASS:    
            writewarpfile<float,2>(nlhs, plhs, nrhs, prhs);
            break;
         case mxDOUBLE_CLASS:
            writewarpfile<double,2>(nlhs, plhs, nrhs, prhs);
            break;
         default:
            mexErrMsgTxt("Pixel type unsupported.");
      }
      break;
   case 3:
      switch ( classID )
      {
         case mxSINGLE_CLASS:    
            writewarpfile<float,3>(nlhs, plhs, nrhs, prhs);
            break;
         case mxDOUBLE_CLASS:
            writewarpfile<double,3>(nlhs, plhs, nrhs, prhs);
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


