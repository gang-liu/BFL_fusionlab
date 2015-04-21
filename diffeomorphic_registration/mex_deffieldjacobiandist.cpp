#include "itkVectorCentralDifferenceImageFunction.h"

#include <itkImageRegionConstIteratorWithIndex.h>

//#include <boost/timer.hpp>

#include <mex.h>

template <class MatlabPixelType, unsigned int Dimension>
void deffieldjacobiandist(int nlhs,
                          mxArray *plhs[],
                          int nrhs,
                          const mxArray *prhs[])
{
   typedef float                                        VectorComponentType;
   typedef itk::Vector<VectorComponentType, Dimension>  VectorPixelType;
   typedef itk::Image<VectorPixelType, Dimension>       DeformationFieldType;

   typedef itk::VectorCentralDifferenceImageFunction
      <DeformationFieldType>                            WarpGradientCalculatorType;

   typedef typename WarpGradientCalculatorType::OutputType
                                                        WarpGradientType;

   typedef itk::ImageRegionConstIteratorWithIndex
      <DeformationFieldType>                            FieldIteratorType;

   //boost::timer timer;

   // Allocate deformation fields
   typename DeformationFieldType::Pointer leftfield =
      DeformationFieldType::New();
   typename DeformationFieldType::Pointer rightfield =
      DeformationFieldType::New();
   
   typename DeformationFieldType::SpacingType spacing;
   spacing.Fill( 1.0 );
   
   typename DeformationFieldType::PointType origin;
   origin.Fill( 0.0 );
   
   typename DeformationFieldType::RegionType     region;
   typename DeformationFieldType::SizeType       size;
   typename DeformationFieldType::IndexType      start;

   unsigned int numPix(1u);
   const MatlabPixelType * leftinptrs[Dimension];
   const MatlabPixelType * rightinptrs[Dimension];
   mwSize matlabdims[Dimension];
   for (unsigned int d=0; d<Dimension; d++)
   {
      matlabdims[d]= mxGetDimensions(prhs[0])[d];
      size[d] = matlabdims[d];
      start[d] = 0;
      numPix *= size[d];

      leftinptrs[d] = static_cast<const MatlabPixelType *>(mxGetData(prhs[d]));
      rightinptrs[d] = static_cast<const MatlabPixelType *>(mxGetData(prhs[Dimension+d]));
   }
   
   region.SetSize( size );
   region.SetIndex( start );
   
   leftfield->SetOrigin( origin );
   leftfield->SetSpacing( spacing );
   leftfield->SetRegions( region );
   leftfield->Allocate();

   rightfield->SetOrigin( origin );
   rightfield->SetSpacing( spacing );
   rightfield->SetRegions( region );
   rightfield->Allocate();

   //mexPrintf("done field->Allocate(); %f sec\n", timer.elapsed());
   //timer.restart();

   VectorPixelType * leftptr = leftfield->GetBufferPointer();
   const VectorPixelType * const leftbuff_end = leftptr + numPix;
   VectorPixelType * rightptr = rightfield->GetBufferPointer();
   
   while ( leftptr != leftbuff_end )
   {
      for (unsigned int d=0; d<Dimension; d++)
      {
         (*leftptr)[d] = *(leftinptrs[d])++;
         (*rightptr)[d] = *(rightinptrs[d])++;
      }
      ++leftptr;
      ++rightptr;
   }

   //mexPrintf("done inputs copy %f sec\n", timer.elapsed());
   //timer.restart();

   
   // Compute distance between jacobians
   typename WarpGradientCalculatorType::Pointer
      leftgradcomp =  WarpGradientCalculatorType::New();
   leftgradcomp->SetInputImage( leftfield );
   
   typename WarpGradientCalculatorType::Pointer
      rightgradcomp =  WarpGradientCalculatorType::New();
   rightgradcomp->SetInputImage( rightfield );

   FieldIteratorType leftIter( leftfield, region );
   leftIter.GoToBegin();
   FieldIteratorType rightIter( rightfield, region );
   rightIter.GoToBegin();

   double fieldGradDist2(0.0), tmp;
   while ( not leftIter.IsAtEnd() )
   {
      tmp = ( ( leftgradcomp->EvaluateAtIndex(leftIter.GetIndex())
                - rightgradcomp->EvaluateAtIndex(rightIter.GetIndex()) )
              .GetVnlMatrix() ).frobenius_norm();
      fieldGradDist2 += tmp*tmp;
      ++leftIter;
      ++rightIter;
   }

   //mexPrintf("done computing output %f sec\n", timer.elapsed());
   //timer.restart();


   plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
   double * outptr = mxGetPr(plhs[0]);
   *outptr = std::sqrt( fieldGradDist2/static_cast<double>(numPix) );

   //mexPrintf("done output copy %f sec\n", timer.elapsed());
}


void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   if (nrhs!=4 and nrhs!=6)
   {
      mexErrMsgTxt("Four or six inputs required.");
   }

   const int dim=nrhs/2;

   const mxClassID classID = mxGetClassID(prhs[0]);

   /* The inputs must be noncomplex double matrices.*/
   for (int n=0; n<nrhs; n++)
   {
      if ( mxGetClassID(prhs[n])!=classID || mxIsComplex(prhs[n]) )
      {
         mexErrMsgTxt("Input must be a noncomplex double.");
      }

      if ( mxGetNumberOfDimensions(prhs[n]) != dim )
      {
         mexErrMsgTxt("The dimension of the inputs must be half the number of inputs.");
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
            deffieldjacobiandist<float,2>(nlhs, plhs, nrhs, prhs);
            break;
         case mxDOUBLE_CLASS:
            deffieldjacobiandist<double,2>(nlhs, plhs, nrhs, prhs);
            break;
         default:
            mexErrMsgTxt("Pixel type unsupported.");
      }
      break;
   case 3:
      switch ( classID )
      {
         case mxSINGLE_CLASS:    
            deffieldjacobiandist<float,3>(nlhs, plhs, nrhs, prhs);
            break;
         case mxDOUBLE_CLASS:
            deffieldjacobiandist<double,3>(nlhs, plhs, nrhs, prhs);
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


