#ifndef __itkPolyaffineLogDomainDeformableRegistrationFilter_txx
#define __itkPolyaffineLogDomainDeformableRegistrationFilter_txx

#include "itkPolyaffineLogDomainDeformableRegistrationFilter.h"

#include "itkExceptionObject.h"
#include "itkImageRegionIterator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkDataObject.h"

#include "itkGaussianOperator.h"
#include "itkVectorNeighborhoodOperatorImageFilter.h"

#include "itkImageDuplicator.h"
#include "itkSubtractImageFilter.h"
#include "itkGradientToMagnitudeImageFilter.h"
#include "itkTransformToVelocityFieldSource.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkMultiplyImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkSquareImageFilter.h"
#include "itkLogImageFilter.h"
#include "itkAddConstantToImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include "vnl/vnl_math.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/algo/vnl_chi_squared.h"
#include "vnl/algo/vnl_cholesky.h"
#include <vnl/vnl_trace.h>

#include <iostream>
#include <fstream>

namespace itk {

// Default constructor
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::PolyaffineLogDomainDeformableRegistrationFilter()
{
 
  this->SetNumberOfRequiredInputs(2);

  //this->SetNumberOfIterations(10);
 
  unsigned int j;
  for( j = 0; j < ImageDimension; j++ )
    {
    m_StandardDeviations[j] = 1.0;
    m_UpdateFieldStandardDeviations[j] = 1.0;
    }

  m_TempField = VelocityFieldType::New();
  m_MaximumError = 0.1;
  m_MaximumKernelWidth = 30;
  m_StopRegistrationFlag = false;

  m_SmoothVelocityField = true;
  m_SmoothUpdateField = false;

  m_Exponentiator = FieldExponentiatorType::New();
  m_Exponentiator->ComputeInverseOff();
  
  m_InverseExponentiator = FieldExponentiatorType::New();
  m_InverseExponentiator->ComputeInverseOn();
  
    m_FixedMaskImage = 0;
    m_MovingMaskImage = 0;
    
    m_ExplainedVelocityField = 0;
		
    //m_AlphaWeightImage = 0;
	
	//m_TotalWeightImage = 0;
  
    m_StartLevel = 0;
    m_EndLevel = 0;
//    m_EstimateLevel = 0;
//    m_CurrentLevel = 0;
//    m_LevelElapsedIterations = 0;
    
    m_VarianceVelocityField = 1;
    m_UnionMaskDilationRadius = 0;
    
//    m_Rank = 0;
    
}


// Set the fixed image.
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::SetFixedImage(
  const FixedImageType * ptr )
{
  this->ProcessObject::SetNthInput( 1, const_cast< FixedImageType * >( ptr ) );
}


// Get the fixed image.
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
const typename PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::FixedImageType *
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GetFixedImage() const
{
  return dynamic_cast< const FixedImageType * >
    ( this->ProcessObject::GetInput( 1 ) );
}


// Set the moving image.
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::SetMovingImage(
  const MovingImageType * ptr )
{
  this->ProcessObject::SetNthInput( 2, const_cast< MovingImageType * >( ptr ) );
}


// Get the moving image.
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
const typename PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::MovingImageType *
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GetMovingImage() const
{
  return dynamic_cast< const MovingImageType * >
    ( this->ProcessObject::GetInput( 2 ) );
}


template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
std::vector<SmartPointer<DataObject> >::size_type
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GetNumberOfValidRequiredInputs() const
{
  typename std::vector<SmartPointer<DataObject> >::size_type num = 0;

  if (this->GetFixedImage())
    {
    num++;
    }

  if (this->GetMovingImage())
    {
    num++;
    }
  
  return num;
}


// Set the standard deviations.
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::SetStandardDeviations(
  double value )
{

  unsigned int j;
  for( j = 0; j < ImageDimension; j++ )
    {
    if( value != m_StandardDeviations[j] )
      {
      break;
      }
    }
  if( j < ImageDimension )
    {
    this->Modified();
    for( j = 0; j < ImageDimension; j++ )
      {
      m_StandardDeviations[j] = value;
      }
    }

}

// Set the standard deviations.
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::SetUpdateFieldStandardDeviations(
  double value )
{

  unsigned int j;
  for( j = 0; j < ImageDimension; j++ )
    {
    if( value != m_UpdateFieldStandardDeviations[j] )
      {
      break;
      }
    }
  if( j < ImageDimension )
    {
    this->Modified();
    for( j = 0; j < ImageDimension; j++ )
      {
      m_UpdateFieldStandardDeviations[j] = value;
      }
    }

}


// Standard PrintSelf method.
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Smooth velocity field: "
     << (m_SmoothVelocityField ? "on" : "off") << std::endl;
  os << indent << "Standard deviations: [";
  unsigned int j;
  for( j = 0; j < ImageDimension - 1; j++ )
    {
    os << m_StandardDeviations[j] << ", ";
    }
  os << m_StandardDeviations[j] << "]" << std::endl;
  os << indent << "Smooth update field: "
     << (m_SmoothUpdateField ? "on" : "off") << std::endl;
  os << indent << "Update field standard deviations: [";
  for( j = 0; j < ImageDimension - 1; j++ )
    {
    os << m_UpdateFieldStandardDeviations[j] << ", ";
    }
  os << m_UpdateFieldStandardDeviations[j] << "]" << std::endl;
  os << indent << "StopRegistrationFlag: ";
  os << m_StopRegistrationFlag << std::endl;
  os << indent << "MaximumError: ";
  os << m_MaximumError << std::endl;
  os << indent << "MaximumKernelWidth: ";
  os << m_MaximumKernelWidth << std::endl;
  os << indent << "Exponentiator: ";
  os << m_Exponentiator << std::endl;
  os << indent << "InverseExponentiator: ";
  os << m_InverseExponentiator << std::endl;

}


// Set the function state values before each iteration
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::InitializeIteration()
{
  //std::cout<<"PolyaffineLogDomainDeformableRegistrationFilter::InitializeIteration"<<std::endl;
  MovingImageConstPointer movingPtr = this->GetMovingImage();
  FixedImageConstPointer fixedPtr = this->GetFixedImage();

  if( !movingPtr || !fixedPtr )
    {
    itkExceptionMacro( << "Fixed and/or moving image not set" );
    }

  // update variables in the equation object
  PDEDeformableRegistrationFunctionType *f = 
    dynamic_cast<PDEDeformableRegistrationFunctionType *>
    (this->GetDifferenceFunction().GetPointer());

  if ( !f )
    {
    itkExceptionMacro(<<"FiniteDifferenceFunction not of type PolyaffineLogDomainDeformableRegistrationFilterFunction");
    }

  f->SetFixedImage( fixedPtr );
  f->SetMovingImage( movingPtr );
  
  this->Superclass::InitializeIteration();

}


/* Override the default implementation for the case when the 
 * initial velocity is not set.
 * If the initial velocity is not set, the output is
 * fill with zero vectors.*/
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::CopyInputToOutput()
{

  typename Superclass::InputImageType::ConstPointer  inputPtr  = this->GetInput();
  
  if( inputPtr )
    {
    this->Superclass::CopyInputToOutput();
    }
  else
    {
    typename Superclass::PixelType zeros;
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      zeros[j] = 0;
      }

    typename OutputImageType::Pointer output = this->GetOutput();
  
    ImageRegionIterator<OutputImageType> out(output, output->GetRequestedRegion());

    while( ! out.IsAtEnd() )
      {
      out.Value() =  zeros;
      ++out;
      }
    }
}


template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GenerateOutputInformation()
{
  //std::cout<<"PolyaffineLogDomainDeformableRegistrationFilter::GenerateOutputInformation"<<std::endl;
  typename DataObject::Pointer output;

  if( this->GetInput(0) )
    {
    // Initial velocity field is set.
    // Copy information from initial field.
    this->Superclass::GenerateOutputInformation();

    }
  else if( this->GetFixedImage() )
    {
    // Initial deforamtion field is not set. 
    // Copy information from the fixed image.
    for (unsigned int idx = 0; idx < 
           this->GetNumberOfOutputs(); ++idx )
      {
      output = this->GetOutput(idx);
      if (output)
        {
        output->CopyInformation(this->GetFixedImage());
        }  
      }

    }

}


template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GenerateInputRequestedRegion()
{
  //std::cout<<"PolyaffineLogDomainDeformableRegistrationFilter::GenerateInputRequestedRegion"<<std::endl;
  // call the superclass's implementation
  Superclass::GenerateInputRequestedRegion();

  // request the largest possible region for the moving image
  MovingImagePointer movingPtr = 
    const_cast< MovingImageType * >( this->GetMovingImage() );
  if( movingPtr )
    {
    movingPtr->SetRequestedRegionToLargestPossibleRegion();
    }
  
  // just propagate up the output requested region for
  // the fixed image and initial velocity field.
  VelocityFieldPointer inputPtr = 
    const_cast< VelocityFieldType * >( this->GetInput() );
  VelocityFieldPointer outputPtr = this->GetOutput();
  FixedImagePointer fixedPtr = 
    const_cast< FixedImageType *>( this->GetFixedImage() );

  if( inputPtr )
    {
    inputPtr->SetRequestedRegion( outputPtr->GetRequestedRegion() );
    }

  if( fixedPtr )
    {
    fixedPtr->SetRequestedRegion( outputPtr->GetRequestedRegion() );
    }
}


// Release memory of internal buffers
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::PostProcessOutput()
{
  this->Superclass::PostProcessOutput();
  m_TempField->Initialize();
}


// Initialize flags
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::Initialize()
{
    
    
	// try union of fixed and mask instead of union of all
	typedef ImageDuplicator< FixedMaskImageType > DuplicatorLabelImageType;
	typename DuplicatorLabelImageType::Pointer imageDuplicator = DuplicatorLabelImageType::New();
	imageDuplicator->SetInputImage(m_FixedMaskImage);
	imageDuplicator->Update();
	typedef itk::ImageRegionIterator< FixedMaskImageType > IteratorLabelType;
	IteratorLabelType iterMaskedUnion(imageDuplicator->GetOutput(), imageDuplicator->GetOutput()->GetRequestedRegion());
	IteratorLabelType iterMaskedMoving(m_MovingMaskImage, m_MovingMaskImage->GetRequestedRegion());
	for(iterMaskedUnion.GoToBegin(), iterMaskedMoving.GoToBegin(); !iterMaskedUnion.IsAtEnd(); ++iterMaskedUnion, ++iterMaskedMoving) {
		if(iterMaskedMoving.Get() || iterMaskedUnion.Get())
			iterMaskedUnion.Set(1);
        else
            iterMaskedUnion.Set(0);
	}
        
//    typedef itk::ImageFileWriter< FixedMaskImageType > MaskWriterType;
//    typename MaskWriterType::Pointer maskWriter = MaskWriterType::New();
//    maskWriter->SetFileName("union_mask_image_before_dilation.mhd");
//    maskWriter->SetInput(imageDuplicator->GetOutput());
//    maskWriter->Update();

    // dilate union image    
    typedef itk::BinaryBallStructuringElement< typename FixedMaskImageType::PixelType,Dim > StructuringElementType;
    StructuringElementType structuringElement;
    structuringElement.SetRadius(m_UnionMaskDilationRadius);
    structuringElement.CreateStructuringElement(); 
    
    typedef itk::BinaryDilateImageFilter< FixedMaskImageType,FixedMaskImageType,StructuringElementType > BinaryDilateImageFilterType;    
    typename BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
    dilateFilter->SetInput(imageDuplicator->GetOutput());
    dilateFilter->SetDilateValue(1);
    dilateFilter->SetKernel(structuringElement);
    dilateFilter->Update();
    
    m_UnionMaskImage = dilateFilter->GetOutput();
    m_UnionMaskImage->DisconnectPipeline();
    
//    maskWriter->SetFileName("union_mask_image_after_dilation.mhd");
//    maskWriter->SetInput(m_UnionMaskImage);
//    maskWriter->Update();
    
    // prepate coordinate image	for Sij	
	m_CoordinatesImage = MatrixImageType::New();
	m_CoordinatesImage->CopyInformation(m_UnionMaskImage);
	m_CoordinatesImage->SetRegions(m_UnionMaskImage->GetRequestedRegion());
	m_CoordinatesImage->Allocate();
	
	itk::ImageRegionIteratorWithIndex<MatrixImageType> coordinatesIt(m_CoordinatesImage, m_CoordinatesImage->GetRequestedRegion());
	for(coordinatesIt.GoToBegin(); !coordinatesIt.IsAtEnd(); ++coordinatesIt) {
		
		vnl_matrix< MatrixElementType > x(DimHom,1);
		for(unsigned int i = 0; i < Dim; ++i) x(i,0) = coordinatesIt.GetIndex()[i] * m_UnionMaskImage->GetSpacing()[i];
		x(Dim,0) = 1.0;
		
		MatrixType matrix = x*x.transpose();
		coordinatesIt.Set(matrix);
	}
    
    std::vector< vnl_vector_fixed< MatrixElementType,Dim > > regionCenters;
    std::vector< vnl_matrix_fixed < MatrixElementType,Dim,Dim > > regionSigmasInverse;
    std::vector< MatrixElementType > regionSigmasDet;
    
    for(unsigned int l = m_StartLevel; l <= m_EndLevel; ++l) {

        std::vector< vnl_matrix_fixed < MatrixElementType,Dim,Dim > > regionSigmas;
        
        // compute total weight image        
        std::vector< PolyaffineTreeTransformType::GaussianTreeNodePointerType > nodeSet;
        m_PolyaffineTree->GetWeightsAtLevel(l, nodeSet);
        for(unsigned int i = 0; i < nodeSet.size(); i++) {
            
            vnl_vector_fixed< MatrixElementType,Dim > regionCenter;
            vnl_matrix_fixed < MatrixElementType,Dim,Dim > regionSigma;
            vnl_matrix< MatrixElementType > regionSigmaInverse;
            
            PolyaffineTreeTransformType::ParametersType param;
            param = nodeSet[i]->Get();
            m_PolyaffineTree->ParamToGaussian(param, regionCenter, regionSigma);
            
            regionCenters.push_back(regionCenter);
            regionSigmas.push_back(regionSigma);
            
            regionSigmaInverse = vnl_matrix_inverse< MatrixElementType >(regionSigma);
            regionSigmasInverse.push_back(regionSigmaInverse);
            
            MatrixElementType det = vnl_determinant(regionSigma);
            if(det>0)
                det = 1.0/std::sqrt(vnl_determinant(regionSigma));
            //MatrixElementType det = 1.0;
            regionSigmasDet.push_back(det);
            
        }
    }
//                
//        typename WeightImageType::IndexType index;
//        typename WeightImageType::PointType point;
//        vnl_matrix< MatrixElementType > x(Dim,1);
//        
//        WeightImagePointer totalWeightImage = WeightImageType::New();
//        totalWeightImage->CopyInformation(m_UnionMaskImage);
//        totalWeightImage->SetRegions(m_UnionMaskImage->GetRequestedRegion());
//        totalWeightImage->Allocate();
//        totalWeightImage->FillBuffer(0.0);
//
//        typedef itk::ImageRegionIteratorWithIndex< WeightImageType > WeightImageIndexIteratorType;
//        WeightImageIndexIteratorType totalWeightIter(totalWeightImage, totalWeightImage->GetRequestedRegion());
//        for(totalWeightIter.GoToBegin(); !totalWeightIter.IsAtEnd(); ++totalWeightIter) {
//        
//            // get physical position
//            index = totalWeightIter.GetIndex();
//            totalWeightImage->TransformIndexToPhysicalPoint(index,point);
//        
//            for(unsigned int i = 0; i < regionCenters.size(); ++i) {
//            
//                for(unsigned int k = 0; k < Dim; ++k)
//                    x(k,0) = point[k]-regionCenters[i][k];
//                
//                // test start
//                //double mahalanobisDistance = (x.transpose()*m_RegionSigmasInverse[i]*x)(0,0);
//                //if(mahalanobisDistance < regionRange) {
//                // test end
//            
//                totalWeightIter.Set( totalWeightIter.Get() + regionSigmasDet[i]*std::exp(-0.5*(x.transpose()*regionSigmasInverse[i]*x)(0,0)) );
//                    
//                // test start
//                //}
//                // test end
//                
//            }
//        }
//        
//        m_TotalWeightImages.push_back(totalWeightImage);
//        
//    }
    
//    // non-null part for rank in pseudo inverse
//    const unsigned int noOfLevelsPrevious = m_EstimateLevel-1;
//    const unsigned int treeDepthPreviousLevel = noOfLevelsPrevious-1;
//    const unsigned int noOfNodesPreviousLevels = std::pow((double)2,(double)(treeDepthPreviousLevel+1))-1;
//    m_Rank = noOfNodesPreviousLevels*Dim*DimHom;
    
    FixedImagePointer fixedPtr = const_cast< FixedImageType *>( this->GetFixedImage() );
    m_PolyaffineTree->SetOutputParametersFromImage(fixedPtr);
    m_TotalWeightImage = m_PolyaffineTree->ComputeTotalWeightImage();    
    ComputeSigma(regionCenters,regionSigmasInverse,regionSigmasDet);
    
	// prepate coordinate image	for Bi
	m_CoorHomImage = CoorHomImageType::New();
	m_CoorHomImage->CopyInformation(m_UnionMaskImage);
	m_CoorHomImage->SetRegions(m_UnionMaskImage->GetRequestedRegion());
	m_CoorHomImage->Allocate();
	itk::ImageRegionIteratorWithIndex<CoorHomImageType> coorHomIt(m_CoorHomImage, m_CoorHomImage->GetRequestedRegion());
	for(coorHomIt.GoToBegin(); !coorHomIt.IsAtEnd(); ++coorHomIt) {
		
		vnl_matrix< MatrixElementType > x(DimHom,1);
		for(unsigned int i = 0; i < Dim; ++i) 
			x(i,0) = coorHomIt.GetIndex()[i] * m_UnionMaskImage->GetSpacing()[i];
		x(Dim,0) = 1.0;
		
		coorHomIt.Set(x);
	}
    	        
  //std::cout<<"PolyaffineLogDomainDeformableRegistrationFilter::Initialize"<<std::endl;
  this->Superclass::Initialize();
  m_StopRegistrationFlag = false;
}


// Smooth velocity using a separable Gaussian kernel
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::SmoothVelocityField()
{
  // The output buffer will be overwritten with new data.
  this->SmoothGivenField(this->GetOutput(), this->m_StandardDeviations);
}

// Smooth update field using a separable Gaussian kernel
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::SmoothUpdateField()
{
  // The update buffer will be overwritten with new data.
  this->SmoothGivenField(this->GetUpdateBuffer(), this->m_UpdateFieldStandardDeviations);
}

// Smooth velocity using a separable Gaussian kernel
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::SmoothGivenField(VelocityFieldType * field, const double StandardDeviations[ImageDimension])
{
    //std::cout << "SmoothGivenField" << std::endl;
    //double startSmoothGivenField = clock();
    //double startBi = 0, endBi = 0, startSolveSystem = 0, endSolveSystem = 0, startCompose= 0, endCompose = 0, startSigma = 0, endSigma = 0;
            
//    typedef itk::ImageDuplicator< VelocityFieldType > DuplicatorType;
//    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
//    duplicator->SetInputImage(field);
//    duplicator->Update();
//    VelocityFieldPointer residualField = duplicator->GetOutput();
    
    VelocityFieldPointer residualField;
    if(m_ExplainedVelocityField) {
        typedef itk::SubtractImageFilter< VelocityFieldType > SubtractImageFilterType;
        typename SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New();
        subtractFilter->SetInput1(field);
        subtractFilter->SetInput2(m_ExplainedVelocityField);
        subtractFilter->Update();
        residualField = subtractFilter->GetOutput();
        residualField->DisconnectPipeline();

    }
    else {
        residualField = field;
    }
    
    /*itk::ImageRegionIterator< VelocityFieldType > iterResidualField(residualField, residualField->GetRequestedRegion());
    float varianceVelocityField = 0;
    for(iterResidualField.GoToBegin(); !iterResidualField.IsAtEnd(); ++iterResidualField) {
        typename VelocityFieldType::PixelType value = iterResidualField.Get();
        for(unsigned int i = 0; i < Dim; ++i)
            varianceVelocityField += value[i]*value[i];
    }
    std::cout << "varianceVelocityField = " << varianceVelocityField << std::endl;*/
    
//    FixedImagePointer fixedPtr = const_cast< FixedImageType *>( this->GetFixedImage() );
//    m_PolyaffineTree->SetOutputParametersFromImage(fixedPtr);
        
    itk::ImageRegionIterator<FixedMaskImageType> iterMaskedUnion(m_UnionMaskImage, m_UnionMaskImage->GetRequestedRegion());
    std::vector< vnl_vector_fixed< MatrixElementType,Dim > > regionCenters;
    std::vector< vnl_matrix_fixed < MatrixElementType,Dim,Dim > > regionSigmasInverse;
    std::vector< MatrixElementType > regionSigmasDet;

    for(unsigned int level = m_StartLevel; level <= m_EndLevel; ++level) {
        
        std::vector< vnl_matrix_fixed < MatrixElementType,Dim,Dim > > regionSigmas;
        
        // compute total weight image        
        std::vector< PolyaffineTreeTransformType::GaussianTreeNodePointerType > weightNodeSet;
        m_PolyaffineTree->GetWeightsAtLevel(level, weightNodeSet);
        for(unsigned int i = 0; i < weightNodeSet.size(); i++) {
            
            vnl_vector_fixed< MatrixElementType,Dim > regionCenter;
            vnl_matrix_fixed < MatrixElementType,Dim,Dim > regionSigma;
            vnl_matrix< MatrixElementType > regionSigmaInverse;
            
            PolyaffineTreeTransformType::ParametersType param;
            param = weightNodeSet[i]->Get();
            m_PolyaffineTree->ParamToGaussian(param, regionCenter, regionSigma);
            
            regionCenters.push_back(regionCenter);
            regionSigmas.push_back(regionSigma);
            
            regionSigmaInverse = vnl_matrix_inverse< MatrixElementType >(regionSigma);
            regionSigmasInverse.push_back(regionSigmaInverse);
            
            MatrixElementType det = vnl_determinant(regionSigma);
            if(det>0)
                det = 1.0/std::sqrt(vnl_determinant(regionSigma));
            //MatrixElementType det = 1.0;
            regionSigmasDet.push_back(det);
            
        }
    }
                    
        // compute B
        //startBi = clock();
        //std::cout << "compute B" << std::endl;
        VelocityHomImagePointer vectorHomField = VelocityHomImageType::New();
        vectorHomField->CopyInformation(residualField);
        vectorHomField->SetRegions(residualField->GetRequestedRegion());
        vectorHomField->Allocate();
        typedef itk::ImageRegionConstIterator< VelocityFieldType > FieldIteratorType;
        typedef itk::ImageRegionIterator< VelocityHomImageType > HomFieldIteratorType;
        FieldIteratorType fieldIter(residualField, residualField->GetRequestedRegion());
        HomFieldIteratorType homFieldIter(vectorHomField, vectorHomField->GetRequestedRegion());	
        for(fieldIter.GoToBegin(), homFieldIter.GoToBegin(); !fieldIter.IsAtEnd(); ++fieldIter, ++homFieldIter) {
            
            vnl_matrix< MatrixElementType > velVnl(Dim,1);		
            typename VelocityFieldType::PixelType vel = fieldIter.Get();
            for(unsigned int i = 0; i < Dim; ++i)
                velVnl(i,0) = vel[i];
            homFieldIter.Set(velVnl);
        }	
        
        vnl_matrix< MatrixElementType > B(Dim*DimHom*regionCenters.size(),1);
        vnl_matrix< MatrixElementType > Bi(Dim*DimHom,1);
        vnl_matrix< MatrixElementType > Id(Dim,Dim);
        Id.set_identity();
        
        for(unsigned int i = 0; i < B.rows(); i += Bi.rows()) {
            
            unsigned int wi = (i-(i%Bi.rows()))/Bi.rows();
            //unsigned int currentLevel = ceil(log2((wi+1)+1)-1);
            //std::cout << "region i = " << wi << " at level = " << currentLevel << std::endl;
            
            typename WeightImageType::IndexType index;
            typename WeightImageType::PointType point;
            
            MatrixElementType weight;
            vnl_matrix< MatrixElementType > x(Dim,1);

            ImageRegionConstIterator<VelocityHomImageType> inputIt(vectorHomField, vectorHomField->GetRequestedRegion());
            //ImageRegionConstIterator<WeightImageType> totalWeightIt(m_TotalWeightImages[currentLevel], m_TotalWeightImages[currentLevel]->GetRequestedRegion());
            ImageRegionConstIterator<WeightImageType> totalWeightIt(m_TotalWeightImage, m_TotalWeightImage->GetRequestedRegion());
            ImageRegionConstIteratorWithIndex<FixedMaskImageType> maskIt(m_UnionMaskImage, m_UnionMaskImage->GetRequestedRegion());
            ImageRegionConstIterator<CoorHomImageType> coorHomIt(m_CoorHomImage, m_CoorHomImage->GetRequestedRegion());
            Bi.fill(0.0);

            for(inputIt.GoToBegin(), totalWeightIt.GoToBegin(), maskIt.GoToBegin(), coorHomIt.GoToBegin(); 
                    !inputIt.IsAtEnd(); ++inputIt, ++totalWeightIt, ++maskIt, ++coorHomIt) {
                
                if(maskIt.Get()) {
                    
                    // get physical position
                    index = maskIt.GetIndex();
                    //m_TotalWeightImages[currentLevel]->TransformIndexToPhysicalPoint(index,point);
                    m_TotalWeightImage->TransformIndexToPhysicalPoint(index,point);
                    
                    // weight
                    for(unsigned int k = 0; k < Dim; ++k)
                        x(k,0) = point[k]-regionCenters[wi][k];
                    
                    weight = regionSigmasDet[wi]*std::exp(-0.5*(x.transpose()*regionSigmasInverse[wi]*x)(0,0)) / totalWeightIt.Get();
                    
                    //Bi += inputIt.Get().GetVnlMatrix()*coorHomIt.Get().GetVnlMatrix()*weight;
                    // Kronecker product
                    Bi += KroneckerProduct(coorHomIt.Get().GetVnlMatrix(),Id)*inputIt.Get().GetVnlMatrix()*weight;
                    
                }
                
            }
                                    
            // fill big matrix
            for(unsigned int m = 0; m < Bi.rows(); ++m)
                B(i+m,0) = Bi(m,0);
            
        }
    
        //endBi = clock();
        //std::cout << "B = \n" << B << std::endl;
        
        // joint estimate of all log matrices
        //startSolveSystem = clock();
        //std::cout << "joint estimate" << std::endl;
        
        //std::cout << "m_Mui = \n" << m_Mui << std::endl;
        vnl_matrix< MatrixElementType > M;
        if(!m_VarianceVelocityField) { // no prior
            M = m_SigmaCheck * B;
        }
        else {
            M = m_SigmaCheck * ((1.0/m_VarianceVelocityField)*B + m_SigmaInverse*m_Mui);
        }
    
//        // compute posterior prior
//        vnl_cholesky sigmaNewCholesky(m_SigmaCheck+M*M.transpose());
//        m_SigmaNewInverse = sigmaNewCholesky.inverse();
//    
//        // optimize \Lambda
//        vnl_cholesky cCholesky(m_CZero);
//        vnl_matrix< MatrixElementType > Id_N(regionCenters.size(),regionCenters.size());
//        Id_N.set_identity();
//        vnl_matrix< MatrixElementType > A = m_SigmaNewInverse * KroneckerProduct(Id_N,cCholesky.inverse());
//        m_LambdaNew.set_size(regionCenters.size(),regionCenters.size());
//        for(unsigned int i = 0; i < regionCenters.size(); ++i) {
//            for(unsigned int j = 0; j < regionCenters.size(); ++j) {
//                m_LambdaNew(i,j) = (1.0/Dim*DimHom) * vnl_trace(A.extract(Dim*DimHom, Dim*DimHom, i*Dim*DimHom, j*Dim*DimHom));
//            }
//        }
//    
//        // optimize C
//        //vnl_cholesky lambdaCholesky(m_LambdaZero);
//        vnl_cholesky lambdaCholesky(m_LambdaNew);
//        vnl_matrix< MatrixElementType > Id_12(Dim*DimHom,Dim*DimHom);
//        Id_12.set_identity();
//        vnl_matrix< MatrixElementType > D = KroneckerProduct(lambdaCholesky.inverse(),Id_12) * m_SigmaNewInverse;
//        m_CNew.set_size(Dim*DimHom,Dim*DimHom);
//        m_CNew.fill(0.0);
//        for(unsigned int i = 0; i < regionCenters.size(); ++i) {
//            m_CNew += D.extract(Dim*DimHom,Dim*DimHom,i*Dim*DimHom,i*Dim*DimHom);
//        }
//        m_CNew /= regionCenters.size();

        //endSolveSystem = clock();
        //std::cout << "M = \n" << M << std::endl;
        
        //startCompose = clock();

        // set nodes
        unsigned int noOfLogs = M.rows()/(Dim*DimHom);
        //std::cout << "noOfLogs = " << noOfLogs << std::endl;
        std::vector< PolyaffineTreeTransformType::NodeValueType > Ms(noOfLogs);

        for(unsigned int k = 0; k < noOfLogs; ++k) {
            
            vnl_matrix< MatrixElementType > Munvectorized(Dim,DimHom);
            for(unsigned int n = 0; n < DimHom; ++n) {
                for(unsigned int m = 0; m < Dim; ++m) {
                    Munvectorized(m,n) = M(k*Dim*DimHom+n*Dim+m,0);
                }
            }
            //std::cout << "Munvectorized = \n" << Munvectorized << std::endl;
            Ms[k].set_size(Dim*DimHom);
            for(unsigned int m = 0; m < Dim; ++m) {
                for(unsigned int n = 0; n < DimHom; ++n) {
                    Ms[k](m*DimHom+n) = Munvectorized(m,n);
                }
            }
            //std::cout << "Ms[" << k << "] = \n" << Munvectorized << std::endl;
            
        }
        
    unsigned int regionCounter = 0;
    for(unsigned int level = m_StartLevel; level <= m_EndLevel; ++level) {
        
        std::vector< PolyaffineTreeTransformType::GaussianTreeNodePointerType > nodeSetWeights;
        m_PolyaffineTree->GetWeightsAtLevel(level, nodeSetWeights);
        //unsigned int noOfRegions = std::pow((double)2.0,(double)level);
        unsigned int noOfRegions = nodeSetWeights.size();

        std::vector< PolyaffineTreeTransformType::GaussianTreeNodePointerType > transformNodeSet(noOfRegions);
        for(unsigned int k = 0; k < transformNodeSet.size(); ++k) {
            PolyaffineTreeTransformType::GaussianTreeNodePointerType newNode = PolyaffineTreeTransformType::GaussianTreeNodeType::New();
            newNode->Set(Ms[regionCounter]);
            transformNodeSet[k] = newNode;
            ++regionCounter;
        }
        m_PolyaffineTree->SetTransformationsAtLevel(level, transformNodeSet);
        
    }
    
    //endCompose = clock();
    VelocityFieldPointer estimatedField = m_PolyaffineTree->GetParametersAsVectorField(m_EndLevel+1);
    
    // field to contain the final smoothed data, do the equivalent of a graft
    field->SetPixelContainer( estimatedField->GetPixelContainer() );
    field->SetRequestedRegion( estimatedField->GetRequestedRegion() );
    field->SetBufferedRegion( estimatedField->GetBufferedRegion() );
    field->SetLargestPossibleRegion( estimatedField->GetLargestPossibleRegion() );
    field->CopyInformation( estimatedField );
        
    //std::cout << "smooth given field: " << (clock()-startSmoothGivenField)/CLOCKS_PER_SEC << std::endl;


}

//	template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
//	vnl_matrix< typename PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>::MatrixElementType > 
//	PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>::ComputeCorrelationsBetweenRegions() {
//		
//		vnl_matrix< MatrixElementType > Sdiag(DimHom*m_RegionCenters.size(),DimHom*m_RegionCenters.size());
//		Sdiag.fill(0.0);
//		vnl_matrix< MatrixElementType > Sij(DimHom,DimHom);
//		
//		for(unsigned int i = 0; i < m_S.rows(); i += Sij.rows()) {
//			for(unsigned int j = 0; j < m_S.cols(); j += Sij.cols()) {
//				
//				unsigned int wi = (i-(i%Sij.rows()))/Sij.rows();
//				unsigned int wj = (j-(j%Sij.cols()))/Sij.cols();
//
//				if(wi==wj) {
//					
//					for(unsigned int m = 0; m < Sij.rows(); ++m)
//						for(unsigned int n = 0; n < Sij.cols(); ++n)
//								Sdiag(i+m,j+n) = m_S(i+m,j+n);
//
//				}
//			}
//		}
//				
//		vnl_symmetric_eigensystem< MatrixElementType > eigensystem(Sdiag);
//		vnl_matrix< MatrixElementType > Sinvsqrt = eigensystem.inverse_square_root();
//		vnl_matrix< MatrixElementType > CorrMatrix = Sinvsqrt * m_S * Sinvsqrt;
//		
//		vnl_matrix< MatrixElementType > pMatrix(m_RegionCenters.size(),m_RegionCenters.size());
//		for(unsigned int i = 0; i < m_S.rows(); i += Sij.rows()) {
//			for(unsigned int j = 0; j < m_S.cols(); j += Sij.cols()) {
//				
//				unsigned int wi = (i-(i%Sij.rows()))/Sij.rows();
//				unsigned int wj = (j-(j%Sij.cols()))/Sij.cols();
//				
//				for(unsigned int m = 0; m < Sij.rows(); ++m)
//					for(unsigned int n = 0; n < Sij.cols(); ++n)
//						Sij(m,n) = CorrMatrix(i+m,j+n);
//				
//				/*
//				// Bartlett-Lawley test [Fujikoshi, Biometrika 66(2), 1979]
//				vnl_svd< MatrixElementType > svd(Sij);
//				double sum = 0;
//				for(unsigned int k = 0; k < DimHom; ++k) {
//					// special case
//					std::cout << svd.W(k) << std::endl;
//					if(svd.W(k) > 0.999)
//						sum += NumericTraits<double>::min();
//					else
//						sum += std::log(1-svd.W(k)*svd.W(k));
//				}
//				std::cout << sum << std::endl;
//				double chisq = -(1-DimHom+0.5) * sum;
//				std::cout << chisq << std::endl;
//				double pValue = 1-vnl_chi_squared_cumulative(chisq, DimHom*DimHom);
//				std::cout << pValue << std::endl;
//				
//				pMatrix(wi,wj) = pValue;
//				*/
//				
//				// maximum coorelation
//				vnl_svd< MatrixElementType > svd(Sij);
//				pMatrix(wi,wj) = svd.W(0);
//				
//			}
//		}
//		
//		return pMatrix;
//	}

    
	template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
	void
	PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>::ComputeSigma(std::vector< vnl_vector_fixed< MatrixElementType,Dim > > regionCenters, std::vector< vnl_matrix_fixed < MatrixElementType,Dim,Dim > > regionSigmasInverse, std::vector< MatrixElementType > regionSigmasDet) {
        
		//startSigma = clock();
		//std::cout << "compute Sigma" << std::endl;
        vnl_matrix< MatrixElementType > S(Dim*DimHom*regionCenters.size(),Dim*DimHom*regionCenters.size());
		vnl_matrix< MatrixElementType > Sij(Dim*DimHom,Dim*DimHom);
        vnl_matrix< MatrixElementType > Id(Dim,Dim);
        Id.set_identity();
        
        typename WeightImageType::IndexType index;
        typename WeightImageType::PointType point;
        
        typename WeightImageType::PixelType wiValue;
        typename WeightImageType::PixelType wjValue;
        
        vnl_matrix< MatrixElementType > xi(Dim,1);
        vnl_matrix< MatrixElementType > xj(Dim,1);
        
        vnl_matrix< MatrixElementType > GramMatrix(regionCenters.size(),regionCenters.size());
        GramMatrix.fill(0.0);
                
		for(unsigned int i = 0; i < S.rows(); i += Sij.rows()) {
			for(unsigned int j = 0; j < S.cols(); j += Sij.cols()) {
				
				unsigned int wi = (i-(i%Sij.rows()))/Sij.rows();
//                unsigned int currentLeveli = ceil(log2((wi+1)+1)-1);
				unsigned int wj = (j-(j%Sij.cols()))/Sij.cols();
//                unsigned int currentLevelj = ceil(log2((wj+1)+1)-1);
				
//                std::cout << "compute S: region=" << wi << "; region=" << wj << std::endl;
				if(wi <= wj) {
					
					ImageRegionConstIterator<MatrixImageType> inputIt(m_CoordinatesImage, m_CoordinatesImage->GetRequestedRegion());
//					ImageRegionConstIterator<WeightImageType> totalWeightItLeveli(m_TotalWeightImages[currentLeveli], m_TotalWeightImages[currentLeveli]->GetRequestedRegion());
//                    ImageRegionConstIterator<WeightImageType> totalWeightItLevelj(m_TotalWeightImages[currentLevelj], m_TotalWeightImages[currentLevelj]->GetRequestedRegion());
                    ImageRegionConstIterator<WeightImageType> totalWeightItLevel(m_TotalWeightImage, m_TotalWeightImage->GetRequestedRegion());
                    
					ImageRegionConstIteratorWithIndex<FixedMaskImageType> maskIt(m_UnionMaskImage, m_UnionMaskImage->GetRequestedRegion());
					MatrixType sum;
					sum.Fill(0.0);
//					for(inputIt.GoToBegin(), totalWeightItLeveli.GoToBegin(), totalWeightItLevelj.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++inputIt, ++totalWeightItLeveli, ++totalWeightItLevelj, ++maskIt) {
                    for(inputIt.GoToBegin(), totalWeightItLevel.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++inputIt, ++totalWeightItLevel, ++maskIt) {

                        if(maskIt.Get()) {

                            // get physical position
                            index = maskIt.GetIndex();
                            m_UnionMaskImage->TransformIndexToPhysicalPoint(index,point);
                            
                            // weight
                            for(unsigned int k = 0; k < Dim; ++k)
                                xi(k,0) = point[k]-regionCenters[wi][k];
                            
                            // test start
                            //if((xi.transpose()*m_RegionSigmasInverse[wi]*xi)(0,0) < 3.84) {
                            // test end
                            
                            for(unsigned int k = 0; k < Dim; ++k)
                                xj(k,0) = point[k]-regionCenters[wj][k];

                            // test start
                            //if((xj.transpose()*m_RegionSigmasInverse[wj]*xj)(0,0) < 3.84) {
                            // test end
                            
//                            wiValue = regionSigmasDet[wi]*std::exp(-0.5*(xi.transpose()*regionSigmasInverse[wi]*xi)(0,0)) / totalWeightItLeveli.Get();
//                            wjValue = regionSigmasDet[wj]*std::exp(-0.5*(xj.transpose()*regionSigmasInverse[wj]*xj)(0,0)) / totalWeightItLevelj.Get();                            
                            wiValue = regionSigmasDet[wi]*std::exp(-0.5*(xi.transpose()*regionSigmasInverse[wi]*xi)(0,0)) / totalWeightItLevel.Get();
                            wjValue = regionSigmasDet[wj]*std::exp(-0.5*(xj.transpose()*regionSigmasInverse[wj]*xj)(0,0)) / totalWeightItLevel.Get();                            
                        
                            GramMatrix(wi,wj) += wiValue*wjValue;
							sum += inputIt.Get()*wiValue*wjValue;
                                
                            // test start
                            //}
                            //}
                            // test end
                            
                        }
                        
					}
					//Sij = sum.GetVnlMatrix();
                    // Kronecker product
                    Sij = KroneckerProduct(sum.GetVnlMatrix(), Id);
                    
				}
				else
					Sij.fill(0.0);
                				
				// fill big matrix				
				for(unsigned int m = 0; m < Sij.rows(); ++m) {
					for(unsigned int n = 0; n < Sij.cols(); ++n) {
						if(wi > wj)
							S(i+m,j+n) = S(j+m,i+n);
						else
							S(i+m,j+n) = Sij(m,n);
					}
				}
				
			}
		}
		//endSigma = clock();
        
        vnl_matrix< MatrixElementType > GramMatrixNoDiag = GramMatrix;
        GramMatrixNoDiag.fill_diagonal(0.0);
        GramMatrix = GramMatrix + GramMatrixNoDiag.transpose();
        //std::cout << "GramMatrix = \n" << GramMatrix << std::endl;
                        
        vnl_matrix< MatrixElementType > PhiTransposePhi = S;
        
        // (9) first line TMI: \check{\Sigma}	& = & ( \Sigma_{0}^{-1} + \frac{1}{\sigma_{v}^{2}} \Phi^{T} \Phi )^{-1}
//        if(m_LambdaZero.empty()) {
//            m_LambdaZero.set_size(regionCenters.size(), regionCenters.size());
//            m_LambdaZero.set_identity();
//            m_CZero.set_size(Dim*DimHom, Dim*DimHom);
//            m_CZero.set_identity();
//        }
        if(m_SigmaInverse.empty()) {
            m_SigmaInverse.set_size(regionCenters.size()*Dim*DimHom, regionCenters.size()*Dim*DimHom);
            m_SigmaInverse.fill(0.0);
            m_Mui.set_size(regionCenters.size()*Dim*DimHom,1);
            m_Mui.fill(0.0);
        }
        
//        vnl_matrix< MatrixElementType > SigmaZero = KroneckerProduct(m_LambdaZero, m_CZero);
        // variance of the noise that we assume in the probabilistic model
//        vnl_cholesky checkSigmaInverse(SigmaZero + (1.0/m_VarianceVelocityField)*PhiTransposePhi);
        
        std::ofstream filePhiTransposePhi;
        filePhiTransposePhi.open ("PhiTransposePhi.txt");
        filePhiTransposePhi << PhiTransposePhi << std::endl;
        filePhiTransposePhi.close();
        
        vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > GammaCheck;
        if(!m_VarianceVelocityField) { // no prior
            GammaCheck = PhiTransposePhi;
//            vnl_matrix_inverse< PolyaffineTreeTransformType::MatrixElementType > checkSigmaInverse(PhiTransposePhi);
//            checkSigmaInverse.zero_out_relative(1.490116e-08);
//            m_SigmaCheck = checkSigmaInverse.pinverse(270);
        }
        else {
            GammaCheck = m_SigmaInverse + (1.0/m_VarianceVelocityField)*PhiTransposePhi;
//            vnl_matrix_inverse< PolyaffineTreeTransformType::MatrixElementType > checkSigmaInverse(m_SigmaInverse + (1.0/m_VarianceVelocityField)*PhiTransposePhi);
//            checkSigmaInverse.zero_out_relative(1.490116e-08);
//            m_SigmaCheck = checkSigmaInverse.pinverse(270);
        }
        
        vnl_svd< PolyaffineTreeTransformType::MatrixElementType > svd(GammaCheck);
        vnl_diag_matrix< PolyaffineTreeTransformType::MatrixElementType > diagD = svd.W();
        for(unsigned int i = 0; i < diagD.diagonal().size(); ++i) {
            if(diagD(i,i) < svd.W()[1]*1.490116e-08)
                diagD(i,i) = 0.0;
            else
                diagD(i,i) = 1.0/diagD(i,i);
        }
        m_SigmaCheck = svd.V() * diagD * svd.U().transpose();
        
        std::ofstream filePhiTransposePhiInv;
        filePhiTransposePhiInv.open ("PhiTransposePhiInv.txt");
        filePhiTransposePhiInv << m_SigmaCheck << std::endl;
        filePhiTransposePhiInv.close();
        
//        vnl_cholesky checkSigmaInverse(m_SigmaInverse + (1.0/m_VarianceVelocityField)*PhiTransposePhi);
//        m_SigmaCheck = checkSigmaInverse.inverse();
                
		//std::cout << "m_SigmaCheck = \n" << m_SigmaCheck << std::endl;
        
	}
    
    
    template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
    vnl_matrix< typename PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
    ::MatrixElementType >
	PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>::KroneckerProduct(vnl_matrix< MatrixElementType > A, vnl_matrix< MatrixElementType > B) {
        
        // C = A (KroneckerProduct) B
        vnl_matrix< MatrixElementType > C(A.rows()*B.rows(),A.cols()*B.cols());
        
        for(unsigned int i = 0; i < A.rows(); ++i) {
            for(unsigned int j = 0; j < A.cols(); ++j) {
                
                for(unsigned int m = 0; m < B.rows(); ++m) {
                    for(unsigned int n = 0; n < B.cols(); ++n) {
                        C(i*B.rows()+m,j*B.cols()+n) = A(i,j)*B(m,n);
                    }
                }

            }
        }
        
        return C;

    }


template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
typename PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::DeformationFieldPointer
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GetDeformationField()
{
  //std::cout<<"LogDomainDeformableRegistration::GetDeformationField"<<std::endl;
  m_Exponentiator->SetInput( this->GetVelocityField() );
  m_Exponentiator->GetOutput()->SetRequestedRegion( this->GetVelocityField()->GetRequestedRegion() );
  m_Exponentiator->Update();
  return m_Exponentiator->GetOutput();
}


template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
typename PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::DeformationFieldPointer
PolyaffineLogDomainDeformableRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GetInverseDeformationField()
{

    //std::cout << "# of children: " << m_PolyaffineTree->GetTransformations() << std::endl;
    if(m_PolyaffineTree->GetTransformations().IsNotNull()) {
        PolyaffineTreeTransformType::Pointer invTransform = PolyaffineTreeTransformType::New();
        m_PolyaffineTree->GetInverse(invTransform);
        FixedImagePointer fixedPtr = const_cast< FixedImageType *>( this->GetFixedImage() );
        invTransform->SetOutputParametersFromImage(fixedPtr);
        PolyaffineTreeTransformType::VectorFieldPointerType vectorField = invTransform->GetDisplacementFieldAsVectorField(m_EndLevel+1);
        return vectorField;
    }
    else {
        m_InverseExponentiator->SetInput( this->GetVelocityField() );
        m_InverseExponentiator->GetOutput()->SetRequestedRegion( this->GetVelocityField()->GetRequestedRegion() );
        m_InverseExponentiator->Update();
        return m_InverseExponentiator->GetOutput();
    }

}


} // end namespace itk

#endif
