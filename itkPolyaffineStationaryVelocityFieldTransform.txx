#ifndef _itkPolyaffineStationaryVelocityFieldTransform_txx_
#define _itkPolyaffineStationaryVelocityFieldTransform_txx_

#include "itkPolyaffineStationaryVelocityFieldTransform.h"

#include <itkNeighborhoodAlgorithm.h>
#include <itkImageRegionIterator.h>
#include <vnl/vnl_det.h>
#include <itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h>
#include <itkMultiplyByConstantImageFilter.h>
#include <itkWarpVectorImageFilter.h>
#include <itkImageRegionIterator.h>  
#include <itkDivideByConstantImageFilter.h>
#include <itkTransformToVelocityFieldSource.h>
#include <itkImageFileWriter.h>

namespace itk
{



template <class TScalarType, class TReferenceType, unsigned int NDimensions>
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
//PolyaffineStationaryVelocityFieldTransform() : Superclass( SpaceDimension, ParametersDimension )
PolyaffineStationaryVelocityFieldTransform() : Superclass( ParametersDimension )
{
    this->m_InterpolateFunction = InterpolateFunctionType::New();
    
    m_TransformationsRootNode = 0;
    m_WeightsRootNode = 0;
    m_TotalWeightImage = 0;
}



template <class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
SetIdentity(void)
{
    TScalarType            value = 0;
    VectorFieldPointerType field = VectorFieldType::New();
    field->Allocate();
    field->FillBuffer(value);
    field->Register();
    this->m_VectorField = field;
}



template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
PrintSelf(std::ostream &os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);
    
    os << indent << "weights = " << this->GetFixedParameters() << std::endl;
    os << indent << "transformations = " << this->GetParameters() << std::endl;
}
    
    
    template<class TScalarType, class TReferenceType, unsigned int NDimensions>
    typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::ScalarImagePointerType
    PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
    ComputeTotalWeightImage()  
    {
        if(m_TotalWeightImage)
            return m_TotalWeightImage;
        
        m_TotalWeightImage = ScalarImageType::New();
        m_TotalWeightImage->CopyInformation(m_ReferenceImage);
        m_TotalWeightImage->SetRegions(m_ReferenceImage->GetRequestedRegion());
        m_TotalWeightImage->Allocate();
        m_TotalWeightImage->FillBuffer(0.0);
        
        std::vector< vnl_vector_fixed< MatrixElementType,NDimensions > > regionCenters;
        std::vector< vnl_matrix_fixed < MatrixElementType,NDimensions,NDimensions > > regionSigmasInverse;
        std::vector< MatrixElementType > regionSigmasDet;
        
        unsigned int noOfLevels = this->GetNumberOfLevels();
        
        for(unsigned int l = 0; l < noOfLevels; ++l) {
            
            std::vector< vnl_matrix_fixed < MatrixElementType,NDimensions,NDimensions > > regionSigmas;
            
            // compute total weight image        
            std::vector< GaussianTreeNodePointerType > nodeSet;
            this->GetWeightsAtLevel(l, nodeSet);
            for(unsigned int i = 0; i < nodeSet.size(); i++) {
                
                vnl_vector_fixed< MatrixElementType,NDimensions > regionCenter;
                vnl_matrix_fixed < MatrixElementType,NDimensions,NDimensions > regionSigma;
                vnl_matrix< MatrixElementType > regionSigmaInverse;
                
                ParametersType param;
                param = nodeSet[i]->Get();
                this->ParamToGaussian(param, regionCenter, regionSigma);
                
                regionCenters.push_back(regionCenter);
                regionSigmas.push_back(regionSigma);
                
                regionSigmaInverse = vnl_matrix_inverse< MatrixElementType >(regionSigma);
                //std::cout << "regionSigmaInverse = " << regionSigmaInverse << std::endl;
                regionSigmasInverse.push_back(regionSigmaInverse);
                
                MatrixElementType det = vnl_determinant(regionSigma);
                if(det>0)
                    det = 1.0/std::sqrt(vnl_determinant(regionSigma));
                //MatrixElementType det = 1.0;
                //std::cout << "det = " << det << std::endl;
                regionSigmasDet.push_back(det);
                
            }
            
            typename ScalarImageType::IndexType index;
            typename ScalarImageType::PointType point;
            vnl_matrix< MatrixElementType > x(NDimensions,1);
                        
            typedef itk::ImageRegionIteratorWithIndex< ScalarImageType > WeightImageIndexIteratorType;
            WeightImageIndexIteratorType totalWeightIter(m_TotalWeightImage, m_TotalWeightImage->GetRequestedRegion());
            for(totalWeightIter.GoToBegin(); !totalWeightIter.IsAtEnd(); ++totalWeightIter) {
                
                // get physical position
                index = totalWeightIter.GetIndex();
                m_TotalWeightImage->TransformIndexToPhysicalPoint(index,point);
                
                for(unsigned int i = 0; i < regionCenters.size(); ++i) {
                    
                    for(unsigned int k = 0; k < NDimensions; ++k)
                        x(k,0) = point[k]-regionCenters[i][k];
                    
                    // test start
                    //double mahalanobisDistance = (x.transpose()*m_RegionSigmasInverse[i]*x)(0,0);
                    //if(mahalanobisDistance < regionRange) {
                    // test end
                    
                    totalWeightIter.Set( totalWeightIter.Get() + regionSigmasDet[i]*std::exp(-0.5*(x.transpose()*regionSigmasInverse[i]*x)(0,0)) );
                    
                    // test start
                    //}
                    // test end
                    
                }
            }
        }
        
        return m_TotalWeightImage;
        
    }

    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::VectorFieldType *
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetParametersAsVectorField(const unsigned int noOfLevels) 
{
    VectorFieldPointerType estimatedField = VectorFieldType::New();
    estimatedField->CopyInformation(m_ReferenceImage);
    estimatedField->SetRegions(m_ReferenceImage->GetRequestedRegion());
    estimatedField->Allocate();
    estimatedField->FillBuffer(0.0);

//    unsigned int counter = 0;
//    this->CountNumberOfNodes(m_WeightsRootNode, counter);
//    unsigned int noOfLevels = log2(counter+1);
    for(unsigned int l = 0; l < noOfLevels; ++l) {
        
        typedef itk::AddImageFilter< VectorFieldType > AddImageFilterType;
        typename AddImageFilterType::Pointer addImage = AddImageFilterType::New();
        addImage->SetInput1(estimatedField);
        addImage->SetInput2(this->GetParametersAsVectorFieldAtLevel(l));
        addImage->Update();
        estimatedField = addImage->GetOutput();
        estimatedField->DisconnectPipeline();
        
    }
    
    this->SetParametersAsVectorField(estimatedField);
    return this->m_VectorField;

}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::VectorFieldType *
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetParametersAsVectorFieldAtLevel(const unsigned int level) 
{
    if (this->m_ReferenceImage.IsNull())
        itkExceptionMacro("No reference image to retrieve output parmeters has been set.");
    
    ScalarImagePointerType totalWeightImage = ComputeTotalWeightImage();
    // debug
    /*typedef itk::ImageFileWriter< ScalarImageType > WeightWriterType;
    typename WeightWriterType::Pointer weightWriter = WeightWriterType::New();
    weightWriter->SetFileName("TotalWeightImage.mhd");
    weightWriter->SetInput(totalWeightImage);
    weightWriter->Update();*/
    
    // create velocity field representation
    VectorFieldPointerType estimatedField = VectorFieldType::New();
    estimatedField->CopyInformation(m_ReferenceImage);
    estimatedField->SetRegions(m_ReferenceImage->GetRequestedRegion());
    estimatedField->Allocate();
    
    // no generation of intermediate weight or velocity images
    typedef itk::ImageRegionIteratorWithIndex< VectorFieldType > FieldIndexIteratorType;
    FieldIndexIteratorType fieldIter(estimatedField, estimatedField->GetRequestedRegion());
    
    std::vector< GaussianTreeNodePointerType > weightSet;
    std::vector< GaussianTreeNodePointerType > transformSet;
    
    this->GetWeightsAtLevel(level, weightSet);
    this->GetTransformationsAtLevel(level, transformSet);

    // some variables defined outside the for loop to save time
    typename VectorFieldType::IndexType index;
    typename VectorFieldType::PointType point;
    typename VectorFieldType::PixelType value;
    vnl_matrix< MatrixElementType > x(NDimensions,1);
    vnl_matrix< MatrixElementType > xHom(NDimensionsHom,1);
    vnl_vector< MatrixElementType > weightAll(weightSet.size());
    vnl_matrix< MatrixElementType > weightedVelocity(NDimensionsHom,1);
    
    // prepare weights for all regions
    std::vector< vnl_vector_fixed< MatrixElementType,NDimensions > > meanAll(weightSet.size());
    std::vector< vnl_matrix< MatrixElementType > > sigmaInverseAll(weightSet.size());
    std::vector< MatrixElementType > detSigma(weightSet.size());
    for(unsigned int i = 0; i < weightSet.size(); ++i) {
        vnl_vector_fixed< MatrixElementType,NDimensions > mean;
        vnl_matrix_fixed< MatrixElementType,NDimensions,NDimensions > covariance;
        ParametersType param;
        
        param = weightSet[i]->Get();
        this->ParamToGaussian(param, mean, covariance);
        
        meanAll[i] = mean;
        sigmaInverseAll[i] = vnl_matrix_inverse< MatrixElementType >(covariance);
        
        detSigma[i] = vnl_determinant(covariance);
        if(detSigma[i]>0)
            detSigma[i] = 1.0/std::sqrt(vnl_determinant(covariance));
        //detSigma[i] = 1.0;
    }
    
    // prepare transfromations for all regions
    std::vector< vnl_matrix_fixed< MatrixElementType,NDimensionsHom,NDimensionsHom > > Mall(weightSet.size());
    for(unsigned int i = 0; i < transformSet.size(); ++i) {
        
        vnl_matrix_fixed< MatrixElementType,NDimensionsHom,NDimensionsHom > Ms;
        ParametersType param;
        param = transformSet[i]->Get();
        this->ParamToLogTransform(param, Ms);
        /*Ms.fill(0.0);
        for(unsigned int m = 0; m < NDimensions; ++m) {
            for(unsigned int n = 0; n < NDimensionsHom; ++n) {
                Ms(m,n) = transformSet[i]->Get()(m*NDimensionsHom+n);
            }
        }*/
        Mall[i] = Ms;
        
    }
    
    // debug
    /*ScalarImagePointerType weightImage = ScalarImageType::New();
    weightImage->CopyInformation(m_ReferenceImage);
    weightImage->SetRegions(m_ReferenceImage->GetRequestedRegion());
    weightImage->Allocate();*/
    
    typedef itk::ImageRegionConstIterator< ScalarImageType > WeightIteratorType;
    WeightIteratorType weightIter(totalWeightImage, totalWeightImage->GetRequestedRegion());
        
    for(fieldIter.GoToBegin(), weightIter.GoToBegin(); !fieldIter.IsAtEnd(); ++fieldIter, ++weightIter) {
        
        // get physical position
        index = fieldIter.GetIndex();
        estimatedField->TransformIndexToPhysicalPoint(index,point);
        value.Fill(0);
        
        // weight
        weightAll.fill(0);
        for(unsigned int i = 0; i < meanAll.size(); ++i) {
			for(unsigned int k = 0; k < NDimensions; ++k)
				x(k,0) = point[k]-meanAll[i][k];
			
            // test start
            //mahalanobisDistance = (x.transpose()*sigmaInverseAll[i]*x)(0,0);
            //if(mahalanobisDistance < 1.07) {
            // test end

            weightAll[i] = detSigma[i]*std::exp(-0.5*(x.transpose()*sigmaInverseAll[i]*x)(0,0));
                
            // test start
            //}
            // test end
            
        }
        
//        totalWeight = weightAll.sum();
        MatrixElementType totalWeight = weightIter.Get();
        if(totalWeight > 0)
            weightAll = weightAll/totalWeight;
        else
            std::cout << "weightAll is below zero" << std::endl;
        
        // debug
        //weightImage->SetPixel(index, totalWeight);
        
        // transform
        for(unsigned int i = 0; i < Mall.size(); ++i) {
            
            // test start
            //if(weightAll[i] > 0) {
            // test end
            
            for(unsigned int k = 0; k < NDimensions; ++k)
                xHom(k,0) = point[k];
            xHom(NDimensions,0) = 1;
            
            weightedVelocity = weightAll[i]*Mall[i]*xHom;
            for(unsigned int k = 0; k < NDimensions; ++k)
                value[k] += weightedVelocity(k,0);
                        
            // test start
            //}
            // test end
            
        }
        
        fieldIter.Set(value);
        
    }
    
    // debug
    /*weightWriter->SetFileName("WeightImage.mhd");
    weightWriter->SetInput(weightImage);
    weightWriter->Update();*/
    
    this->SetParametersAsVectorField(estimatedField);
    return this->m_VectorField;
}



template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
SetParametersAsVectorField(VectorFieldType * field)
{
    this->m_VectorField = field;
    this->m_InterpolateFunction->SetInputImage(this->m_VectorField);
    this->Modified();
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void 
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
SetParametersToNode(const GaussianTreeNodePointerType node, const ParametersType& allParameters, unsigned int& currentIndex, unsigned int currentLevel) 
{        
    ParametersType oneParameter(NDimensions*(NDimensions+1));
    for(unsigned int i = 0; i < oneParameter.GetSize(); ++i, ++currentIndex)
        oneParameter[i] = allParameters[currentIndex];
    node->Set(oneParameter);
    
    unsigned int lastLevel = log2( allParameters.size()/(NDimensions*(NDimensions+1)) + 1) - 1;    
    //std::cout << "currentLevel = " << currentLevel << " oneParameter = " << oneParameter << std::endl;
    
    if(currentLevel < lastLevel) {
        
        GaussianTreeNodeType::Pointer child1 = GaussianTreeNodeType::New();
        node->AddChild(child1);
        SetParametersToNode(child1, allParameters, currentIndex, currentLevel+1);

        GaussianTreeNodeType::Pointer child2 = GaussianTreeNodeType::New();
        node->AddChild(child2);
        SetParametersToNode(child2, allParameters, currentIndex, currentLevel+1);
        
    }    
}
    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
SetParametersByValue( const ParametersType & allParameters )
{        
    this->m_TransformationsRootNode = GaussianTreeNodeType::New();
    unsigned int currentIndex = 0;
    unsigned int currentLevel = 0;
    this->SetParametersToNode(this->m_TransformationsRootNode, allParameters, currentIndex, currentLevel);
    this->UpdateParameters(this->m_TransformationsRootNode, this->m_Parameters);
}
    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
SetFixedParameters( const ParametersType & allParameters )
{
    this->m_WeightsRootNode = GaussianTreeNodeType::New();
    unsigned int currentIndex = 0;
    unsigned int currentLevel = 0;
    this->SetParametersToNode(this->m_WeightsRootNode, allParameters, currentIndex, currentLevel);
    this->UpdateParameters(this->m_WeightsRootNode, this->m_FixedParameters);
}    

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void 
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
SetOutputParametersFromImage( const ReferenceImagePointerType referenceImage )
{
    this->m_ReferenceImage = referenceImage;
}

/*template<class TScalarType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, NDimensions>::VectorFieldType *
PolyaffineStationaryVelocityFieldTransform<TScalarType, NDimensions>::ConcatVelocityField(const VectorFieldType * field)
{
    typedef itk::MultiplyByConstantImageFilter< VectorFieldType, float, VectorFieldType > MultiplyFilterType;
    MultiplyFilterType::Pointer multiplyFilter1 = MultiplyFilterType::New();
    multiplyFilter1->SetInput( field );
    //multiplyFilter1->SetConstant( args.coefA );
    multiplyFilter1->Update();

    MultiplyFilterType::Pointer multiplyFilter2 = MultiplyFilterType::New();
    multiplyFilter2->SetInput( this->m_VectorField );
    //multiplyFilter2->SetConstant( args.coefB );
    multiplyFilter2->Update();

    typedef itk::VelocityFieldBCHCompositionFilter< VectorFieldType, VectorFieldType > BCHFilterType;
    BCHFilterType::Pointer bchFilter = BCHFilterType::New();
    bchFilter->SetInput( 0, multiplyFilter1->GetOutput() );
    bchFilter->SetInput( 1, multiplyFilter2->GetOutput() );
    //bchFilter->SetNumberOfApproximationTerms( args.nbBCHTerms );
    bchFilter->InPlaceOn();
    try
    {
        bchFilter->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
        std::cerr << excep << std::endl;

    }
     return bchFilter->GetOutput();

}*/

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
bool
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetInverse( Self* inverse ) const
{
    // set output parameters
    inverse->SetOutputParametersFromImage(this->m_ReferenceImage);
    
    // copy weights
    GaussianTreeNodePointerType weightNode = GaussianTreeNodeType::New();
    this->CopyTreeContent(this->GetWeights(),weightNode);
    inverse->SetWeights(weightNode);
    
    // copy transformations
    GaussianTreeNodePointerType transformNode = GaussianTreeNodeType::New();
    this->CopyTreeContent(this->GetTransformations(),transformNode);
    this->MultplyTreeByFactor(transformNode,-1);
    inverse->SetTransformations(transformNode);

    return true;
}



template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::InverseTransformBasePointer
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetInverseTransform() const
{
    Pointer inv = New();
    return GetInverse(inv) ? inv.GetPointer() : NULL;
}



template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::OutputPointType
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
TransformPoint(const InputPointType & point) const
{
    double minpixelspacing = m_VectorField->GetSpacing()[0];
    for (unsigned int i = 0; i<3; ++i)
    {
        if ( m_VectorField->GetSpacing()[i] < minpixelspacing/2 )
        {
            minpixelspacing = m_VectorField->GetSpacing()[i];
        }
    }

    typedef typename itk::ImageRegionIterator<VectorFieldType> IteratorType;

    IteratorType InputIt(const_cast<VectorFieldType*>(m_VectorField.GetPointer()), m_VectorField->GetRequestedRegion());

    float norm2,maxnorm2=0;
    for( InputIt.GoToBegin(); !InputIt.IsAtEnd(); ++InputIt )
    {
        norm2 = InputIt.Get().GetSquaredNorm();
        if (maxnorm2<norm2) maxnorm2=norm2;
    }

    maxnorm2 /= vnl_math_sqr(minpixelspacing);
    
    float numiterfloat = 2.0 +
            0.5 * vcl_log(maxnorm2)/vnl_math::ln2;
    std::cout<<"Num iter float: "<<numiterfloat<<std::endl;

    unsigned int constant=static_cast<unsigned int>(1<<static_cast<unsigned int>(abs(numiterfloat)));

    typedef typename itk::DivideByConstantImageFilter<VectorFieldType,float,VectorFieldType> DividerType;
    typename DividerType::Pointer Divider=DividerType::New();
    Divider->SetInput(m_VectorField);
    Divider->SetConstant( constant );
    std::cout<<"divider "<<  static_cast<unsigned int>(1<<static_cast<unsigned int>(abs(numiterfloat))) <<std::endl;

    Divider->Update();

    VectorFieldPointerType UpdatedVector = Divider->GetOutput();

    typedef typename itk::WarpVectorImageFilter<VectorFieldType,VectorFieldType,VectorFieldType> VectorWarperType;
    typename VectorWarperType::Pointer VectorWarper=VectorWarperType::New();

    VectorWarper->SetOutputOrigin(UpdatedVector->GetOrigin());
    VectorWarper->SetOutputSpacing(UpdatedVector->GetSpacing());
    VectorWarper->SetOutputDirection(UpdatedVector->GetDirection());


    OutputPointType output = point;

    for (unsigned int i =0; i < constant;++i)
    { 
        this->m_InterpolateFunction->SetInputImage(UpdatedVector);
        typename InterpolateFunctionType::OutputType vector = m_InterpolateFunction->Evaluate(output);

        for (unsigned int i=0; i<NDimensions; i++)
            output[i] += vector[i];

        VectorWarper->SetInput(UpdatedVector);
        VectorWarper->SetDisplacementField(UpdatedVector);
        //VectorWarper->SetDeformationField(UpdatedVector);
        VectorWarper->Update();

        UpdatedVector=VectorWarper->GetOutput();
        UpdatedVector->DisconnectPipeline();
    }

    return output;
}



template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::OutputVectorType
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
TransformVector(const InputVectorType & vector) const
{
    // Convert vector into point
    InputPointType point_0;
    for (unsigned int i=0; i<NDimensions; i++)
        point_0[i] = vector[i];

    // Transform point
    InputPointType point_1 = TransformPoint(point_0);

    // Convert point into vector
    OutputVectorType vector_1;
    for (unsigned int i=0; i<NDimensions; i++)
        vector_1[i] = point_1[i];
    return vector_1;
}



template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::OutputVnlVectorType
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
TransformVector(const InputVnlVectorType & vector) const
{
    // Convert vector into point
    InputPointType point_0;
    for (unsigned int i=0; i<NDimensions; i++)
        point_0[i] = vector[i];

    // Transform point
    InputPointType point_1 = TransformPoint(point_0);
    return point_1.GetVnlVector();
}



template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::VectorFieldPointerType
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetDisplacementFieldAsVectorField(const unsigned int noOfLevels /*typename SVFExponentialType::NumericalScheme scheme*/) 
{
    this->GetParametersAsVectorField(noOfLevels);
    
    // integrate to obtain the displacment vector field
    
    typename SVFExponentialType::Pointer exp = SVFExponentialType::New();
    exp->SetInput( this->m_VectorField );
    //exp->SetIterativeScheme( scheme );
    exp->Update();
    VectorFieldPointerType expField = exp->GetOutput();
    expField->DisconnectPipeline();
    return expField;
}



/*template<class TScalarType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, NDimensions>::ScalarImagePointerType
PolyaffineStationaryVelocityFieldTransform<TScalarType, NDimensions>::
GetLogSpatialJacobianDeterminant(typename SVFExponentialType::NumericalScheme scheme) const
{
    typename SVFExponentialType::Pointer exp = SVFExponentialType::New();
    exp->SetInput( this->m_VectorField );
    exp->SetIterativeScheme( scheme );
    exp->UpdateLargestPossibleRegion();
    return exp->GetLogJacobianDeterminant();
}*/



template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::OriginType
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetOrigin(void) const
{
    if (this->m_VectorField.IsNull())
        itkExceptionMacro("No field has been set.");
    return this->m_VectorField->GetOrigin();
}



template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::SpacingType
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetSpacing(void) const
{
    if (this->m_VectorField.IsNull())
        itkExceptionMacro("No field has been set.");
    return this->m_VectorField->GetSpacing();
}



template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::DirectionType
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetDirection(void) const
{
    if (this->m_VectorField.IsNull())
        itkExceptionMacro("No field has been set.");
    return this->m_VectorField->GetDirection();
}



template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::RegionType
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetLargestPossibleRegion(void) const
{
    if (this->m_VectorField.IsNull())
        itkExceptionMacro("No field has been set.");
    return this->m_VectorField->GetLargestPossibleRegion();
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
std::string 
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetTransformTypeAsString(void) const
{    
    std::ostringstream n;
    n << GetNameOfClass();
    n << "_";
    if ( typeid ( TScalarType ) == typeid ( float ) )
    {
        n << "float";
    }
    else if ( typeid ( TScalarType ) == typeid ( double ) )
    {
        n << "double";
    }
    else
    {
        n << "other";
    }
    n << "_" << this->GetOutputSpaceDimension()*(this->GetOutputSpaceDimension()+1);
    return n.str();
}
    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
CopyTreeStructure(GaussianTreeNodePointerType sourceNode, GaussianTreeNodePointerType targetNode) const
{
    // initialize current node value to zero
    NodeValueType value(sourceNode->Get().size());
    value.fill(0);
    targetNode->Set(value);

    // recursively go through tree
    GaussianTreeNodeType::ChildrenListType sourceChildren = sourceNode->GetChildrenList();
    for(unsigned int i = 0; i < sourceChildren.size(); ++i) {
        GaussianTreeNodePointerType targetChild = GaussianTreeNodeType::New();
        targetNode->AddChild(targetChild);
        this->CopyTreeStructure(sourceChildren[i], targetChild);
    }
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
CopyTreeContent(GaussianTreeNodePointerType sourceNode, GaussianTreeNodePointerType targetNode) const
{
    // copy current node value
    NodeValueType value = sourceNode->Get();
    targetNode->Set(value);

    // recursively go through tree
    GaussianTreeNodeType::ChildrenListType sourceChildren = sourceNode->GetChildrenList();
    for(unsigned int i = 0; i < sourceChildren.size(); ++i) {
        GaussianTreeNodePointerType targetChild = GaussianTreeNodeType::New();
        targetNode->AddChild(targetChild);
        this->CopyTreeContent(sourceChildren[i], targetChild);
    }
}
    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
MultplyTreeByFactor(GaussianTreeNodePointerType sourceNode, const MatrixElementType factor) const
{
    // multiply current node value by factor
    sourceNode->Set(factor*sourceNode->Get());
    
    // recursively go through tree
    GaussianTreeNodeType::ChildrenListType sourceChildren = sourceNode->GetChildrenList();
    for(unsigned int i = 0; i < sourceChildren.size(); ++i)
        this->MultplyTreeByFactor(sourceChildren[i], factor);
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
SetWeightsAtLevel(const unsigned int level, std::vector< GaussianTreeNodePointerType >& nodeSet) 
{
    std::vector< GaussianTreeNodePointerType > localNodeSet;
    this->GetWeightsAtLevel(level, localNodeSet);
    
    for(unsigned int i = 0; i < localNodeSet.size(); ++i) {
        localNodeSet[i]->Set(nodeSet[i]->Get());
    }
    
    this->UpdateParameters(this->m_WeightsRootNode, this->m_FixedParameters);
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
SetTransformationsAtLevel(const unsigned int level, std::vector< GaussianTreeNodePointerType >& nodeSet) 
{
    if(!m_TransformationsRootNode) {
        m_TransformationsRootNode = GaussianTreeNodeType::New();
        this->CopyTreeStructure(m_WeightsRootNode, m_TransformationsRootNode);
    }
    
    std::vector< GaussianTreeNodePointerType > localNodeSet;
    this->GetTransformationsAtLevel(level, localNodeSet);
    
    for(unsigned int i = 0; i < localNodeSet.size(); ++i) {
        localNodeSet[i]->Set(nodeSet[i]->Get());
    }
    
    this->UpdateParameters(this->m_TransformationsRootNode, this->m_Parameters);
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void 
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetNodesAtLevel(unsigned int currentLevel, const unsigned int lastLevel, GaussianTreeNodePointerType node, std::vector< GaussianTreeNodePointerType >& nodeSet) const
{
    if(currentLevel < lastLevel) {
        GaussianTreeNodeType::ChildrenListType children = node->GetChildrenList();
        for(unsigned int i = 0; i < children.size(); ++i)
            this->GetNodesAtLevel(currentLevel+1, lastLevel, children[i], nodeSet);
    }
    else
        nodeSet.push_back(node);
}
    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void 
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetWeightsAtLevel(const unsigned int level, std::vector< GaussianTreeNodePointerType >& nodeSet) const 
{
    unsigned int zeroLevel = 0;
    this->GetNodesAtLevel(zeroLevel, level, m_WeightsRootNode, nodeSet);
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void 
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetTransformationsAtLevel(const unsigned int level, std::vector< GaussianTreeNodePointerType >& nodeSet) const 
{
    unsigned int zeroLevel = 0;
    this->GetNodesAtLevel(zeroLevel, level, this->m_TransformationsRootNode, nodeSet);
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void 
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetResidualTransformationsAtLevel(const unsigned int level, std::vector< GaussianTreeNodePointerType >& nodeSet)
{
    unsigned int zeroLevel = 0;
    this->GetNodesAtLevel(zeroLevel, level, this->m_TransformationsRootNode, nodeSet);
    for(unsigned int i = 0; i < nodeSet.size(); ++i) {
        if(nodeSet[i]->HasParent())
            nodeSet[i]->Set( nodeSet[i]->Get() - nodeSet[i]->GetParent()->Get() );            
    }
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::GaussianTreeNodePointerType
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetWeights(void) const
{
    return this->m_WeightsRootNode;
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
SetWeights(const GaussianTreeNodePointerType rootNode)
{
    this->m_WeightsRootNode = rootNode;
    this->UpdateParameters(this->m_WeightsRootNode, this->m_FixedParameters);
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::GaussianTreeNodePointerType
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetTransformations(void) const
{
    return this->m_TransformationsRootNode;
}
    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
SetTransformations(const GaussianTreeNodePointerType rootNode)
{
    this->m_TransformationsRootNode = rootNode;
    this->UpdateParameters(this->m_TransformationsRootNode, this->m_Parameters);
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void 
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetParametersFromNode(const GaussianTreeNodePointerType node, ParametersType& parameters, unsigned int& currentIndex) const
{
    vnl_vector< double > current = node->Get();
    for(unsigned int i = 0; i < current.size(); ++i, ++currentIndex)
        parameters[currentIndex] = current[i];
    
    GaussianTreeNodeType::ChildrenListType children = node->GetChildrenList();
    for(unsigned int i = 0; i < children.size(); ++i)
        this->GetParametersFromNode(children[i], parameters, currentIndex);
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void 
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
CountNumberOfNodes(const GaussianTreeNodePointerType node, unsigned int& counter) const
{
    ++counter;
    GaussianTreeNodeType::ChildrenListType children = node->GetChildrenList();
    for(unsigned int i = 0; i < children.size(); ++i)
        this->CountNumberOfNodes(children[i], counter);
}
    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
unsigned int 
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetNumberOfLevels() const
{
    if(!m_WeightsRootNode) {
        std::cout << "no weights are defined, returning number of level = 0";
        return 0;
    }
    unsigned int counter = 0;
    this->CountNumberOfNodes(m_WeightsRootNode, counter);
    return log2(counter+1);
}
    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
const typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::ParametersType & 
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetParameters(void) const
{
    return m_Parameters;
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
const typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::ParametersType & 
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetFixedParameters(void) const
{
    return m_FixedParameters;
}
    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
UpdateParameters(const GaussianTreeNodePointerType node, ParametersType& parameters) 
{
    unsigned int numberOfNodes = 0;
    this->CountNumberOfNodes(node, numberOfNodes);
    parameters = ParametersType(numberOfNodes*NDimensions*(NDimensions+1));
    unsigned int currentIndex = 0;
    this->GetParametersFromNode(node, parameters, currentIndex);
}
    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
ParamToLogTransform(ParametersType params, vnl_matrix_fixed< MatrixElementType,NDimensionsHom,NDimensionsHom >& Ms) {
    
    Ms.fill(0.0);
    for(unsigned int m = 0; m < NDimensions; ++m) {
        for(unsigned int n = 0; n < NDimensionsHom; ++n) {
            Ms(m,n) = params(m*NDimensionsHom+n);
        }
    }
    
}
    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::ParametersType
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::LogTransformToParam(vnl_matrix_fixed< MatrixElementType,NDimensionsHom,NDimensionsHom >& Ms) {
    
    ParametersType params(NDimensions*NDimensions+NDimensions);
    
    for(unsigned int m = 0; m < NDimensions; ++m) {
        for(unsigned int n = 0; n < NDimensionsHom; ++n) {
            params(m*NDimensionsHom+n) = Ms(m,n);
        }
    }

    return params;
    
}
    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
ParamToGaussian(ParametersType params, vnl_vector_fixed< MatrixElementType,NDimensions >& regionCenter, vnl_matrix_fixed < MatrixElementType,NDimensions,NDimensions >& regionSigma) const {
    
    for(unsigned int i = 0; i < NDimensions; ++i)
        regionCenter[i] = params[i];
    
    for(unsigned int i = 0; i < NDimensions; ++i)
        for(unsigned int j = 0; j < NDimensions; ++j)
            regionSigma(i,j) = params[NDimensions+NDimensions*i+j];
    
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::ParametersType
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GaussianToParam(vnl_vector_fixed< MatrixElementType,NDimensions >& regionCenter, vnl_matrix_fixed < MatrixElementType,NDimensions,NDimensions >& regionSigma) {
        
        ParametersType params(NDimensions*NDimensions+NDimensions);
        
        for(unsigned int i = 0; i < NDimensions; ++i)
            params[i] = regionCenter[i];
        
        for(unsigned int i = 0; i < NDimensions; ++i)
            for(unsigned int j = 0; j < NDimensions; ++j)
                params[NDimensions+NDimensions*i+j] = regionSigma(i,j);
        
        return params;
}
    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
typename PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::ScalarImagePointerType
PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::    
GaussianToImage(vnl_vector_fixed< MatrixElementType,NDimensions >& mean, vnl_matrix_fixed < MatrixElementType,NDimensions,NDimensions >& covariance, ScalarImagePointerType referenceImage) {
    
    vnl_matrix< MatrixElementType > sigmaInverse = vnl_matrix_inverse< MatrixElementType >(covariance);
    
    ScalarImagePointerType weightImage = ScalarImageType::New();
    weightImage->CopyInformation(referenceImage);
    weightImage->SetRegions(referenceImage->GetRequestedRegion());
    weightImage->Allocate();
    
    typedef itk::ImageRegionIteratorWithIndex< ScalarImageType > ImageWithIndexIteratorType;
    ImageWithIndexIteratorType weightIter(weightImage, weightImage->GetRequestedRegion());
    for(weightIter.GoToBegin(); !weightIter.IsAtEnd(); ++weightIter) {
        
        // get physical position
        typename ScalarImageType::IndexType index = weightIter.GetIndex();
        typename ScalarImageType::PointType point;
        weightImage->TransformIndexToPhysicalPoint(index,point);
        
        vnl_matrix< MatrixElementType > x(NDimensions,1);
        for(unsigned int k = 0; k < NDimensions; ++k)
            x(k,0) = point[k]-mean[k];
        
        typename ScalarImageType::PixelType weight = std::exp(-0.5*(x.transpose()*sigmaInverse*x)(0,0));
        weightIter.Set(weight);
    }
    
    return weightImage;
}
    
    
} // namespace

#endif // _itkPolyaffineStationaryVelocityFieldTransform_txx_
