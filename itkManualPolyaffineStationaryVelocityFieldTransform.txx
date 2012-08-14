#ifndef _itkManualPolyaffineStationaryVelocityFieldTransform_txx_
#define _itkManualPolyaffineStationaryVelocityFieldTransform_txx_

#include "itkManualPolyaffineStationaryVelocityFieldTransform.h"

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
ManualPolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
ManualPolyaffineStationaryVelocityFieldTransform() : Superclass( )
{
}

template<class TScalarType, class TReferenceType, unsigned int NDimensions>
void 
ManualPolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
SetParametersToNode(const GaussianTreeNodePointerType node, const ParametersType& allParameters, unsigned int& currentIndex, unsigned int currentLevel) 
{
    ParametersType oneParameter(NDimensions*(NDimensions+1));
    
    // get root
    for(unsigned int i = 0; i < oneParameter.GetSize(); ++i, ++currentIndex)
        oneParameter[i] = allParameters[currentIndex];
    node->Set(oneParameter);
    
    // first level: manual regions
    unsigned int noOfManualRegions = (allParameters.GetSize() / oneParameter.GetSize()) - 1;
    for(unsigned int i = 0; i < noOfManualRegions; ++i) {
        
        GaussianTreeNodeType::Pointer child = GaussianTreeNodeType::New();
        node->AddChild(child);
        
        for(unsigned int i = 0; i < oneParameter.GetSize(); ++i, ++currentIndex)
            oneParameter[i] = allParameters[currentIndex];
        child->Set(oneParameter);

    }    
}
    
template<class TScalarType, class TReferenceType, unsigned int NDimensions>
unsigned int 
ManualPolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>::
GetNumberOfLevels() const
{
    return 2;
}
    
} // namespace

#endif // _itkManualPolyaffineStationaryVelocityFieldTransform_txx_
