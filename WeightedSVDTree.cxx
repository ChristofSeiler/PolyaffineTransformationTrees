/*
 * ISTB & INRIA
 * Christof Seiler
 */

#include "WeightedSVDTree.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"

void WeightedSVDTree::computeGaussianWeight(OutputImageType::Pointer featureImage, double weightScaleParameter, vnl_vector_fixed< PolyaffineTreeTransformType::MatrixElementType,Dimension >& meanPoint, vnl_matrix_fixed < PolyaffineTreeTransformType::MatrixElementType,Dimension,Dimension >& covarianceMatrix) {
    
    meanPoint.fill(0);
    covarianceMatrix.fill(0);
    float informationMass = 0;
    
    typedef itk::ImageRegionIteratorWithIndex< OutputImageType > FeatureImageIterator;
    FeatureImageIterator featureIter(featureImage, featureImage->GetRequestedRegion());
    typedef itk::ImageRegionIterator< LabelImageType > MaskImageIterator;
    
    // compute information mass
    for(featureIter.GoToBegin(); !featureIter.IsAtEnd(); ++featureIter)
        informationMass += featureIter.Get();
    //std::cout << "informationMass = " << informationMass << std::endl;
    
    // compute weighted mean position
    for(featureIter.GoToBegin(); !featureIter.IsAtEnd(); ++featureIter) {
        if(featureIter.Get()) {
            
            OutputImageType::IndexType index = featureIter.GetIndex();
            OutputImageType::PointType point;
            featureImage->TransformIndexToPhysicalPoint(index,point);
            
            for(unsigned int i = 0; i < Dimension; ++i)
                meanPoint[i] += (featureIter.Get()/informationMass) * point[i];
            
        }
    }
    //std::cout << "meanPoint = " << meanPoint << std::endl;
    
    typedef itk::ImageFileWriter< OutputImageType > FeatureWriterType;
    FeatureWriterType::Pointer featureWriter = FeatureWriterType::New();
    featureWriter->SetInput(featureImage);
    std::ostringstream osstrFeature;
    osstrFeature << "CenterX" << meanPoint[0] << "Y" << meanPoint[1] << "Z" << meanPoint[2] << ".mhd";
    featureWriter->SetFileName(osstrFeature.str().c_str());
    //featureWriter->Update();
    
    // compute weighed covariance matrix
    for(featureIter.GoToBegin(); !featureIter.IsAtEnd(); ++featureIter) {
        if(featureIter.Get()) {
            
            OutputImageType::IndexType index = featureIter.GetIndex();
            OutputImageType::PointType point;
            featureImage->TransformIndexToPhysicalPoint(index,point);
            
            vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > centeredPoint(Dimension,1);
            for(unsigned int i = 0; i < Dimension; ++i) {
                centeredPoint(i,0) = (featureIter.Get()/informationMass) * (point[i] - meanPoint[i]);
            }
            
            covarianceMatrix += centeredPoint * centeredPoint.transpose();
            
        }
    }
    
    vnl_svd< PolyaffineTreeTransformType::MatrixElementType > svd(covarianceMatrix);
    //std::cout << "svd = " << svd << std::endl;
    
    //compute extend of boxes
    OutputImageType::PointType maxExtent[3];
    OutputImageType::PointType minExtent[3];
    OutputImageType::PixelType maxValue[3];
    OutputImageType::PixelType minValue[3];
    for(unsigned int i = 0; i < Dimension; ++i) {
        maxExtent[i].Fill(0);
        minExtent[i].Fill(0);
        maxValue[i] = 0;
        minValue[i] = 0;
    }
    for(featureIter.GoToBegin(); !featureIter.IsAtEnd(); ++featureIter) {
        if(featureIter.Get()) {
            
            OutputImageType::IndexType index = featureIter.GetIndex();
            OutputImageType::PointType point;
            featureImage->TransformIndexToPhysicalPoint(index,point);
            
            vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > centeredPoint(Dimension,1);
            for(unsigned int i = 0; i < Dimension; ++i)
                centeredPoint(i,0) = (featureIter.Get()/informationMass) * (point[i] - meanPoint[i]);
            
            vnl_matrix< PolyaffineTreeTransformType::MatrixElementType> projection = centeredPoint.transpose() * svd.U();
            
            for(unsigned int i = 0; i < Dimension; ++i) {
                if(projection(0,i) > maxValue[i]) {
                    maxValue[i] = projection(0,i);
                	maxExtent[i] = point;
                }
                if(projection(0,i) < minValue[i]) {
                    minValue[i] = projection(0,i);
                    minExtent[i] = point;
                }
            }
            
        }
    }
    for(unsigned int i = 0; i < Dimension; ++i) {
        float stdev = weightScaleParameter*0.5*std::sqrt(minExtent[i].SquaredEuclideanDistanceTo(maxExtent[i]));
        //std::cout << "stdev " << i << " = " << stdev << std::endl;
        svd.W(i,i) = stdev*stdev;
    }
    covarianceMatrix = svd.recompose();
    //std::cout << "covarianceMatrix = " << covarianceMatrix << std::endl;

}

void WeightedSVDTree::weightTreeEstimation(unsigned int currentLevel, const unsigned int lastLevel, double weightScaleParameter, OutputImageType::Pointer featureImage, PolyaffineTreeTransformType::Pointer polyaffineTransform, GaussianTreeNodeType::Pointer currNode) {
        
    // weighted PCA
    vnl_vector_fixed< PolyaffineTreeTransformType::MatrixElementType,Dimension > meanPoint;
    vnl_matrix_fixed< PolyaffineTreeTransformType::MatrixElementType,Dimension,Dimension > covarianceMatrix;

    // comptue Gaussian weight
    this->computeGaussianWeight(featureImage, weightScaleParameter, meanPoint, covarianceMatrix);
    
    //std::cout << "scaledCovarianceMatrix = " << scaledCovarianceMatrix << std::endl;
    //std::cout << "scaled svd = " << svd << std::endl;
    
    NodeValueType estimatedParams = polyaffineTransform->GaussianToParam(meanPoint, covarianceMatrix);
    currNode->Set(estimatedParams);
    
    OutputImageType::Pointer weightedCovarianceImage = polyaffineTransform->GaussianToImage(meanPoint, covarianceMatrix, featureImage);
    typedef itk::ImageFileWriter< OutputImageType > WeightWriterType;
    WeightWriterType::Pointer weightWriter = WeightWriterType::New();
    weightWriter->SetInput(weightedCovarianceImage);
    std::ostringstream osstr;
    osstr << "WeightedCovarianceImage" << currentLevel << "CenterX" << meanPoint[0] << "Y" << meanPoint[1] << "Z" << meanPoint[2] << ".mhd";
    weightWriter->SetFileName(osstr.str().c_str());
    //weightWriter->Update();
        
    if(currentLevel < lastLevel) {
    
        // cluster        
        OutputImageType::Pointer featureImage1 = OutputImageType::New();
        featureImage1->CopyInformation(featureImage);
        featureImage1->SetRegions(featureImage->GetRequestedRegion());
        featureImage1->Allocate();
        featureImage1->FillBuffer(0);

        OutputImageType::Pointer featureImage2 = OutputImageType::New();
        featureImage2->CopyInformation(featureImage);
        featureImage2->SetRegions(featureImage->GetRequestedRegion());
        featureImage2->Allocate();
        featureImage2->FillBuffer(0);
        
        typedef itk::ImageRegionIterator< OutputImageType > OutputImageIterator;
        OutputImageIterator featureImageIter1(featureImage1, featureImage1->GetRequestedRegion());
        OutputImageIterator featureImageIter2(featureImage2, featureImage2->GetRequestedRegion());
        
        typedef itk::ImageRegionIteratorWithIndex< OutputImageType > FeatureImageIterator;
        FeatureImageIterator featureIter(featureImage, featureImage->GetRequestedRegion());

        for(featureIter.GoToBegin(), featureImageIter1.GoToBegin(), featureImageIter2.GoToBegin(); !featureIter.IsAtEnd(); ++featureIter, ++featureImageIter1, ++featureImageIter2) {
            
            if(featureIter.Get()) {
                
                OutputImageType::IndexType index = featureIter.GetIndex();
                OutputImageType::PointType point;
                featureImage->TransformIndexToPhysicalPoint(index,point);
                
                vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > centeredPoint(Dimension,1);
                for(unsigned int i = 0; i < Dimension; ++i)
                    centeredPoint(i,0) = point[i] - meanPoint[i];
                
                vnl_svd< PolyaffineTreeTransformType::MatrixElementType > svd(covarianceMatrix);
                vnl_vector< PolyaffineTreeTransformType::MatrixElementType > projection = centeredPoint.transpose() * svd.U().get_column(0);
                if( projection(0) < 0.0 )
                    featureImageIter1.Set(featureIter.Get());
                else
                    featureImageIter2.Set(featureIter.Get());
                
            }
        }
        
        GaussianTreeNodeType::Pointer child1 = GaussianTreeNodeType::New();
        currNode->AddChild(child1);
        weightTreeEstimation(currentLevel+1, lastLevel, weightScaleParameter, featureImage1, polyaffineTransform, child1);
        
        GaussianTreeNodeType::Pointer child2 = GaussianTreeNodeType::New();
        currNode->AddChild(child2);
        weightTreeEstimation(currentLevel+1, lastLevel, weightScaleParameter, featureImage2, polyaffineTransform, child2);
        
    }
    
}
