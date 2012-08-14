/*
 * ISTB & INRIA
 * Christof Seiler
 */

#ifndef WEIGHTEDSVDTREE_H
#define WEIGHTEDSVDTREE_H

#include "itkPolyaffineStationaryVelocityFieldTransform.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTransformFileReader.h>
#include <itkWarpImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkVector.h>
#include <itkListSample.h>
#include <itkGaussianMixtureModelComponent.h>
#include <itkExpectationMaximizationMixtureModelEstimator.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkDivideByConstantImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkTreeNode.h>
#include <itkGaussianSpatialFunction.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkLabelOverlayImageFilter.h>
#include <itkBoxMeanImageFilter.h>

class WeightedSVDTree {

public:
    
    static const unsigned int Dimension = 3;
    
    typedef signed short PixelType;
    typedef itk::Image< PixelType, Dimension > ImageType;
    typedef signed short LabelPixelType;
    typedef itk::Image< LabelPixelType, Dimension > LabelImageType;
    typedef itk::Image< float,Dimension > OutputImageType;
    
    typedef itk::Array< double > ParametersType;
    static const unsigned int numberOfClasses = 2;
    
    typedef itk::PolyaffineStationaryVelocityFieldTransform<> PolyaffineTreeTransformType;
    typedef vnl_vector< double > NodeValueType;
    typedef itk::TreeNode< NodeValueType > GaussianTreeNodeType;
                
    void weightTreeEstimation(unsigned int currentLevel, const unsigned int lastLevel, double weightScaleParameter, OutputImageType::Pointer featureImage, PolyaffineTreeTransformType::Pointer polyaffineTransform, GaussianTreeNodeType::Pointer currNode);
    
    void computeGaussianWeight(OutputImageType::Pointer featureImage, double weightScaleParameter, vnl_vector_fixed< PolyaffineTreeTransformType::MatrixElementType,Dimension >& meanPoint, vnl_matrix_fixed < PolyaffineTreeTransformType::MatrixElementType,Dimension,Dimension >& scaledCovarianceMatrix);

};

#endif
