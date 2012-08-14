/**
 * This is a small tool that shows how to use the polyaffine log-domain demons algorithm.
 * Christof Seiler, ISTB, University of Bern and Asclepios Team, INRIA Sophia Antipolis
 */

#include "itkPolyaffineLogDomainDeformableRegistrationFilter.h"
#include "itkPolyaffineLogDomainDemonsRegistrationFilter.h"
#include "itkPolyaffineSymmetricLogDomainDemonsRegistrationFilter.h"
//#include "vtkOBBTreeCustomized.h"
#include "itkPolyaffineStationaryVelocityFieldTransform.h"
#include "itkManualPolyaffineStationaryVelocityFieldTransform.h"
#include "WeightedSVDTree.h"

#include <itkCommand.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTransformFileReader.h>
#include <itkTransformToVelocityFieldSource.h>
#include <itkWarpImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkWarpHarmonicEnergyCalculator.h>
#include <itkMaskImageFilter.h>
#include <itkLabelContourImageFilter.h>
#include <itkTransformFileWriter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkMRIBiasFieldCorrectionFilter.h>
#include <itkSquareImageFilter.h>
#include <itkAddConstantToImageFilter.h>
#include <itkAffineTransform.h>
#include <itkHistogramMatchingImageFilter.h>
#include <itkRegularExpressionSeriesFileNames.h>
#include <itkTransformFactory.h>
#include <itkMinimumMaximumImageFilter.h>

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#define itkProbesCreate()  \
itk::TimeProbesCollectorBase chronometer; \
itk::MemoryProbesCollectorBase memorymeter
#define itkProbesStart( text ) memorymeter.Start( text ); chronometer.Start( text )
#define itkProbesStop( text )  chronometer.Stop( text ); memorymeter.Stop( text  )
#define itkProbesReport( stream )  chronometer.Report( stream ); memorymeter.Report( stream  )

//#include <vtkSmartPointer.h>
//#include <vtkMarchingCubes.h>
//#include <vtkFloatArray.h>
//#include <vtkPolyDataWriter.h>
//#include <vtkPointData.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkLookupTable.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkProperty.h>
//#include <vtkTextProperty.h>
//#include <vtkTextActor.h>
//#include <vtkCamera.h>
//#include <vtkWindowToImageFilter.h>
//#include <vtkPNGWriter.h>
//#include <vtkLine.h>
//#include <vtkCellArray.h>
//#include <vtkCellData.h>
//#include <vtkSTLWriter.h>

#include <metaCommand.h>

#include <errno.h>
#include <limits.h>

#include <vnl/vnl_trace.h>

const unsigned int Dimension = 3;

typedef signed short PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::Vector<float,Dimension> VectorType;
typedef itk::Image<VectorType,Dimension> FieldType;
typedef signed short LabelPixelType;
typedef itk::Image< LabelPixelType, Dimension > LabelImageType;
typedef itk::AffineTransform< double,Dimension > AffineTransformType;

//const ImageType::PixelType OUTSIDE_VOXEL = -1000;
//const LabelImageType::PixelType OUTSIDE_LABEL = 0;
//const LabelImageType::PixelType INSIDE_LABEL = 1;

typedef itk::Image< float,Dimension > WeightImageType;

//const double currweight = 0.75;

struct arguments
{
	std::string  fixedImageFile;  /* -f option */
	std::string  movingImageFile; /* -m option */
	std::string  fixedmaskImageFile; /* -F option */
    std::string  manualRegionsLabelFile;
	std::string  movingmaskImageFile; /* -M option */
	std::string  inputFieldFile;  /* -b option */
    std::string  priorMui;
    std::string  priorSigmaInverse;
//    std::string  priorTransformDirectoryC;
//    std::string  priorTransformDirectoryLambda;
	std::string  outputImageFile; /* -o option */
	std::string  outputDeformationFieldFile;
	std::string  outputInverseDeformationFieldFile;
	std::string  outputVelocityFieldFile;
    std::string  outputPolyaffineTreeTransformFile;
//    std::string  outputLambdaNewFile;
//    std::string  outputCNewFile;
	std::string  outputCorrelationMatrixFile;
	std::string  outputMetricFile;
	std::string  outputHarmonicEnergyFile;
	float maxStepLength;          /* -l option */
	unsigned int updateRule;      /* -a option */
	float weightScalingValue;	/* -d option */
	unsigned int numberOfLevels; /* -s option */
	double varianceVelocityField; /* -w option */
    double unionMaskDilationRadius; /* -r option */
//    double scaleLambdaPrior; /* -z option */
 	unsigned int gradientType;    /* -t option */
	unsigned int numberOfBCHApproximationTerms; /* -c option */
	unsigned int leastSquaresType; /* -x option */
	
	friend std::ostream& operator<< (std::ostream& o, const arguments& args)
    {				
		std::string gtypeStr;
		switch (args.gradientType)
		{
			case 0:
				gtypeStr = "symmetrized (ESM for diffeomorphic and compositive)";
				break;
			case 1:
				gtypeStr = "fixed image (Thirion's vanilla forces)";
				break;
			case 2:
				gtypeStr = "warped moving image (Gauss-Newton for diffeomorphic and compositive)";
				break;
			case 3:
				gtypeStr = "mapped moving image (Gauss-Newton for additive)";
				break;
			default:
				gtypeStr = "unsuported";
		}

		std::string uruleStr;
		switch (args.updateRule)
		{
			case 0:
				uruleStr = "BCH approximation on velocity fields (log-domain)";
				break;
			case 1:
				uruleStr = "Symmetrized BCH approximation on velocity fields (symmetric log-domain)";
				break;
			default:
				uruleStr = "unsuported";
		}
		
		std::string lstypeStr;
		switch (args.leastSquaresType)
		{
			case 0:
				lstypeStr = "log(gradient + 1)";
				break;
			case 1:
				lstypeStr = "contour";
				break;
            case 2:
                lstypeStr = "mask";
				break;
			default:
				lstypeStr = "unsuported";
		}		
						
		return o
		<<"Arguments structure:"<<std::endl
		<<"  Fixed image file: "<<args.fixedImageFile<<std::endl
		<<"  Moving image file: "<<args.movingImageFile<<std::endl
		<<"  Fixed mask image file: "<<args.fixedmaskImageFile<<std::endl
		<<"  Moving mask image file: "<<args.movingmaskImageFile<<std::endl
        <<"  Label image defining manual regions: "<<args.manualRegionsLabelFile<<std::endl
		<<"  Input velocity field file: "<<args.inputFieldFile<<std::endl
        <<"  Prior mean: "<<args.priorMui<<std::endl 
        <<"  Prior SigmaInverse: "<<args.priorSigmaInverse<<std::endl 
//        <<"  Prior transform directory for C: "<<args.priorTransformDirectoryC<<std::endl 
//        <<"  Prior transform directory for Lambda: "<<args.priorTransformDirectoryLambda<<std::endl 
		<<"  Output image file: "<<args.outputImageFile<<std::endl
		<<"  Output deformation field file: "<<args.outputDeformationFieldFile<<std::endl
		<<"  Output inverse deformation field file: "<<args.outputInverseDeformationFieldFile<<std::endl
		<<"  Output velocity field file: "<<args.outputVelocityFieldFile<<std::endl
        <<"  Output polyaffine tree transform file: "<<args.outputPolyaffineTreeTransformFile<<std::endl
//        <<"  Output new Lambda file: "<<args.outputLambdaNewFile<<std::endl
//        <<"  Output new C file: "<<args.outputCNewFile<<std::endl
		<<"  Output correlation matrix file: "<<args.outputCorrelationMatrixFile<<std::endl
		<<"  Output metric file: "<<args.outputMetricFile<<std::endl
		<<"  Output harmonic energy file: "<<args.outputHarmonicEnergyFile<<std::endl
		<<"  Maximum update step length: "<<args.maxStepLength<<std::endl
		<<"  Update rule: "<<uruleStr<<std::endl
		<<"  Type of gradient: "<<gtypeStr<<std::endl
		<<"  Weight scaling factor: "<<args.weightScalingValue<<std::endl
		<<"  Number of terms in the BCH expansion: "<<args.numberOfBCHApproximationTerms<<std::endl
		<<"  Number of tree levels: "<<args.numberOfLevels<<std::endl
		<<"  Type of least squares: "<<lstypeStr<<std::endl
		<<"  Variance velocity field: "<<args.varianceVelocityField<<std::endl
        <<"  Union mask dilation radius: "<<args.unionMaskDilationRadius<<std::endl;
//        <<"  Scale lambda prior: "<<args.scaleLambdaPrior<<std::endl;
    }
};

void help_callback()
{
	std::cout<<std::endl;
	std::cout<<"Copyright (c) 2011 University of Bern and INRIA"<<std::endl;
	std::cout<<"Code: Christof Seiler"<<std::endl;
	std::cout<<"Report bugs to <christof.seiler \\at istb.unibe.ch>"<<std::endl;
	
	exit( EXIT_FAILURE );
};


void parseOpts (int argc, char **argv, struct arguments & args)
{
	// Command line parser
	MetaCommand command;
	command.SetParseFailureOnUnrecognizedOption( true );
	command.SetHelpCallBack(help_callback);
	
	// Fill some information about the software
	command.SetAuthor("Christof Seiler");
	
	command.SetAcknowledgments("This work is supported by University of Bern (ISTB) and INRIA (Asclepios team)");
	
	command.SetDescription("Basic image registration tool with the polyaffine log-domain demons algorithm.");
	
	// Define parsing options
	command.SetOption("FixedImageFile","f",true,"Fixed image filename");
	command.SetOptionLongTag("FixedImageFile","fixed-image");
	command.AddOptionField("FixedImageFile","filename",MetaCommand::STRING,true);
	
	command.SetOption("MovingImageFile","m",true,"Moving image filename");
	command.SetOptionLongTag("MovingImageFile","moving-image");
	command.AddOptionField("MovingImageFile","filename",MetaCommand::STRING,true);
	
	command.SetOption("FixedMaskImageFile","F",true,"Fixed mask image filename");
	command.SetOptionLongTag("FixedMaskImageFile","fixed-mask-image");
	command.AddOptionField("FixedMaskImageFile","filename",MetaCommand::STRING,true);
	
    command.SetOption("MovingMaskImageFile","M",true,"Moving mask image filename");
	command.SetOptionLongTag("MovingMaskImageFile","moving-mask-image");
	command.AddOptionField("MovingMaskImageFile","filename",MetaCommand::STRING,true);

    command.SetOption("ManualRegionsLabelFile","R",false,"Label image defining manual regions");
	command.SetOptionLongTag("ManualRegionsLabelFile","label-image-manual-regions");
	command.AddOptionField("ManualRegionsLabelFile","filename",MetaCommand::STRING,true);
		
	command.SetOption("InputFieldFile","b",false,"Input velocity field filename");
	command.SetOptionLongTag("InputFieldFile","input-field");
	command.AddOptionField("InputFieldFile","filename",MetaCommand::STRING,true);

    command.SetOption("PriorMui","",false,"mu_i file name");
    command.SetOptionLongTag("PriorMui","prior-mu_i");
    command.AddOptionField("PriorMui","filename",MetaCommand::STRING,true);
    
    command.SetOption("PriorSigmaInverse","",false,"SigmaInverse file name");
    command.SetOptionLongTag("PriorSigmaInverse","prior-SigmaInverse");
    command.AddOptionField("PriorSigmaInverse","filename",MetaCommand::STRING,true);
    
//    command.SetOption("PriorTransformDirectoryC","",false,"Transformations from prior C registration directory name");
//	command.SetOptionLongTag("PriorTransformDirectoryC","prior-C-directory");
//	command.AddOptionField("PriorTransformDirectoryC","filename",MetaCommand::STRING,true);
//
//    command.SetOption("PriorTransformDirectoryLambda","",false,"Transformations from prior Lambda registration directory name");
//	command.SetOptionLongTag("PriorTransformDirectoryLambda","prior-Lambda-directory");
//	command.AddOptionField("PriorTransformDirectoryLambda","filename",MetaCommand::STRING,true);
    
	command.SetOption("OutputImageFile","o",false,"Output image filename");
	command.SetOptionLongTag("OutputImageFile","output-image");
	command.AddOptionField("OutputImageFile","filename",MetaCommand::STRING,true,"output.mha");
	
	command.SetOption("OutputDeformationFieldFile","",false,"Output deformation field filename");
	command.SetOptionLongTag("OutputDeformationFieldFile","outputDef-field");
	command.AddOptionField("OutputDeformationFieldFile","filename",MetaCommand::STRING,false,"OUTPUTIMAGENAME-deformationField.mha");
	
	command.SetOption("OutputInverseDeformationFieldFile","",false,"Output inverse deformation field filename");
	command.SetOptionLongTag("OutputInverseDeformationFieldFile","outputInvDef-field");
	command.AddOptionField("OutputInverseDeformationFieldFile","filename",MetaCommand::STRING,false,"OUTPUTIMAGENAME-inverseDeformationField.mha");
	
	command.SetOption("OutputVelocityFieldFile","",false,"Output velocity field filename");
	command.SetOptionLongTag("OutputVelocityFieldFile","outputVel-field");
	command.AddOptionField("OutputVelocityFieldFile","filename",MetaCommand::STRING,false,"OUTPUTIMAGENAME-velocityField.mha");

    command.SetOption("OutputPolyaffineTreeTransformFile","",false,"Output polyaffine tree transform filename");
	command.SetOptionLongTag("OutputPolyaffineTreeTransformFile","output-polyaffine-tree-transform");
	command.AddOptionField("OutputPolyaffineTreeTransformFile","filename",MetaCommand::STRING,false,"OUTPUTIMAGENAME-polyaffineTreeTransform.tfm");

//    command.SetOption("OutputLambdaNewFile","",false,"Output new Lambda filename");
//	command.SetOptionLongTag("OutputLambdaNewFile","output-Lambda-new");
//	command.AddOptionField("OutputLambdaNewFile","filename",MetaCommand::STRING,false,"OUTPUTIMAGENAME-LambdaNew.txt");
//
//    command.SetOption("OutputCNewFile","",false,"Output new C filename");
//	command.SetOptionLongTag("OutputCNewFile","output-C-new");
//	command.AddOptionField("OutputCNewFile","filename",MetaCommand::STRING,false,"OUTPUTIMAGENAME-CNew.txt");
		
	command.SetOption("OutputCorrelationMatrixFile","",false,"Output correlation matrix filename");
	command.SetOptionLongTag("OutputCorrelationMatrixFile","output-correlation-matrix");
	command.AddOptionField("OutputCorrelationMatrixFile","filename",MetaCommand::STRING,false,"OUTPUTIMAGENAME-correlationMatrix-LEVEL.txt");

	command.SetOption("OutputMetricFile","",false,"Output metric filename");
	command.SetOptionLongTag("OutputMetricFile","output-metric");
	command.AddOptionField("OutputMetricFile","filename",MetaCommand::STRING,false,"OUTPUTIMAGENAME-metric.txt");
	
	command.SetOption("OutputHarmonicEnergyFile","",false,"Output harmonic energy filename");
	command.SetOptionLongTag("OutputHarmonicEnergyFile","output-harmonic-energy");
	command.AddOptionField("OutputHarmonicEnergyFile","filename",MetaCommand::STRING,false,"OUTPUTIMAGENAME-harmonic-energy.txt");

	command.SetOption("MaximumUpdateStepLength","l",false,"Maximum length of an update vector (pixel units). Setting it to 0 implies no restrictions will be made on the step length");
	command.SetOptionLongTag("MaximumUpdateStepLength","max-step-length");
	command.AddOptionField("MaximumUpdateStepLength","floatval",MetaCommand::FLOAT,true,"2.0");		
	
	command.SetOption("UpdateRule","a",false,"Type of update rule. 0: exp(v) <- exp(v) o exp(u) (log-domain), 1: exp(v) <- symmetrized( exp(v) o exp(u) ) (symmetric log-domain)");
	command.SetOptionLongTag("UpdateRule","update-rule");
	command.AddOptionField("UpdateRule","type",MetaCommand::INT,true,"1");
	command.SetOptionRange("UpdateRule","type","0","1");
	
	command.SetOption("GradientType","t",false,"Type of gradient used for computing the demons force. 0 is symmetrized, 1 is fixed image, 2 is warped moving image, 3 is mapped moving image");
	command.SetOptionLongTag("GradientType","gradient-type");
	command.AddOptionField("GradientType","type",MetaCommand::INT,true,"0");
	command.SetOptionRange("GradientType","type","0","3");
	
	command.SetOption("NumberOfBCHApproximationTerms","c",false,"Number of terms in the BCH expansion");
	command.SetOptionLongTag("NumberOfBCHApproximationTerms","num-bch-terms");
	command.AddOptionField("NumberOfBCHApproximationTerms","intval",MetaCommand::INT,true,"2");
	command.SetOptionRange("NumberOfBCHApproximationTerms","intval","2","4");

	command.SetOption("NumberOfLevels","s",false,"Number of levels");
	command.SetOptionLongTag("NumberOfLevels","number-of-level");
	command.AddOptionField("NumberOfLevels","intval",MetaCommand::INT,true,"1");
	
	command.SetOption("WeightScalingValue","d",false,"Weight scaling value");
	command.SetOptionLongTag("WeightScalingValue","weight-scaling-value");
	command.AddOptionField("WeightScalingValue","floatval",MetaCommand::FLOAT,true,"1");	
	
	command.SetOption("LeastSquaresType","x",false,"Type of feature. (0) log(gardient + 1), (1) contour, (2) mask");
	command.SetOptionLongTag("LeastSquaresType","feature-type");
	command.AddOptionField("LeastSquaresType","type",MetaCommand::INT,true,"0");
	command.SetOptionRange("LeastSquaresType","type","0","2");
	
	command.SetOption("VarianceVelocityField","w",false,"Variance velocity field");
	command.SetOptionLongTag("VarianceVelocityField","variance-velocity-field");
	command.AddOptionField("VarianceVelocityField","floatval",MetaCommand::FLOAT,true,"1");
    
	command.SetOption("UnionMaskDilationRadius","r",false,"Union mask dilation radius");
	command.SetOptionLongTag("UnionMaskDilationRadius","union-mask-dilation-radius");
	command.AddOptionField("UnionMaskDilationRadius","floatval",MetaCommand::FLOAT,true,"0");

//    command.SetOption("ScaleLambdaPrior","z",false,"Scale lambda coefficients for prior");
//	command.SetOptionLongTag("ScaleLambdaPrior","scale-lambda-prior");
//	command.AddOptionField("ScaleLambdaPrior","floatval",MetaCommand::FLOAT,true,"0.01");
    
	// Actually parse the command line
	if (!command.Parse(argc,argv))
    {
		exit( EXIT_FAILURE );
    }
	
	
	// Store the parsed information into a struct
	args.fixedImageFile = command.GetValueAsString("FixedImageFile","filename");
	args.movingImageFile = command.GetValueAsString("MovingImageFile","filename");
	args.fixedmaskImageFile = command.GetValueAsString("FixedMaskImageFile","filename");
	args.movingmaskImageFile = command.GetValueAsString("MovingMaskImageFile","filename");
    args.manualRegionsLabelFile = command.GetValueAsString("ManualRegionsLabelFile","filename");

	args.inputFieldFile = command.GetValueAsString("InputFieldFile","filename");
    args.priorMui = command.GetValueAsString("PriorMui","filename");
    args.priorSigmaInverse = command.GetValueAsString("PriorSigmaInverse","filename");
//    args.priorTransformDirectoryC = command.GetValueAsString("PriorTransformDirectoryC","filename");
//    args.priorTransformDirectoryLambda = command.GetValueAsString("PriorTransformDirectoryLambda","filename");
	args.outputImageFile = command.GetValueAsString("OutputImageFile","filename");
    	
	args.outputDeformationFieldFile = command.GetValueAsString("OutputDeformationFieldFile","filename");
	args.outputInverseDeformationFieldFile = command.GetValueAsString("OutputInverseDeformationFieldFile","filename");
	args.outputVelocityFieldFile = command.GetValueAsString("OutputVelocityFieldFile","filename");
    args.outputPolyaffineTreeTransformFile = command.GetValueAsString("OutputPolyaffineTreeTransformFile","filename");
//    args.outputLambdaNewFile = command.GetValueAsString("OutputLambdaNewFile","filename");
//    args.outputCNewFile = command.GetValueAsString("OutputCNewFile","filename");
	args.outputCorrelationMatrixFile = command.GetValueAsString("OutputCorrelationMatrixFile","filename");
	args.outputMetricFile = command.GetValueAsString("OutputMetricFile","filename");
	args.outputHarmonicEnergyFile = command.GetValueAsString("OutputHarmonicEnergyFile", "filename");
	
	unsigned int pos = args.outputImageFile.rfind(".");
	
	// Change the extension by -deformationField.mha
	if ( args.outputDeformationFieldFile == "OUTPUTIMAGENAME-deformationField.mha" )
	{
		if ( pos < args.outputDeformationFieldFile.size() )
		{
			args.outputDeformationFieldFile = args.outputImageFile;
			args.outputDeformationFieldFile.replace(pos, args.outputDeformationFieldFile.size(), "-deformationField.mha");
		}
		else
		{
			args.outputDeformationFieldFile = args.outputImageFile + "-deformationField.mha";
		}
	}
	
	// Change the extension by -inverseDeformationField.mha
	if ( args.outputInverseDeformationFieldFile == "OUTPUTIMAGENAME-inverseDeformationField.mha" )
	{
		if ( pos < args.outputInverseDeformationFieldFile.size() )
		{
			args.outputInverseDeformationFieldFile = args.outputImageFile;
			args.outputInverseDeformationFieldFile.replace(pos, args.outputInverseDeformationFieldFile.size(), "-inverseDeformationField.mha");
		}
		else
		{
			args.outputInverseDeformationFieldFile = args.outputImageFile + "-inverseDeformationField.mha";
		}
	}
	
	// Change the extension by -velocityField.mha
	if ( args.outputVelocityFieldFile == "OUTPUTIMAGENAME-velocityField.mha" )
    {
		if ( pos < args.outputVelocityFieldFile.size() )
		{
			args.outputVelocityFieldFile = args.outputImageFile;
			args.outputVelocityFieldFile.replace(pos, args.outputVelocityFieldFile.size(), "-velocityField.mha");
		}
		else
		{
			args.outputVelocityFieldFile = args.outputImageFile + "-velocityField.mha";
		}
	}

    // Change the extension by -polyaffineTreeTransform.tfm
	if ( args.outputPolyaffineTreeTransformFile == "OUTPUTIMAGENAME-polyaffineTreeTransform.tfm" )
    {
		if ( pos < args.outputPolyaffineTreeTransformFile.size() )
		{
			args.outputPolyaffineTreeTransformFile = args.outputImageFile;
			args.outputPolyaffineTreeTransformFile.replace(pos, args.outputPolyaffineTreeTransformFile.size(), "-polyaffineTreeTransform.tfm");
		}
		else
		{
			args.outputPolyaffineTreeTransformFile = args.outputImageFile + "-polyaffineTreeTransform.tfm";
		}
	}
    
//    // Change the extension by -Lambda.txt
//	if ( args.outputLambdaNewFile == "OUTPUTIMAGENAME-LambdaNew.txt" )
//    {
//		if ( pos < args.outputLambdaNewFile.size() )
//		{
//			args.outputLambdaNewFile = args.outputImageFile;
//			args.outputLambdaNewFile.replace(pos, args.outputLambdaNewFile.size(), "-LambdaNew.txt");
//		}
//		else
//		{
//			args.outputLambdaNewFile = args.outputImageFile + "-LambdaNew.txt";
//		}
//	}
//
//    // Change the extension by -C.txt
//	if ( args.outputCNewFile == "OUTPUTIMAGENAME-CNew.txt" )
//    {
//		if ( pos < args.outputCNewFile.size() )
//		{
//			args.outputCNewFile = args.outputImageFile;
//			args.outputCNewFile.replace(pos, args.outputCNewFile.size(), "-CNew.txt");
//		}
//		else
//		{
//			args.outputCNewFile = args.outputImageFile + "-CNew.txt";
//		}
//	}

	// Change the extension by -correlationMatrix.txt
	if ( args.outputCorrelationMatrixFile == "OUTPUTIMAGENAME-correlationMatrix-LEVEL.txt" )
    {
		if ( pos < args.outputCorrelationMatrixFile.size() )
		{
			args.outputCorrelationMatrixFile = args.outputImageFile;
			args.outputCorrelationMatrixFile.replace(pos, args.outputCorrelationMatrixFile.size(), "-correlationMatrix.txt");
		}
		else
		{
			args.outputCorrelationMatrixFile = args.outputImageFile + "-correlationMatrix.txt";
		}
	}	
	
	// Change the extension by -metric.txt
	if ( args.outputMetricFile == "OUTPUTIMAGENAME-metric.txt" )
    {
		if ( pos < args.outputMetricFile.size() )
		{
			args.outputMetricFile = args.outputImageFile;
			args.outputMetricFile.replace(pos, args.outputMetricFile.size(), "-metric.txt");
		}
		else
		{
			args.outputMetricFile = args.outputImageFile + "-metric.txt";
		}
	}
	
	// Change the extension by -harmonic-energy.txt
	if ( args.outputHarmonicEnergyFile == "OUTPUTIMAGENAME-harmonic-energy.txt" )
    {
		if ( pos < args.outputHarmonicEnergyFile.size() )
		{
			args.outputHarmonicEnergyFile = args.outputImageFile;
			args.outputHarmonicEnergyFile.replace(pos, args.outputHarmonicEnergyFile.size(), "-harmonic-energy.txt");
		}
		else
		{
			args.outputHarmonicEnergyFile = args.outputImageFile + "-harmonic-energy.txt";
		}
	}
	
	args.numberOfLevels = command.GetValueAsInt("NumberOfLevels","intval");
	args.weightScalingValue = command.GetValueAsFloat("WeightScalingValue","floatval");
	args.maxStepLength = command.GetValueAsFloat("MaximumUpdateStepLength","floatval");
	args.varianceVelocityField = command.GetValueAsFloat("VarianceVelocityField","floatval");
    args.unionMaskDilationRadius = command.GetValueAsFloat("UnionMaskDilationRadius","floatval");
//    args.scaleLambdaPrior = command.GetValueAsFloat("ScaleLambdaPrior","floatval");
	args.updateRule = command.GetValueAsInt("UpdateRule","type");
	args.gradientType = command.GetValueAsInt("GradientType","type");
	args.numberOfBCHApproximationTerms = command.GetValueAsInt("NumberOfBCHApproximationTerms","intval");
	args.leastSquaresType = command.GetValueAsInt("LeastSquaresType","type");

}

template<typename TRegistration>
class ShowProgressObject
{
public:
	ShowProgressObject(TRegistration* o, double ratio)
    {m_Process = o; m_PreviousMetric = 0.0; m_Ratio = ratio;}
	void ShowProgress() {
		std::cout << "Iter: " << m_Process->GetElapsedIterations() << "  ";
		std::cout << "Metric: "   << m_Process->GetMetric()   << "  ";
		std::cout << "RMSChange: " << m_Process->GetRMSChange() << "  ";
		std::cout << std::endl;
		
		if(m_Process->GetMetric()==0)
			m_Process->StopRegistration();
        
		if( m_Process->GetElapsedIterations() > 2) {
			if( (m_PreviousMetric - m_Process->GetMetric())/m_Process->GetMetric() < m_Ratio )
			{ 
                m_Process->StopRegistration(); 
            }
		}
		m_PreviousMetric = m_Process->GetMetric();
	}
	
	typename TRegistration::Pointer m_Process;
	double m_PreviousMetric;
	double m_Ratio;
};

template <unsigned int Dimension>
void LogDomainDemonsRegistrationFunction( arguments args )
{	
    itkProbesCreate();
    itkProbesStart( "Registration" );
    
	// Images we use
	typename ImageType::Pointer fixedImage = 0;
	typename ImageType::Pointer movingImage = 0;
	typename LabelImageType::Pointer fixedMaskImage = 0;
	typename LabelImageType::Pointer movingMaskImage = 0;
    typename LabelImageType::Pointer manualRegionsImage = 0;
	typename FieldType::Pointer inputVelField = 0;
	
	
	// Set up the file readers
	typedef itk::ImageFileReader< ImageType >         FixedImageReaderType;
	typedef itk::ImageFileReader< ImageType >         MovingImageReaderType;
	typedef itk::ImageFileReader< LabelImageType > FixedMaskReaderType;
	typedef itk::ImageFileReader< LabelImageType > MovingMaskReaderType;
	typedef itk::ImageFileReader< FieldType > VelocityFieldReaderType;
	typedef itk::TransformFileReader                  TransformReaderType;
	

    typename FixedImageReaderType::Pointer fixedImageReader = FixedImageReaderType::New();
    typename MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
    
    fixedImageReader->SetFileName( args.fixedImageFile.c_str() );
    movingImageReader->SetFileName( args.movingImageFile.c_str() );
    
    typename FixedMaskReaderType::Pointer fixedMaskImageReader = FixedMaskReaderType::New();
    typename MovingMaskReaderType::Pointer movingMaskImageReader = MovingMaskReaderType::New();
    
    fixedMaskImageReader->SetFileName(  args.fixedmaskImageFile.c_str() );
    movingMaskImageReader->SetFileName(  args.movingmaskImageFile.c_str() );
    
    // Update the reader
    try
    {
        fixedImageReader->Update();
        fixedImage = fixedImageReader->GetOutput();
        fixedImage->DisconnectPipeline();
        
        movingImageReader->Update();
        movingImage = movingImageReader->GetOutput();
        movingImage->DisconnectPipeline();
        
        fixedMaskImageReader->Update();
        fixedMaskImage = fixedMaskImageReader->GetOutput();
        fixedMaskImage->DisconnectPipeline();

        movingMaskImageReader->Update();
        movingMaskImage = movingMaskImageReader->GetOutput();
        movingMaskImage->DisconnectPipeline();
        
        typedef itk::HistogramMatchingImageFilter< ImageType,ImageType > HistogramMatchingFilterType;
        HistogramMatchingFilterType::Pointer histogramMatchingFilter = HistogramMatchingFilterType::New();
        histogramMatchingFilter->SetReferenceImage( fixedImage );
        histogramMatchingFilter->SetSourceImage( movingImage );
        histogramMatchingFilter->SetNumberOfHistogramLevels( 50 );
        histogramMatchingFilter->SetNumberOfMatchPoints( 8 );
        histogramMatchingFilter->ThresholdAtMeanIntensityOn();
        histogramMatchingFilter->Update();
        
        movingImage = histogramMatchingFilter->GetOutput();
        movingImage->DisconnectPipeline();
        
        typedef itk::ImageFileWriter< ImageType > ImageWriterType;
        ImageWriterType::Pointer imageWriter = ImageWriterType::New();
        imageWriter->SetInput(movingImage);
        imageWriter->SetFileName("HistMatchMovingImage.mhd");
        imageWriter->Update();
                    
//            // masking for LONI images
//            typedef itk::MaskImageFilter< ImageType,LabelImageType,LabelImageType > MaskImageFilterType;
//            MaskImageFilterType::Pointer maskFixedFilter = MaskImageFilterType::New();
//            maskFixedFilter->SetInput1(fixedMaskImageReader->GetOutput());
//            maskFixedFilter->SetInput2(fixedImage);
//            maskFixedFilter->Update();
//			fixedMaskImage = maskFixedFilter->GetOutput();
//			fixedMaskImage->DisconnectPipeline();
//
//            MaskImageFilterType::Pointer maskMovingFilter = MaskImageFilterType::New();
//            maskMovingFilter->SetInput1(movingMaskImageReader->GetOutput());
//            maskMovingFilter->SetInput2(movingImage);
//            maskMovingFilter->Update();
//			movingMaskImage = maskMovingFilter->GetOutput();
//			movingMaskImage->DisconnectPipeline();
//            
//            typedef itk::ImageFileWriter< LabelImageType > LabelWriterType;
//            LabelWriterType::Pointer labelWriter = LabelWriterType::New();
//            labelWriter->SetInput(fixedMaskImage);
//            labelWriter->SetFileName("FixedMask.mhd");
//            labelWriter->Update();
//            
//            labelWriter->SetInput(movingMaskImage);
//            labelWriter->SetFileName("MovingMask.mhd");
//            labelWriter->Update();
                    
        // multiresolution
//            const unsigned int shrinkFactor = 4;
//            typedef itk::MultiResolutionPyramidImageFilter< ImageType,ImageType > MultiResolutionFilter;
//            MultiResolutionFilter::Pointer multiResolutionFixedFilter = MultiResolutionFilter::New();
//            multiResolutionFixedFilter->SetInput(fixedImage);
//            multiResolutionFixedFilter->SetNumberOfLevels(1);
//            multiResolutionFixedFilter->SetStartingShrinkFactors(shrinkFactor);
//            multiResolutionFixedFilter->Update();
//            fixedImage = multiResolutionFixedFilter->GetOutput(0);
//            fixedImage->DisconnectPipeline();
//            
//            MultiResolutionFilter::Pointer multiResolutionMovingFilter = MultiResolutionFilter::New();
//            multiResolutionMovingFilter->SetInput(movingImage);
//            multiResolutionMovingFilter->SetNumberOfLevels(1);
//            multiResolutionMovingFilter->SetStartingShrinkFactors(shrinkFactor);
//            multiResolutionMovingFilter->Update();
//            movingImage = multiResolutionMovingFilter->GetOutput(0);
//            movingImage->DisconnectPipeline();
//            
//            typedef itk::MultiResolutionPyramidImageFilter< LabelImageType,LabelImageType > MultiResolutionMaskFilter;
//            MultiResolutionMaskFilter::Pointer multiResolutionFixedMaskFilter = MultiResolutionMaskFilter::New();
//            multiResolutionFixedMaskFilter->SetInput(fixedMaskImage);
//            multiResolutionFixedMaskFilter->SetNumberOfLevels(1);
//            multiResolutionFixedMaskFilter->SetStartingShrinkFactors(shrinkFactor);
//            multiResolutionFixedMaskFilter->Update();
//            fixedMaskImage = multiResolutionFixedMaskFilter->GetOutput(0);
//            fixedMaskImage->DisconnectPipeline();
//            
//            MultiResolutionMaskFilter::Pointer multiResolutionMovingMaskFilter = MultiResolutionMaskFilter::New();
//            multiResolutionMovingMaskFilter->SetInput(movingMaskImage);
//            multiResolutionMovingMaskFilter->SetNumberOfLevels(1);
//            multiResolutionMovingMaskFilter->SetStartingShrinkFactors(shrinkFactor);
//            multiResolutionMovingMaskFilter->Update();
//            movingMaskImage = multiResolutionMovingMaskFilter->GetOutput(0);
//            movingMaskImage->DisconnectPipeline();
        
        // bias correction
//            typedef itk::MRIBiasFieldCorrectionFilter< ImageType,ImageType,LabelImageType > BiasCorrectionFilter;
//            BiasCorrectionFilter::Pointer biasCorrectionFixedFilter = BiasCorrectionFilter::New();
//            biasCorrectionFixedFilter->SetInput(fixedImage);
//            biasCorrectionFixedFilter->Update();
//            fixedImage = biasCorrectionFixedFilter->GetOutput();
//            fixedImage->DisconnectPipeline();
//            
//            BiasCorrectionFilter::Pointer biasCorrectionMovingFilter = BiasCorrectionFilter::New();
//            biasCorrectionMovingFilter->SetInput(movingImage);
//            biasCorrectionMovingFilter->Update();
//            movingImage = biasCorrectionMovingFilter->GetOutput();
//            movingImage->DisconnectPipeline();
        
        // smoothing
        //typedef itk::RecursiveGaussianImageFilter< ImageType,ImageType > SmoothingFilterType;
//            typedef itk::DiscreteGaussianImageFilter< ImageType,ImageType > SmoothingFilterType;
//            SmoothingFilterType::Pointer smoothingFixedFilter = SmoothingFilterType::New();
//            smoothingFixedFilter->SetInput(fixedImage);
//            smoothingFixedFilter->SetVariance(20);
//            smoothingFixedFilter->Update();
//            fixedImage = smoothingFixedFilter->GetOutput();
//            fixedImage->DisconnectPipeline();
//            
//            SmoothingFilterType::Pointer smoothingMovingFilter = SmoothingFilterType::New();
//            smoothingMovingFilter->SetInput(movingImage);
//            smoothingMovingFilter->SetVariance(20);
//            smoothingMovingFilter->Update();
//            movingImage = smoothingMovingFilter->GetOutput();
//            movingImage->DisconnectPipeline();

    }
    catch( itk::ExceptionObject& err )
    {
        std::cout << "Could not read one of the input images." << std::endl;
        std::cout << err << std::endl;
        exit( EXIT_FAILURE );
    }
            
    if ( ! args.inputFieldFile.empty() )
    {
        // Set up the file readers
        typename VelocityFieldReaderType::Pointer fieldReader = VelocityFieldReaderType::New();
        fieldReader->SetFileName(  args.inputFieldFile.c_str() );
        
        // Update the reader
        try
        {
            fieldReader->Update();
        }
        catch( itk::ExceptionObject& err )
        {
            std::cout << "Could not read the input field." << std::endl;
            std::cout << err << std::endl;
            exit( EXIT_FAILURE );
        }
        
        inputVelField = fieldReader->GetOutput();
        inputVelField->DisconnectPipeline();
    }
    
    if ( ! args.manualRegionsLabelFile.empty() ) {
        // Set up the file readers
        typename FixedMaskReaderType::Pointer manualRegionsReader = FixedMaskReaderType::New();
        manualRegionsReader->SetFileName(  args.manualRegionsLabelFile.c_str() );
        
        // Update the reader
        try
        {
            manualRegionsReader->Update();
        }
        catch( itk::ExceptionObject& err )
        {
            std::cout << "Could not read the manual regions label image." << std::endl;
            std::cout << err << std::endl;
            exit( EXIT_FAILURE );
        }
        
        manualRegionsImage = manualRegionsReader->GetOutput();
        manualRegionsImage->DisconnectPipeline();
        
    }

    // Set up the demons filter output deformation field
    typename FieldType::Pointer defField = 0;
    typename FieldType::Pointer invDefField = 0;
    
    typedef itk::ManualPolyaffineStationaryVelocityFieldTransform<> ManaualPolyaffineTreeTransformType;
    typedef itk::PolyaffineStationaryVelocityFieldTransform<> PolyaffineTreeTransformType;
    PolyaffineTreeTransformType::Pointer polyaffineTreeTransform = 0;
    
    std::vector<double> metricVec;
    std::vector<double> harmonicEnergyVec;

    WeightImageType::Pointer featureImage = 0;
    const unsigned int zeroLevel = 0;
    
    if(!manualRegionsImage) {
        
        // gaussian mixture tree
        //GaussianMixtureModelTree gaussianTree;
        
        // structure tensor
        /*double derivationSigma = 2;
         double integrationSigma = 2;
         double meanThresholdSigma = 0.95;
         featureImage = gaussianTree.computeHeterogeneousFeatureImage(fixedImage, fixedMaskImage, derivationSigma, integrationSigma, meanThresholdSigma);*/
        
        if(args.leastSquaresType == 0) {
            // log(gradient + 1)
            
            // gradient
            typedef itk::GradientMagnitudeImageFilter< ImageType,WeightImageType > GradientMagnitudeImageFilterType;
            GradientMagnitudeImageFilterType::Pointer gradMagFilter = GradientMagnitudeImageFilterType::New();
            gradMagFilter->SetInput(fixedImage);
            
            typedef itk::SquareImageFilter< WeightImageType,WeightImageType > SquareImageFilterType;
            typename SquareImageFilterType::Pointer squareImageFilter = SquareImageFilterType::New();
            squareImageFilter->SetInput(gradMagFilter->GetOutput());
            
            typedef itk::AddConstantToImageFilter< WeightImageType,float,WeightImageType > AddConstantFilterType;
            typename AddConstantFilterType::Pointer addConstantFilter = AddConstantFilterType::New();
            addConstantFilter->SetInput(squareImageFilter->GetOutput());
            addConstantFilter->SetConstant(1);
            
            typedef itk::LogImageFilter< WeightImageType,WeightImageType > LogImageFilterType;
            typename LogImageFilterType::Pointer logFilter = LogImageFilterType::New();
            logFilter->SetInput(addConstantFilter->GetOutput());
            
            typedef itk::MaskImageFilter< WeightImageType,LabelImageType > MaskImageFilterType;
            MaskImageFilterType::Pointer maskFeatureImageFilter = MaskImageFilterType::New();
            maskFeatureImageFilter->SetInput1(logFilter->GetOutput());
            maskFeatureImageFilter->SetInput2(fixedMaskImage);
            maskFeatureImageFilter->Update();
            
            featureImage = maskFeatureImageFilter->GetOutput();
            featureImage->DisconnectPipeline();
                        
            /*gradMagFilter->Update();
             typedef itk::StatisticsImageFilter< WeightImageType > StatisticsImageFilterType;
             StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
             statisticsImageFilter->SetInput(gradMagFilter->GetOutput());
             statisticsImageFilter->Update();
             StatisticsImageFilterType::PixelType maxGradMag = statisticsImageFilter->GetMaximum();
             std::cout << "max = " << maxGradMag << std::endl;
             
             typedef itk::ImageRegionIterator< WeightImageType > ImageIterator;
             ImageIterator iter(gradMagFilter->GetOutput(), gradMagFilter->GetOutput()->GetRequestedRegion());
             for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
             double gradValue = iter.Get();
             //iter.Set(vcl_exp(-(maxGradMag*maxGradMag)/(gradValue*gradValue)));
             iter.Set((gradValue*gradValue)/(maxGradMag*maxGradMag));
             }*/            
            
        }
        else if(args.leastSquaresType == 1) {
            // contour
            
            typedef itk::MinimumMaximumImageFilter< LabelImageType > MinimumMaximumImageFilterType;
            MinimumMaximumImageFilterType::Pointer minMaxFilter = MinimumMaximumImageFilterType::New();
            minMaxFilter->SetInput(fixedMaskImage);
            minMaxFilter->Update();
            MinimumMaximumImageFilterType::PixelType minLabel = minMaxFilter->GetMinimum();
            std::cout << "minLabel = " << minLabel << std::endl;

            typedef itk::LabelContourImageFilter< LabelImageType,WeightImageType > ContourFilterType;
            ContourFilterType::Pointer contourFilter = ContourFilterType::New();
            contourFilter->SetBackgroundValue(minLabel+0.1);
            contourFilter->SetInput(fixedMaskImage);
            contourFilter->Update();
             
            featureImage = contourFilter->GetOutput();
            featureImage->DisconnectPipeline();

        }
        else if(args.leastSquaresType == 2) {
            typedef itk::CastImageFilter< LabelImageType,WeightImageType > CastImageFilterType;
            CastImageFilterType::Pointer castImageFilter = CastImageFilterType::New();
            castImageFilter->SetInput(fixedMaskImage);
            featureImage = castImageFilter->GetOutput();
            featureImage->Update();
        }
        
        // write feature image for debugging
        typedef itk::ImageFileWriter< WeightImageType > FeatureWriterType;
        FeatureWriterType::Pointer featureWriter = FeatureWriterType::New();
        featureWriter->SetInput(featureImage);
        featureWriter->SetFileName("AnchorImage.mhd");
        featureWriter->Update();
                
        // hierarchical regions
        polyaffineTreeTransform = PolyaffineTreeTransformType::New();
        
        PolyaffineTreeTransformType::GaussianTreeNodePointerType rootNode = PolyaffineTreeTransformType::GaussianTreeNodeType::New();
        WeightedSVDTree weightedTree;
        weightedTree.weightTreeEstimation(zeroLevel, args.numberOfLevels-1, args.weightScalingValue, featureImage, polyaffineTreeTransform, rootNode);
        polyaffineTreeTransform->SetWeights(rootNode);
        
    }    
    
    else {
            
        // manual regions
        polyaffineTreeTransform = ManaualPolyaffineTreeTransformType::New();
        
        featureImage = WeightImageType::New();
        featureImage->CopyInformation(manualRegionsImage);
        featureImage->SetRegions(manualRegionsImage->GetRequestedRegion());
        featureImage->Allocate();
                
        // root node is zero (if we don't want to do a global affine prior to the registration)
        typedef itk::ImageRegionConstIterator< LabelImageType > MaskImageIterator;
        MaskImageIterator labelIter(manualRegionsImage, manualRegionsImage->GetRequestedRegion());
        
        typedef itk::ImageRegionIterator< WeightImageType > FeatureImageIterator;
        FeatureImageIterator featureIter(featureImage, featureImage->GetRequestedRegion());
        
        featureImage->FillBuffer(0);
        for(featureIter.GoToBegin(), labelIter.GoToBegin(); !featureIter.IsAtEnd(); ++featureIter, ++labelIter) {
            if(labelIter.Get()) {
                featureIter.Set(1);
            }
        }
        
        WeightedSVDTree weightedTree;
        vnl_vector_fixed< PolyaffineTreeTransformType::MatrixElementType,Dimension > meanPoint;
        vnl_matrix_fixed< PolyaffineTreeTransformType::MatrixElementType,Dimension,Dimension > covarianceMatrix;
        weightedTree.computeGaussianWeight(featureImage, args.weightScalingValue, meanPoint, covarianceMatrix);
        std::cout << "meanPoint = " << meanPoint << std::endl;
        std::cout << "covarianceMatrix = " << covarianceMatrix << std::endl;
        
//        vnl_vector_fixed< PolyaffineTreeTransformType::MatrixElementType,Dimension > zeroMeanPoint;
//        zeroMeanPoint.fill(0);
//        vnl_matrix_fixed< PolyaffineTreeTransformType::MatrixElementType,Dimension,Dimension > zeroCovarianceMatrix;
//        zeroCovarianceMatrix.fill(0);
//        PolyaffineTreeTransformType::NodeValueType zeroParams = polyaffineTreeTransform->GaussianToParam(zeroMeanPoint, zeroCovarianceMatrix);
        PolyaffineTreeTransformType::NodeValueType zeroParams = polyaffineTreeTransform->GaussianToParam(meanPoint, covarianceMatrix);
        PolyaffineTreeTransformType::GaussianTreeNodePointerType manualRootNode = PolyaffineTreeTransformType::GaussianTreeNodeType::New();
        manualRootNode->Set(zeroParams);
        
        typedef itk::MinimumMaximumImageFilter< LabelImageType > MinimumMaximumImageFilterType;
        MinimumMaximumImageFilterType::Pointer minMaxFilter = MinimumMaximumImageFilterType::New();
        minMaxFilter->SetInput(manualRegionsImage);
        minMaxFilter->Update();
        MinimumMaximumImageFilterType::PixelType maximumLabel = minMaxFilter->GetMaximum();
        std::cout << "maximumLabel = " << maximumLabel << std::endl;
        for(unsigned int i = 1; i <= maximumLabel; ++i) {
            
            MinimumMaximumImageFilterType::PixelType currentLabel = i;            
            featureImage->FillBuffer(0);
            
            MaskImageIterator labelIter(manualRegionsImage, manualRegionsImage->GetRequestedRegion());
            FeatureImageIterator featureIter(featureImage, featureImage->GetRequestedRegion());
            for(featureIter.GoToBegin(), labelIter.GoToBegin(); !featureIter.IsAtEnd(); ++featureIter, ++labelIter) {
                if(labelIter.Get() == currentLabel) {
                    featureIter.Set(1);
                }
            }

            weightedTree.computeGaussianWeight(featureImage, args.weightScalingValue, meanPoint, covarianceMatrix);
            std::cout << "meanPoint = " << meanPoint << std::endl;
            std::cout << "covarianceMatrix = " << covarianceMatrix << std::endl;            
            
            // first level are the regions defined by the label image
            PolyaffineTreeTransformType::GaussianTreeNodeType::Pointer child = PolyaffineTreeTransformType::GaussianTreeNodeType::New();
            manualRootNode->AddChild(child);
            PolyaffineTreeTransformType::NodeValueType estimatedParams = polyaffineTreeTransform->GaussianToParam(meanPoint, covarianceMatrix);
            child->Set(estimatedParams);
            
        }

        polyaffineTreeTransform->SetWeights(manualRootNode);
    }
  
    /*
    vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > MeanSetAllLevels;
    vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > CovSetAllLevels;
    //const unsigned int noOfLevels = args.numberOfLevels;
    //const unsigned int noOfLevelsPrevious = noOfLevels-1;
    //const unsigned int treeDepthPreviousLevel = noOfLevelsPrevious-1;
    const unsigned int treeDepth = args.numberOfLevels-1;
    //unsigned int noOfNodesPreviousLevels = std::pow((double)2,(double)(treeDepthPreviousLevel+1))-1;
    //std::cout << "noOfNodesPreviousLevels = " << noOfNodesPreviousLevels << std::endl;
    unsigned int noOfNodes = std::pow((double)2,(double)(treeDepth+1))-1;
    //std::cout << "noOfNodes = " << noOfNodes << std::endl;

    // covariance for all nodes plus one (prior)    
    if ( ! args.priorTransformDirectory.empty() )
    {
        // reading transforms
        itk::RegularExpressionSeriesFileNames::Pointer directory = itk::RegularExpressionSeriesFileNames::New();
        directory->SetDirectory(args.priorTransformDirectory.c_str());
        directory->SetRegularExpression("_polyaffine_transform.tfm");
        std::vector< std::string > fileNames = directory->GetFileNames();
        
        //for(unsigned int i = 0; i < fileNames.size(); ++i)
        //    std::cout << i << " " << fileNames[i] << std::endl;
                        
        std::vector< PolyaffineTreeTransformType::Pointer > transformSet;
        for(unsigned int i = 0; i < fileNames.size(); ++i) {
            
            typedef itk::TransformFileReader::TransformListType* TransformListType;
            itk::TransformFactory< PolyaffineTreeTransformType >::RegisterTransform();
            itk::TransformFileReader::Pointer transformReader = itk::TransformFileReader::New();
            transformReader->SetFileName(fileNames[i]);
            transformReader->Update();
            PolyaffineTreeTransformType::Pointer transform = static_cast< PolyaffineTreeTransformType* >(transformReader->GetTransformList()->begin()->GetPointer());
            transformSet.push_back(transform);
            
        }

        // filling the data matrix
        //vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > dataMatrix(Dimension*(Dimension+1)*noOfNodesPreviousLevels,transformSet.size());
        vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > dataMatrix(Dimension*(Dimension+1)*noOfNodes,transformSet.size());
        
        unsigned int rowIndex = 0;
        //for(unsigned int level = 0; level < noOfLevelsPrevious; ++level) {
        for(unsigned int level = 0; level < args.numberOfLevels; ++level) {

            std::vector< PolyaffineTreeTransformType::GaussianTreeNodePointerType > nodeSet;
            transformSet[0]->GetTransformationsAtLevel(level, nodeSet);
            
            for(unsigned int j = 0; j < nodeSet.size(); ++j, rowIndex += Dimension*(Dimension+1)) {
                                                
                for(unsigned int i = 0; i < transformSet.size(); ++i) {
                    std::vector< PolyaffineTreeTransformType::GaussianTreeNodePointerType > currentNodeSet;
                    transformSet[i]->GetTransformationsAtLevel(level, currentNodeSet);
                    
                    vnl_matrix_fixed< PolyaffineTreeTransformType::MatrixElementType,Dimension+1,Dimension+1 > M;
                    PolyaffineTreeTransformType::ParametersType param;
                    param = currentNodeSet[j]->Get();
                    transformSet[0]->ParamToLogTransform(param, M);
                    
                    vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > MnonNull = M.extract(Dimension, Dimension+1);
                    
                    // vectorize M
                    vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > Mvect(MnonNull.rows()*MnonNull.cols(),1);
                    for(unsigned int n = 0; n < MnonNull.cols(); ++n) {
                        for(unsigned int m = 0; m < MnonNull.rows(); ++m) {
                            Mvect(n*MnonNull.rows()+m,0) = MnonNull(m,n);
                        }
                    }
                    
                    for(unsigned int m = 0; m < Mvect.size(); ++m)
                        dataMatrix(rowIndex+m,i) = Mvect(m,0);                                        
                    
                }
            }
        }
        //std::cout << "dataMatrix = \n" << dataMatrix << std::endl;
        std::ofstream fileDataMatrix;
        fileDataMatrix.open ("DataMatrix.txt");
        fileDataMatrix << dataMatrix << std::endl;
        fileDataMatrix.close();
        
        // compute the mean
        vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > Mmean(dataMatrix.rows(), 1);
        Mmean.fill(0.0);
        for(unsigned int j = 0; j < dataMatrix.cols(); ++j)
            Mmean += (1.0/dataMatrix.cols())*dataMatrix.get_n_columns(j,1);
        //std::cout << "Mmean = \n" << Mmean << std::endl;

        // center the data matrix
        for(unsigned int j = 0; j < dataMatrix.cols(); ++j)
            dataMatrix.set_columns(j,dataMatrix.get_n_columns(j,1) - Mmean);
        //std::cout << "center data matrix = \n" << dataMatrix << std::endl;
        
        vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > Mtest(dataMatrix.rows(), 1);
        Mtest.fill(0.0);
        for(unsigned int j = 0; j < dataMatrix.cols(); ++j)
            Mtest += dataMatrix.get_n_columns(j,1);
        //std::cout << "Mtest = \n" << Mtest << std::endl;
        
        // compute sample covariance matrix
        vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > Cov = (1.0/dataMatrix.cols())*dataMatrix*dataMatrix.transpose();
        //std::cout << "Cov = \n" << Cov << std::endl;
        std::ofstream fileCovMatrix;
        fileCovMatrix.open ("Cov.txt");
        fileCovMatrix << Cov << std::endl;
        fileCovMatrix.close();

        // block pca
        std::vector< vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > > CovRegions(Cov.rows());
        const unsigned int noOfParameters = 12;
        const unsigned int noOfTransformations = Cov.rows()/noOfParameters;
        vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > C(noOfParameters,noOfParameters);
        C.fill(0.0);
        std::cout << "no of transforms = " << noOfTransformations << std::endl;
        for(unsigned int i = 0; i < noOfTransformations; ++i) {
            unsigned int currentRegionStartIndex = i*noOfParameters;
            CovRegions[i] = Cov.extract(noOfParameters,noOfParameters,currentRegionStartIndex,currentRegionStartIndex);
            C += CovRegions[i];
        }
        C *= (1.0/noOfTransformations);

        std::ofstream cMatrixFile;
        cMatrixFile.open ("CMatrix.txt");
        cMatrixFile << C << std::endl;
        cMatrixFile.close();
        
        vnl_vector< double > lambdas(noOfTransformations);
        for(unsigned int i = 0; i < noOfTransformations; ++i)
            lambdas[i] = vnl_trace(CovRegions[i]*C)/vnl_trace(vnl_transpose(C)*C);

        std::cout << "lambdas = \n" << lambdas << std::endl;
        std::ofstream lambdasFile;
        lambdasFile.open ("lambdas.txt");
        lambdasFile << lambdas << std::endl;
        lambdasFile.close();
        
        for(unsigned int i = 0; i < noOfTransformations; ++i)
            lambdas[i] *= args.scaleLambdaPrior;
        std::cout << "level scaled lambdas = \n" << lambdas << std::endl;

        const PolyaffineTreeTransformType::MatrixElementType lamdbaMax = 1e-8;
        vnl_vector< double > pseudoInverseLambdas(noOfTransformations);
        pseudoInverseLambdas.fill(0.0);
        for(unsigned int i = 0; i < noOfTransformations; ++i) {
            if(lambdas[i] > lamdbaMax)
                pseudoInverseLambdas[i] = 1.0/lambdas[i];
        }
        std::cout << "pseudoInverseLambdas = \n" << pseudoInverseLambdas << std::endl;
        
        vnl_cholesky CInverse(C);
        std::cout << "CInverse.inverse() = \n" << CInverse.inverse() << std::endl;
        vnl_diag_matrix< PolyaffineTreeTransformType::MatrixElementType > diag(pseudoInverseLambdas);
        //std::cout << "diag.asMatrix() = \n" << diag.asMatrix() << std::endl;
        vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > BlockWisePCA = itk::PolyaffineLogDomainDemonsRegistrationFilter < ImageType, ImageType, FieldType, LabelImageType, LabelImageType, WeightImageType >::KroneckerProduct(diag.asMatrix(), CInverse.inverse());

        std::ofstream blockWisePCAFile;
        blockWisePCAFile.open ("BlockWisePCA.txt");
        blockWisePCAFile << BlockWisePCA << std::endl;
        blockWisePCAFile.close();
        
        //exit(EXIT_SUCCESS);
        
//        ifstream fileCovInverse("CovInverse.txt");
//        vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > CovInverse;
//        CovInverse.read_ascii(fileCovInverse);
//        fileCovInverse.close();
        
        MeanSetAllLevels = Mmean;
        //CovSetAllLevels = Cov;
        //CovSetAllLevels = CovInverse;
        CovSetAllLevels = BlockWisePCA;
                
//        // paste into bigger matrix with one additional level with zeros (no prior for that level available yet)
//        MeanSetAllLevels.set_size(Dimension*(Dimension+1)*noOfNodes, 1);
//        MeanSetAllLevels.fill(0.0);
//        for(unsigned int i = 0; i < Mmean.rows(); ++i) {
//            for(unsigned int j = 0; j < Mmean.cols(); ++j) {
//                 MeanSetAllLevels(i,j) = Mmean(i,j);
//            }
//        }
//        //std::cout << "MeanSetAllLevels = \n" << MeanSetAllLevels << std::endl;
//        
//        CovSetAllLevels.set_size(Dimension*(Dimension+1)*noOfNodes, Dimension*(Dimension+1)*noOfNodes);
//        CovSetAllLevels.fill(0.0);
//        for(unsigned int i = 0; i < Cov.rows(); ++i) {
//            for(unsigned int j = 0; j < Cov.cols(); ++j) {
//                CovSetAllLevels(i,j) = Cov(i,j);
//            }
//        }
        
        //std::cout << "CovSetAllLevels = \n" << CovSetAllLevels << std::endl;
        //std::cout << "CovSetAllLevels rank = \n" << inverse.rank() << std::endl;
        //std::cout << "CovSetAllLevels pseudo inverse = \n" << inverse.pinverse() << std::endl;
        
//        vnl_matrix_inverse< PolyaffineTreeTransformType::MatrixElementType > CovInverse(Cov);
//         PolyaffineTreeTransformType::MatrixElementType mahalanobisSquared = (Mmean.transpose()*CovInverse*Mmean)(0,0);
//         PolyaffineTreeTransformType::MatrixElementType p = Mmean.rows()*Mmean.cols();
//         PolyaffineTreeTransformType::MatrixElementType n = transformSet.size();
//         PolyaffineTreeTransformType::MatrixElementType t2 = n*mahalanobisSquared;
//         std::cout << "Mahalanobis distance squared of the mean = " << mahalanobisSquared << std::endl;
//         std::cout << "Standard deviations from zero = " << std::sqrt(mahalanobisSquared) << std::endl;
//         std::cout << "Hotelling's T-squared statistic of the mean, t^2 = " << t2 << std::endl;
//         PolyaffineTreeTransformType::MatrixElementType Fdist = ((n-p)/(p*(n-1)))*t2;
//         std::cout << "F distribution F_{" << p << "," << n << "-" << p << "} of the mean = " << Fdist << std::endl;
        
    }
    */

    vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > SigmaInverse;
    if ( ! args.priorSigmaInverse.empty() )
    {
        std::ifstream fileSigmaInverse(args.priorSigmaInverse.c_str());
        SigmaInverse.read_ascii(fileSigmaInverse);
        fileSigmaInverse.close();
    }
    
    vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > mui;
    if ( ! args.priorMui.empty() )
    {
        std::ifstream fileMui(args.priorMui.c_str());
        mui.read_ascii(fileMui);
        fileMui.close();
    }    
    //std::cout << "mui = " << mui << std::endl;
    
//    vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > LambdaZeroMean;
//    // covariance for all nodes plus one (prior)    
//    if ( ! args.priorTransformDirectoryLambda.empty() )
//    {        
//        // reading transforms
//        itk::RegularExpressionSeriesFileNames::Pointer directory = itk::RegularExpressionSeriesFileNames::New();
//        directory->SetDirectory(args.priorTransformDirectoryLambda.c_str());
//        directory->SetRegularExpression("LambdaNew.txt");
//        std::vector< std::string > fileNamesLambdaZero = directory->GetFileNames();
//        for(unsigned int i = 0; i < fileNamesLambdaZero.size(); ++i) {
//            std::cout << fileNamesLambdaZero[i] << std::endl;
//            
//            std::ifstream fileLambdaZero(fileNamesLambdaZero[i].c_str());
//            vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > LambdaZero;
//            LambdaZero.read_ascii(fileLambdaZero);
//            fileLambdaZero.close();
//            
//            if(LambdaZeroMean.empty())
//                LambdaZeroMean = LambdaZero;
//            else
//                LambdaZeroMean += LambdaZero;
//            
//        }
//        LambdaZeroMean *= 1.0/fileNamesLambdaZero.size();
//        std::cout << "LambdaZeroMean = \n" << LambdaZeroMean << std::endl;
//        
//    }
//
//    vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > CZeroMean;
//    if ( ! args.priorTransformDirectoryC.empty() )
//    {   
//        itk::RegularExpressionSeriesFileNames::Pointer directory = itk::RegularExpressionSeriesFileNames::New();
//        directory->SetDirectory(args.priorTransformDirectoryC.c_str());
//        directory->SetRegularExpression("CNew.txt");
//        std::vector< std::string > fileNamesCZero = directory->GetFileNames();
//        for(unsigned int i = 0; i < fileNamesCZero.size(); ++i) {
//            std::cout << fileNamesCZero[i] << std::endl;
//            
//            std::ifstream fileCZero(fileNamesCZero[i].c_str());
//            vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > CZero;
//            CZero.read_ascii(fileCZero);
//            fileCZero.close();
//            
//            if(CZeroMean.empty())
//                CZeroMean = CZero;
//            else
//                CZeroMean += CZero;
//        }
//        CZeroMean *= 1.0/fileNamesCZero.size();
//        std::cout << "CZeroMean = \n" << CZeroMean << std::endl;
//        
//    }
    
    // Set up the demons filter
    typedef typename itk::PolyaffineLogDomainDeformableRegistrationFilter< ImageType, ImageType, FieldType, LabelImageType, LabelImageType, WeightImageType> BaseRegistrationFilterType;
    typename BaseRegistrationFilterType::Pointer filter;
    
    switch (args.updateRule)
    {
        case 0:
        {
            // exp(v) <- exp(v) o exp(u) (log-domain demons)
            typedef typename itk::PolyaffineLogDomainDemonsRegistrationFilter < ImageType, ImageType, FieldType, LabelImageType, LabelImageType, WeightImageType > ActualRegistrationFilterType;
            typedef typename ActualRegistrationFilterType::GradientType GradientType;
            
            typename ActualRegistrationFilterType::Pointer actualfilter = ActualRegistrationFilterType::New();
            
            actualfilter->SetMaximumUpdateStepLength( args.maxStepLength );
            actualfilter->SetUseGradientType(static_cast<GradientType>(args.gradientType) );
            actualfilter->SetNumberOfBCHApproximationTerms(args.numberOfBCHApproximationTerms);
            actualfilter->SetVarianceVelocityField(args.varianceVelocityField);
            actualfilter->SetUnionMaskDilationRadius(args.unionMaskDilationRadius);
            
            typedef ShowProgressObject<ActualRegistrationFilterType> RegistrationProgressType;
            RegistrationProgressType* progressWatch = new RegistrationProgressType(actualfilter, 0.01);
            typedef typename itk::SimpleMemberCommand<RegistrationProgressType> CommandType;
            typename CommandType::Pointer command = CommandType::New();
            command->SetCallbackFunction(progressWatch, &RegistrationProgressType::ShowProgress);
            actualfilter->AddObserver( itk::ProgressEvent(), command);
            
            filter = actualfilter;
            
            break;
        }
        case 1:
        {
            // exp(v) <- Symmetrized( exp(v) o exp(u) ) (symmetriclog-domain demons)
            typedef typename itk::PolyaffineSymmetricLogDomainDemonsRegistrationFilter< ImageType, ImageType, FieldType, LabelImageType, LabelImageType, WeightImageType > ActualRegistrationFilterType;
            typedef typename ActualRegistrationFilterType::GradientType GradientType;
            
            typename ActualRegistrationFilterType::Pointer actualfilter = ActualRegistrationFilterType::New();
            
            actualfilter->SetMaximumUpdateStepLength( args.maxStepLength );
            actualfilter->SetUseGradientType(static_cast<GradientType>(args.gradientType) );
            actualfilter->SetNumberOfBCHApproximationTerms(args.numberOfBCHApproximationTerms);
            actualfilter->SetVarianceVelocityField(args.varianceVelocityField);
            actualfilter->SetUnionMaskDilationRadius(args.unionMaskDilationRadius);
            
            typedef ShowProgressObject<ActualRegistrationFilterType> RegistrationProgressType;
            RegistrationProgressType* progressWatch = new RegistrationProgressType(actualfilter, 0.01);
            typedef typename itk::SimpleMemberCommand<RegistrationProgressType> CommandType;
            typename CommandType::Pointer command = CommandType::New();
            command->SetCallbackFunction(progressWatch, &RegistrationProgressType::ShowProgress);
            actualfilter->AddObserver( itk::ProgressEvent(), command);
            
            filter = actualfilter;
            
            break;
        }
        default:
        {
            std::cout << "Unsupported update rule." << std::endl;
            exit( EXIT_FAILURE );
        }
    }

    std::vector< vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > > LambdaNewAll;
    std::vector< vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > > CNewAll;
    std::vector< vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > > SigmaNewInverseAll;
    
			// multiple levels
			for(unsigned int level = zeroLevel; level < args.numberOfLevels; ++level) {
                std::ostringstream currLevelStr;
                currLevelStr << "Number of levels " << args.numberOfLevels;
                itkProbesStart(currLevelStr.str().c_str());
				
                // obtree
				/*typedef itk::ImageToVTKImageFilter< LabelImageType > VTKExporterType;
				VTKExporterType::Pointer exporterFixed = VTKExporterType::New();
				exporterFixed->SetInput(fixedMaskImage);
				vtkSmartPointer< vtkMarchingCubes > marchingFixed = vtkSmartPointer< vtkMarchingCubes >::New();
				marchingFixed->SetInput(exporterFixed->GetOutput());
				marchingFixed->SetValue(0, 0.1*(INSIDE_LABEL - OUTSIDE_LABEL));
				marchingFixed->Update();

				vtkSmartPointer< vtkOBBTreeCustomized > obbTree = vtkSmartPointer< vtkOBBTreeCustomized >::New();
				obbTree->SetDataSet(marchingFixed->GetOutput());
				obbTree->BuildLocator();
				
				vtkSmartPointer< vtkPolyData > obbCenters = vtkSmartPointer< vtkPolyData >::New();
				
				std::vector< vnl_vector_fixed< BaseRegistrationFilterType::MatrixElementType,3 > > regionCenters;
				std::vector< vnl_matrix_fixed < BaseRegistrationFilterType::MatrixElementType,3,3 > > regionSigmas;

				obbTree->GenerateRepresentationMeanPoints(level,obbCenters,regionCenters,regionSigmas,args.boxOverlapValues[level-args.startLevel]);

                filter->SetRegionCenters(regionCenters);
                filter->SetRegionSigmas(regionSigmas);
                */

                // covariance for all nodes plus one (prior)    
                if ( ! args.priorSigmaInverse.empty() ) {
                    filter->SetStartLevel(level);
                    filter->SetEndLevel(args.numberOfLevels-1);
                    level = args.numberOfLevels;
                }
                else {
                    filter->SetStartLevel(level);
                    filter->SetEndLevel(level);
                    
                    if(level > 0)
                        filter->SetExplainedVelocityField(polyaffineTreeTransform->GetParametersAsVectorField(level));
                }
                
                filter->SetPolyaffineTransform(polyaffineTreeTransform);                
//                filter->SetMeanSetAllLevels(MeanSetAllLevels);
//                filter->SetCovSetAllLevels(CovSetAllLevels);
                filter->SetMui(mui);
                filter->SetSigmaInverse(SigmaInverse);
//                filter->SetLambdaZero(LambdaZeroMean);
//                filter->SetCZero(CZeroMean);
				filter->SmoothVelocityFieldOn();
				filter->SmoothUpdateFieldOff();				
				//filter->SetIntensityDifferenceThreshold( 0.001 );
				
				filter->SetNumberOfIterations( 100 );
				filter->SetFixedImage( fixedImage );
				filter->SetMovingImage( movingImage );
				
				filter->SetFixedMaskImage(fixedMaskImage);
				filter->SetMovingMaskImage(movingMaskImage);
                
				if( inputVelField )
					filter->SetInitialVelocityField( inputVelField );
                else {
                    FieldType::Pointer zeroField = FieldType::New();
                    zeroField->CopyInformation(fixedImage);
                    zeroField->SetRegions(fixedImage->GetRequestedRegion());
                    zeroField->Allocate();
                    zeroField->FillBuffer(0.0);
                    
//                    typedef itk::ImageFileWriter< FieldType > FieldWriterType;
//                    typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
//                    fieldWriter->SetFileName(  "ZeroField.mhd" );
//                    fieldWriter->SetInput( zeroField );
//                    fieldWriter->Update();
                    
                    filter->SetInitialVelocityField(zeroField);
                }
                
				// Compute the deformation field
				try
				{
					//filter->UpdateLargestPossibleRegion();
                    filter->Update();
				}
				catch( itk::ExceptionObject& err )
				{
					std::cout << "Unexpected error." << std::endl;
					std::cout << err << std::endl;
					exit( EXIT_FAILURE );
				}
								
				std::cout << "Iter: " << filter->GetElapsedIterations() << "  ";
				std::cout << "Metric: "   << filter->GetMetric()   << "  ";
				std::cout << "RMSChange: " << filter->GetRMSChange() << "  ";
				std::cout << std::endl;
                
				// Get various outputs
				
				// Final deformation field
				defField = filter->GetDeformationField();
				defField->DisconnectPipeline();
				
				// Inverse final deformation field
				invDefField = filter->GetInverseDeformationField();
				invDefField->DisconnectPipeline();
				
				// Final velocity field
				inputVelField =  filter->GetVelocityField();
				inputVelField->DisconnectPipeline();
				
				metricVec.push_back(filter->GetMetric());
				
				if(!args.outputHarmonicEnergyFile.empty())
				{
					typedef itk::WarpHarmonicEnergyCalculator<FieldType> HarmonicEnergyCalculatorType;
					typename HarmonicEnergyCalculatorType::Pointer harmonicEnergyCalculator	= HarmonicEnergyCalculatorType::New();
					harmonicEnergyCalculator->SetImage( defField );
					harmonicEnergyCalculator->Compute();
					const double harmonicEnergy = harmonicEnergyCalculator->GetHarmonicEnergy();
					harmonicEnergyVec.push_back(harmonicEnergy);
				}
                
                // Write output inverse deformation field
                if (!args.outputPolyaffineTreeTransformFile.empty())
                {            
                    itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
                    transformWriter->SetFileName(args.outputPolyaffineTreeTransformFile.c_str());
                    transformWriter->SetInput(polyaffineTreeTransform);
                    
                    try
                    {
                        transformWriter->Update();
                    }
                    catch( itk::ExceptionObject& err )
                    {
                        std::cout << "Unexpected error." << std::endl;
                        std::cout << err << std::endl;
                        exit( EXIT_FAILURE );
                    }
                }

                // collect new prior
//                LambdaNewAll.push_back(filter->GetLambdaNew());
//                CNewAll.push_back(filter->GetCNew());
//                SigmaNewInverseAll.push_back(filter->GetSigmaNewInverse());
                				
				if(!args.outputCorrelationMatrixFile.empty())
				{
                    std::cout << "not yet implemented for gaussian mixture tree" << std::endl;
                    /*
					QString fileName = args.outputCorrelationMatrixFile.c_str();
					fileName.remove(".txt");
					fileName = fileName + "-level" + QString::number(level) + ".txt";	
					QFile fileCorrelationMatrix(fileName);
					if( !fileCorrelationMatrix.open(QIODevice::WriteOnly) )
						std::cout << "Failed to open file." << std::endl;
					QTextStream streamCorrelationMatrix(&fileCorrelationMatrix);
					vnl_matrix< BaseRegistrationFilterType::MatrixElementType > corr = filter->ComputeCorrelationsBetweenRegions();
					for(unsigned int m = 0; m < corr.rows(); ++m) {
						for(unsigned int n = 0; n < corr.cols(); ++n)
							streamCorrelationMatrix << corr(m,n) << " ";
						streamCorrelationMatrix << endl;
					}
					fileCorrelationMatrix.close();
					
					// visualize correlation using skeleton					
					vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
					//pts->Allocate(regionCenters.size());
					for(unsigned int i = 0; i < regionCenters.size(); ++i) {
						double point[3];
						regionCenters[i].copy_out(point);
						pts->InsertNextPoint(point);
					}

					//vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
					vtkSmartPointer<vtkFloatArray> colors = vtkSmartPointer<vtkFloatArray>::New();
					colors->SetNumberOfComponents(1);
					colors->SetName("correlation");
					
					vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
					
					for(unsigned int m = 0; m < corr.rows(); ++m) {
						for(unsigned int n = 0; n < corr.cols(); ++n) {
							
							//unsigned char grey[3] = {corr(m,n), corr(m,n), corr(m,n)};
							float c[1] = {corr(m,n)};
							colors->InsertNextTupleValue(c);
							
							vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
							line->GetPointIds()->SetId(0,m);
							line->GetPointIds()->SetId(1,n);
							
							lines->InsertNextCell(line);
						}
					}
					
					vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
					linesPolyData->SetPoints(pts);
					linesPolyData->SetLines(lines);
					linesPolyData->GetCellData()->SetScalars(colors);
					
					vtkSmartPointer< vtkPolyDataWriter > writer = vtkSmartPointer< vtkPolyDataWriter >::New();
					QString filename = QString("corr-level") + QString::number(level) + QString(".vtk");
					writer->SetFileName(filename.toAscii());
					writer->SetFileTypeToBinary();
					writer->SetInput(linesPolyData);
					writer->Update();
                     */
					
				}
							
                itkProbesStop(currLevelStr.str().c_str());
                itkProbesReport( std::cout );
    
        }
    
//        if (!args.outputCNewFile.empty() && !args.outputLambdaNewFile.empty())
//        {
//            vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > combinedCNew;
//            vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > combinedLambdaNew;
//            
//            if(CNewAll.size() > 1) {
//               
//                combinedCNew.set_size(CNewAll[0].rows(),CNewAll[0].cols());
//                combinedCNew.fill(0.0);
//                for(unsigned int i = 0; i < CNewAll.size(); ++i) {
//                    combinedCNew += CNewAll[i];
//                }
//                combinedCNew *= 1.0/CNewAll.size();
//                
//                unsigned int combinedNoOfRows = 0;
//                unsigned int combinedNoOfCols = 0;
//                for(unsigned int i = 0; i < LambdaNewAll.size(); ++i) {
//                    combinedNoOfRows += LambdaNewAll[i].rows();
//                    combinedNoOfCols += LambdaNewAll[i].cols();
//                }
//                
//                vnl_cholesky combinedCNewCholesky(combinedCNew);
//                vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > Id_N(combinedNoOfRows,combinedNoOfCols);
//                Id_N.set_identity();
//                
//                vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > combinedSigmaNewInverse(combinedNoOfRows*combinedCNew.rows(), combinedNoOfCols*combinedCNew.cols());
//                combinedSigmaNewInverse.fill(0.0);
//                unsigned int combinedSigmaNewInverseIndex = 0;
//                for(unsigned int i = 0; i < SigmaNewInverseAll.size(); ++i) {
//                    combinedSigmaNewInverse.update(SigmaNewInverseAll[i], combinedSigmaNewInverseIndex, combinedSigmaNewInverseIndex);
//                    combinedSigmaNewInverseIndex += SigmaNewInverseAll[i].rows();
//                }
//                
//                vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > A = combinedSigmaNewInverse * BaseRegistrationFilterType::KroneckerProduct(Id_N,combinedCNewCholesky.inverse());
//                
//                combinedLambdaNew.set_size(combinedNoOfRows,combinedNoOfCols);
//                combinedLambdaNew.fill(0.0);
//                for(unsigned int i = 0; i < combinedNoOfRows; ++i) {
//                    for(unsigned int j = 0; j < combinedNoOfCols; ++j) {
//                        combinedLambdaNew(i,j) = (1.0/(Dimension)*(Dimension+1)) * vnl_trace(A.extract((Dimension)*(Dimension+1), (Dimension)*(Dimension+1), i*(Dimension)*(Dimension+1), j*(Dimension)*(Dimension+1)));
//                    }
//                }
//                
////                combinedLambdaNew *= 0.01;
//                
//            }
//            else {
//                combinedCNew = CNewAll[0];
//                combinedLambdaNew = LambdaNewAll[0];
//            }
//            
//            combinedLambdaNew *= args.scaleLambdaPrior;
//            
//            vnl_matrix< PolyaffineTreeTransformType::MatrixElementType > CNew = filter->GetCNew();
//            std::ofstream cNewFile;
//            cNewFile.open (args.outputCNewFile.c_str());
//            cNewFile << combinedCNew << std::endl;
//            cNewFile.close();
//            
//            std::ofstream lambdaNewFile;
//            lambdaNewFile.open (args.outputLambdaNewFile.c_str());
//            lambdaNewFile << combinedLambdaNew << std::endl;
//            lambdaNewFile.close();
//        }
    
		// warp the result
		typedef itk::WarpImageFilter< ImageType, ImageType, FieldType >  WarperType;
		typename WarperType::Pointer warper = WarperType::New();
		warper->SetInput( movingImage );
		warper->SetOutputSpacing( fixedImage->GetSpacing() );
		warper->SetOutputOrigin( fixedImage->GetOrigin() );
		warper->SetOutputDirection( fixedImage->GetDirection() );
		warper->SetDisplacementField( defField );
        //warper->SetDeformationField( defField );
    
		
		// Write warped image out to file
		typedef PixelType                                OutputPixelType;
		typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
		typedef itk::CastImageFilter
		< ImageType, OutputImageType >                  CastFilterType;
		typedef itk::ImageFileWriter< OutputImageType >  WriterType;
		
		typename WriterType::Pointer      writer =  WriterType::New();
		typename CastFilterType::Pointer  caster =  CastFilterType::New();
		writer->SetFileName( args.outputImageFile.c_str() );
		caster->SetInput( warper->GetOutput() );
		writer->SetInput( caster->GetOutput()   );
		writer->SetUseCompression( true );
		
		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject& err )
		{
			std::cout << "Unexpected error." << std::endl;
			std::cout << err << std::endl;
			exit( EXIT_FAILURE );
		}
		
		
		// Write output deformation field
		if (!args.outputDeformationFieldFile.empty())
		{
			// Write the deformation field as an image of vectors.
			// Note that the file format used for writing the deformation field must be
			// capable of representing multiple components per pixel. This is the case
			// for the MetaImage and VTK file formats for example.
			typedef itk::ImageFileWriter< FieldType > FieldWriterType;
			typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
			fieldWriter->SetFileName(  args.outputDeformationFieldFile.c_str() );
			fieldWriter->SetInput( defField );
			fieldWriter->SetUseCompression( true );
			
			try
			{
				fieldWriter->Update();
			}
			catch( itk::ExceptionObject& err )
			{
				std::cout << "Unexpected error." << std::endl;
				std::cout << err << std::endl;
				exit( EXIT_FAILURE );
			}
		}
		
		// Write output inverse deformation field
		if (!args.outputInverseDeformationFieldFile.empty())
		{
			// Write the inverse deformation field as an image of vectors.
			// Note that the file format used for writing the inverse deformation field must be
			// capable of representing multiple components per pixel. This is the case
			// for the MetaImage and VTK file formats for example.
			typedef itk::ImageFileWriter< FieldType > FieldWriterType;
			typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
			fieldWriter->SetFileName(  args.outputInverseDeformationFieldFile.c_str() );
			fieldWriter->SetInput( invDefField );
			fieldWriter->SetUseCompression( true );
			
			try
			{
				fieldWriter->Update();
			}
			catch( itk::ExceptionObject& err )
			{
				std::cout << "Unexpected error." << std::endl;
				std::cout << err << std::endl;
				exit( EXIT_FAILURE );
			}
		}
		
		// Write output inverse deformation field
		if (!args.outputVelocityFieldFile.empty())
		{
			// Write the velocity field as an image of vectors.
			// Note that the file format used for writing the velocity field must be
			// capable of representing multiple components per pixel. This is the case
			// for the MetaImage and VTK file formats for example.
			typedef itk::ImageFileWriter< FieldType > FieldWriterType;
			typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
			fieldWriter->SetFileName(  args.outputVelocityFieldFile.c_str() );
			fieldWriter->SetInput( inputVelField );
			fieldWriter->SetUseCompression( true );
			
			try
			{
				fieldWriter->Update();
			}
			catch( itk::ExceptionObject& err )
			{
				std::cout << "Unexpected error." << std::endl;
				std::cout << err << std::endl;
				exit( EXIT_FAILURE );
			}
		}
		
		if(!args.outputMetricFile.empty())
		{
            std::ofstream metricFile;
            metricFile.open (args.outputMetricFile.c_str());
			for(unsigned int i = 0; i < metricVec.size(); ++i)
				metricFile << metricVec[i] << " ";
			metricFile << std::endl;
			metricFile.close();
		}
		
		if(!args.outputHarmonicEnergyFile.empty())
		{
            std::ofstream harmonicFile;
            harmonicFile.open (args.outputHarmonicEnergyFile.c_str());
			for(unsigned int i = 0; i < harmonicEnergyVec.size(); ++i)
				harmonicFile << harmonicEnergyVec[i] << " ";
			harmonicFile << std::endl;
			harmonicFile.close();
		}

    
    itkProbesStop( "Registration" );
    itkProbesReport( std::cout );
}
	
int main( int argc, char *argv[] )
{
  struct arguments args;
  parseOpts (argc, argv, args);

  std::cout<<"Starting demons registration with the following arguments:"<<std::endl;
  std::cout<<args<<std::endl<<std::endl;

  // FIXME uncomment for debug only
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);

  // Get the image dimension
  itk::ImageIOBase::Pointer imageIO;
  try
    {
    imageIO = itk::ImageIOFactory::CreateImageIO(
       args.fixedImageFile.c_str(), itk::ImageIOFactory::ReadMode);
    if ( imageIO )
      {
      imageIO->SetFileName(args.fixedImageFile.c_str());
      imageIO->ReadImageInformation();
      }
    else
      {
      std::cout << "Could not read the fixed image information." << std::endl;
      exit( EXIT_FAILURE );
      }
    }
  catch( itk::ExceptionObject& err )
    {
    std::cout << "Could not read the fixed image information." << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
    }
  
  switch ( imageIO->GetNumberOfDimensions() )
  {
  case 3:
    LogDomainDemonsRegistrationFunction<3>(args);
    break;
  default:
    std::cout << "Unsuported dimension" << std::endl;
    exit( EXIT_FAILURE );
  }
  
  return EXIT_SUCCESS;
}

#undef itkProbesCreate
#undef itkProbesStart
#undef itkProbesStop
#undef itkProbesReport
