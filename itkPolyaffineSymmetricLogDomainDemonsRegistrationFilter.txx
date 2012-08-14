#ifndef __itkPolyaffineSymmetricLogDomainDemonsRegistrationFilter_txx
#define __itkPolyaffineSymmetricLogDomainDemonsRegistrationFilter_txx

#include "itkPolyaffineSymmetricLogDomainDemonsRegistrationFilter.h"

#include "itkOppositeImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkVelocityFieldBCHCompositionFilter.h"

namespace itk {

// Default constructor
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::PolyaffineSymmetricLogDomainDemonsRegistrationFilter()
{
  DemonsRegistrationFunctionPointer drfpf = DemonsRegistrationFunctionType::New();
  this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
                                 drfpf.GetPointer() ) );

  DemonsRegistrationFunctionPointer drfpb = DemonsRegistrationFunctionType::New();
  this->SetBackwardDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
                                          drfpb.GetPointer() ) );
  
  m_Multiplier = MultiplyByConstantType::New();
  m_Multiplier->InPlaceOn();
 
  m_Adder = AdderType::New();
  m_Adder->InPlaceOn();

  // Set number of terms in the BCH approximation to default value
  m_NumberOfBCHApproximationTerms = 2;

  m_BackwardUpdateBuffer = 0;
}


// Checks whether the DifferenceFunction is of type DemonsRegistrationFunction.
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
typename PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::DemonsRegistrationFunctionType*
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GetForwardRegistrationFunctionType()
{
  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>(this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
    {
    itkExceptionMacro( << "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
    }

  return drfp;
}


// Checks whether the DifferenceFunction is of type DemonsRegistrationFunction.
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
const typename PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::DemonsRegistrationFunctionType*
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GetForwardRegistrationFunctionType() const
{
  const DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<const DemonsRegistrationFunctionType *>(this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
    {
    itkExceptionMacro( << "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
    }

  return drfp;
}


// Checks whether the DifferenceFunction is of type DemonsRegistrationFunction.
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
typename PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::DemonsRegistrationFunctionType*
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GetBackwardRegistrationFunctionType()
{
  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>(this->GetBackwardDifferenceFunction().GetPointer());
 
  if( !drfp )
    {
    itkExceptionMacro( << "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
    }

  return drfp;
}


// Checks whether the DifferenceFunction is of type DemonsRegistrationFunction.
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
const typename PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::DemonsRegistrationFunctionType*
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GetBackwardRegistrationFunctionType() const
{
  const DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<const DemonsRegistrationFunctionType *>(this->GetBackwardDifferenceFunction().GetPointer());
 
  if( !drfp )
    {
    itkExceptionMacro( << "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
    }

  return drfp;
}


// Set the function state values before each iteration
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::InitializeIteration()
{
  // update variables in the equation object
  DemonsRegistrationFunctionType *f = this->GetForwardRegistrationFunctionType();
  f->SetDisplacementField( this->GetDeformationField() );
  //f->SetDeformationField( this->GetDeformationField() );

  DemonsRegistrationFunctionType *b = this->GetBackwardRegistrationFunctionType();
  b->SetFixedImage( this->GetMovingImage() );
  b->SetMovingImage( this->GetFixedImage() );
  b->SetDisplacementField( this->GetInverseDeformationField() );
  //b->SetDeformationField( this->GetInverseDeformationField() );
  b->InitializeIteration();
  
  // call the superclass  implementation ( initializes f )
  Superclass::InitializeIteration();
}


// Initialize flags
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::Initialize()
{
  DemonsRegistrationFunctionType *f = this->GetForwardRegistrationFunctionType();
  f->SetFixedMaskImage(this->GetFixedMaskImage());
  f->SetMovingMaskImage(this->GetMovingMaskImage());
	
  DemonsRegistrationFunctionType *b = this->GetBackwardRegistrationFunctionType();
  b->SetFixedMaskImage(this->GetMovingMaskImage());
  b->SetMovingMaskImage(this->GetFixedMaskImage());
	
  //std::cout<<"PolyaffineSymmetricLogDomainDemonsRegistrationFilter::Initialize"<<std::endl;
  this->Superclass::Initialize();

  const FixedImageType * fixim = this->GetFixedImage();
  const MovingImageType * movim = this->GetMovingImage();

  if ( fixim==0 || movim==0 )
    {
    itkExceptionMacro( << "A fixed and a moving image are required" );
    }

  if ( fixim->GetLargestPossibleRegion() != movim->GetLargestPossibleRegion() )
    {
    itkExceptionMacro( << "Registering images that have diffent sizes is not supported yet." );
    }

  if ( (fixim->GetSpacing() - movim->GetSpacing()).GetNorm() > 1e-10 )
    {
    itkExceptionMacro( << "Registering images that have diffent spacing is not supported yet." );
    }

  if ( (fixim->GetOrigin() - movim->GetOrigin()).GetNorm() > 1e-10 )
    {
    itkExceptionMacro( << "Registering images that have diffent origins is not supported yet." );
    }
}


// Get the metric value from the difference function
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
double
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GetMetric() const
{
  const DemonsRegistrationFunctionType *drfpf = this->GetForwardRegistrationFunctionType();
  const DemonsRegistrationFunctionType *drfpb = this->GetBackwardRegistrationFunctionType();
  return 0.5*(drfpf->GetMetric() + drfpb->GetMetric());
}


// Get the metric value from the difference function
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
double
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GetPreviousMetric() const
{
    const DemonsRegistrationFunctionType *drfpf = this->GetForwardRegistrationFunctionType();
    const DemonsRegistrationFunctionType *drfpb = this->GetBackwardRegistrationFunctionType();
    return 0.5*(drfpf->GetPreviousMetric() + drfpb->GetPreviousMetric());
}
    

// Get Intensity Difference Threshold
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
double
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GetIntensityDifferenceThreshold() const
{
  const DemonsRegistrationFunctionType *drfpf = this->GetForwardRegistrationFunctionType();
  const DemonsRegistrationFunctionType *drfpb = this->GetBackwardRegistrationFunctionType();

  if ( drfpf->GetIntensityDifferenceThreshold() != drfpb->GetIntensityDifferenceThreshold())
  {
     itkExceptionMacro(<<"Forward and backward FiniteDifferenceFunctions not in sync");
  }
  return drfpf->GetIntensityDifferenceThreshold();
}


// Set Intensity Difference Threshold
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::SetIntensityDifferenceThreshold(double threshold) 
{
  DemonsRegistrationFunctionType *drfpf = this->GetForwardRegistrationFunctionType();
  DemonsRegistrationFunctionType *drfpb = this->GetBackwardRegistrationFunctionType();
  drfpf->SetIntensityDifferenceThreshold(threshold);
  drfpb->SetIntensityDifferenceThreshold(threshold);
}


// Set Maximum Update Step Length
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::SetMaximumUpdateStepLength(double step)
{
  DemonsRegistrationFunctionType *drfpf = this->GetForwardRegistrationFunctionType();
  DemonsRegistrationFunctionType *drfpb = this->GetBackwardRegistrationFunctionType();
  drfpf->SetMaximumUpdateStepLength(step);
  drfpb->SetMaximumUpdateStepLength(step);
}


// Get Maximum Update Step Length
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
double
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GetMaximumUpdateStepLength() const
{
  const DemonsRegistrationFunctionType *drfpf = this->GetForwardRegistrationFunctionType();
  const DemonsRegistrationFunctionType *drfpb = this->GetBackwardRegistrationFunctionType();

  if ( drfpf->GetMaximumUpdateStepLength() != drfpb->GetMaximumUpdateStepLength())
  {
     itkExceptionMacro(<<"Forward and backward FiniteDifferenceFunctions not in sync");
  }
  return drfpf->GetMaximumUpdateStepLength();
}


// Get gradient type  
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
typename PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>::GradientType
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::GetUseGradientType() const
{
  const DemonsRegistrationFunctionType *drfpf = this->GetForwardRegistrationFunctionType();
  const DemonsRegistrationFunctionType *drfpb = this->GetBackwardRegistrationFunctionType();

  if ( drfpf->GetUseGradientType() != drfpb->GetUseGradientType())
  {
     itkExceptionMacro(<<"Forward and backward FiniteDifferenceFunctions not in sync");
  }
  return drfpf->GetUseGradientType();
}


// Set gradient type  
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::SetUseGradientType(GradientType gtype) 
{
  DemonsRegistrationFunctionType *drfpf = this->GetForwardRegistrationFunctionType();
  DemonsRegistrationFunctionType *drfpb = this->GetBackwardRegistrationFunctionType();
  drfpf->SetUseGradientType(gtype);
  drfpb->SetUseGradientType(gtype);
}


// Allocate storage in m_UpdateBuffer
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::AllocateUpdateBuffer()
{
  Superclass::AllocateUpdateBuffer();

  this->AllocateBackwardUpdateBuffer();
}


// Allocates storage in m_BackwardUpdateBuffer
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::AllocateBackwardUpdateBuffer()
{
  if ( m_NumberOfBCHApproximationTerms < 3 )
    {
    m_BackwardUpdateBuffer = 0;
    return;
    }
  
  // The backward update buffer looks just like the output.
  VelocityFieldPointer output = this->GetOutput();

  if ( !m_BackwardUpdateBuffer ) m_BackwardUpdateBuffer = VelocityFieldType::New();
  m_BackwardUpdateBuffer->SetOrigin(output->GetOrigin());
  m_BackwardUpdateBuffer->SetSpacing(output->GetSpacing());
  m_BackwardUpdateBuffer->SetDirection(output->GetDirection());
  m_BackwardUpdateBuffer->SetLargestPossibleRegion(output->GetLargestPossibleRegion());
  m_BackwardUpdateBuffer->SetRequestedRegion(output->GetRequestedRegion());
  m_BackwardUpdateBuffer->SetBufferedRegion(output->GetBufferedRegion());
  m_BackwardUpdateBuffer->Allocate();
}


// Smooth the backward update field using a separable Gaussian kernel
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::SmoothBackwardUpdateField()
{
  // The update buffer will be overwritten with new data.
  this->SmoothGivenField(this->GetBackwardUpdateBuffer(), this->GetUpdateFieldStandardDeviations());
}


template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
typename
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>::TimeStepType
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::ThreadedCalculateChange(const ThreadRegionType &regionToProcess, int)
{
  typedef typename VelocityFieldType::RegionType      RegionType;
  typedef typename VelocityFieldType::SizeType        SizeType;
  typedef typename VelocityFieldType::SizeValueType   SizeValueType;
  typedef typename VelocityFieldType::IndexType       IndexType;
  typedef typename VelocityFieldType::IndexValueType  IndexValueType;
  typedef typename 
    FiniteDifferenceFunctionType::NeighborhoodType    NeighborhoodIteratorType;
  typedef ImageRegionIterator<VelocityFieldType>      UpdateIteratorType;

  VelocityFieldPointer output = this->GetOutput();

  // Get the FiniteDifferenceFunction to use in calculations.
  const typename FiniteDifferenceFunctionType::Pointer dff
    = this->GetDifferenceFunction();
  const typename FiniteDifferenceFunctionType::Pointer dfb
    = this->GetBackwardDifferenceFunction();

  if ( dff->GetRadius() != dfb->GetRadius())
  {
     itkExceptionMacro(<<"Forward and backward FiniteDifferenceFunctions not in sync");
  }
  
  const SizeType  radius = dff->GetRadius();
  
  // Break the input into a series of regions.  The first region is free
  // of boundary conditions, the rest with boundary conditions.  We operate
  // on the output region because input has been copied to output.
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<VelocityFieldType>
    FaceCalculatorType;

  typedef typename FaceCalculatorType::FaceListType FaceListType;

  FaceCalculatorType faceCalculator;

  FaceListType faceList = faceCalculator(output, regionToProcess, radius);
  typename FaceListType::iterator fIt = faceList.begin();

  // Ask the function object for a pointer to a data structure it
  // will use to manage any global values it needs.  We'll pass this
  // back to the function object at each calculation and then
  // again so that the function object can use it to determine a
  // time step for this iteration.
  void *globalDataf = dff->GetGlobalDataPointer();
  void *globalDatab = dfb->GetGlobalDataPointer();

  // Process the non-boundary region.
  NeighborhoodIteratorType nD(radius, output, *fIt);
  if ( m_NumberOfBCHApproximationTerms == 2 )
    {
    UpdateIteratorType       nU(this->GetUpdateBuffer(),  *fIt);
    while( !nD.IsAtEnd() )
      {
      nU.Value() = (dff->ComputeUpdate(nD, globalDataf)-dfb->ComputeUpdate(nD, globalDatab))*0.5;
      ++nD;
      ++nU;
      }
    
    // Process each of the boundary faces.
    NeighborhoodIteratorType bD;
    UpdateIteratorType   bU;
    for (++fIt; fIt != faceList.end(); ++fIt)
      {
      bD = NeighborhoodIteratorType(radius, output, *fIt);
      bU = UpdateIteratorType  (this->GetUpdateBuffer(), *fIt);
      while ( !bD.IsAtEnd() )
        {
        bU.Value() = (dff->ComputeUpdate(bD, globalDataf)-dfb->ComputeUpdate(bD, globalDatab))*0.5;
        ++bD;
        ++bU;
        }
      }
    }
  else
    {
    UpdateIteratorType       nUF(this->GetUpdateBuffer(),  *fIt);
    UpdateIteratorType       nUB(this->GetBackwardUpdateBuffer(),  *fIt);
    while( !nD.IsAtEnd() )
      {
      nUF.Value() = dff->ComputeUpdate(nD, globalDataf);
      nUB.Value() = dfb->ComputeUpdate(nD, globalDatab);
      ++nD;
      ++nUF;
      ++nUB;
      }
    
    // Process each of the boundary faces.
    NeighborhoodIteratorType bD;
    UpdateIteratorType   bUF;
    UpdateIteratorType   bUB;
    for (++fIt; fIt != faceList.end(); ++fIt)
      {
      bD = NeighborhoodIteratorType(radius, output, *fIt);
      bUF = UpdateIteratorType  (this->GetUpdateBuffer(), *fIt);
      bUB = UpdateIteratorType  (this->GetBackwardUpdateBuffer(), *fIt);
      while ( !bD.IsAtEnd() )
        {
        bUF.Value() = dff->ComputeUpdate(bD, globalDataf);
        bUB.Value() = dfb->ComputeUpdate(bD, globalDatab);
        ++bD;
        ++bUF;
        ++bUB;
        }
      }
    }

  // Ask the finite difference function to compute the time step for
  // this iteration.  We give it the global data pointer to use, then
  // ask it to free the global data memory.
  TimeStepType timeStep = 0.5*( dff->ComputeGlobalTimeStep(globalDataf)
                                + dfb->ComputeGlobalTimeStep(globalDatab) );
  dff->ReleaseGlobalDataPointer(globalDataf);
  dfb->ReleaseGlobalDataPointer(globalDatab);

  return timeStep;
}


// Get the metric value from the difference function
template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::ApplyUpdate(const TimeStepType& dt)
{
  const DemonsRegistrationFunctionType *drfpf = this->GetForwardRegistrationFunctionType();
  const DemonsRegistrationFunctionType *drfpb = this->GetBackwardRegistrationFunctionType();
  this->SetRMSChange( 0.5*(drfpf->GetRMSChange() + drfpb->GetRMSChange()) );

  if ( this->m_NumberOfBCHApproximationTerms < 3 )
    {
    // If we smooth the update buffer before applying it, then the are
    // approximating a viscuous problem as opposed to an elastic problem
    if ( this->GetSmoothUpdateField() )
      {
      this->SmoothUpdateField();
      }
  
    // Use time step if necessary. In many cases
    // the time step is one so this will be skipped
    if ( fabs(dt - 1.0)>1.0e-4 )
      {
      itkDebugMacro( "Using timestep: " << dt );
      m_Multiplier->SetConstant( dt );
      m_Multiplier->SetInput( this->GetUpdateBuffer() );
      m_Multiplier->GraftOutput( this->GetUpdateBuffer() );
      // in place update
      m_Multiplier->Update();
      // graft output back to this->GetUpdateBuffer()
      this->GetUpdateBuffer()->Graft( m_Multiplier->GetOutput() );
      }


    // Apply update
    m_Adder->SetInput( 0, this->GetOutput() );
    m_Adder->SetInput( 1, this->GetUpdateBuffer() );
    m_Adder->GraftOutput( this->GetOutput() );
    m_Adder->GetOutput()->SetRequestedRegion( this->GetOutput()->GetRequestedRegion() );
    
    // Triggers in place update
    m_Adder->Update();
    
    // Region passing stuff
    this->GraftOutput( m_Adder->GetOutput() );
    }
  else
    {
    // If we smooth the update buffer before applying it, then the are
    // approximating a viscuous problem as opposed to an elastic problem
    if ( this->GetSmoothUpdateField() )
      {
      this->SmoothUpdateField();
      this->SmoothBackwardUpdateField();
      }
  
    // Use time step if necessary. In many cases
    // the time step is one so this will be skipped
    if ( fabs(dt - 1.0)>1.0e-4 )
      {
      itkDebugMacro( "Using timestep: " << dt );
      m_Multiplier->SetConstant( dt );
      
      m_Multiplier->SetInput( this->GetUpdateBuffer() );
      m_Multiplier->GraftOutput( this->GetUpdateBuffer() );
      // in place update
      m_Multiplier->Update();
      // graft output back to this->GetUpdateBuffer()
      this->GetUpdateBuffer()->Graft( m_Multiplier->GetOutput() );

      m_Multiplier->SetInput( this->GetBackwardUpdateBuffer() );
      m_Multiplier->GraftOutput( this->GetBackwardUpdateBuffer() );
      // in place update
      m_Multiplier->Update();
      // graft output back to this->GetUpdateBuffer()
      this->GetBackwardUpdateBuffer()->Graft( m_Multiplier->GetOutput() );
      }

    // Apply update (declare the filters here as efficiency is not critical
    // with "high" order BCH approximations)
    typedef VelocityFieldBCHCompositionFilter<
       VelocityFieldType, VelocityFieldType>   BCHFilterType;
    
    typename BCHFilterType::Pointer bchfilter = BCHFilterType::New();
    bchfilter->SetNumberOfApproximationTerms( this->m_NumberOfBCHApproximationTerms );

    // First get Z( v, K_fluid * u_forward )
    bchfilter->SetInput( 0, this->GetOutput() );
    bchfilter->SetInput( 1, this->GetUpdateBuffer() );

    bchfilter->GetOutput()->SetRequestedRegion( this->GetOutput()->GetRequestedRegion() );
    bchfilter->Update();
    VelocityFieldPointer Zf = bchfilter->GetOutput();
    Zf->DisconnectPipeline();

    // Now get Z( -v, K_fluid * u_backward )
    typedef OppositeImageFilter<
       VelocityFieldType, VelocityFieldType>  OppositeFilterType;

    typename OppositeFilterType::Pointer oppositefilter = OppositeFilterType::New();
    oppositefilter->SetInput( this->GetOutput() );
    oppositefilter->InPlaceOn();

    bchfilter->SetInput( 0, oppositefilter->GetOutput() );
    bchfilter->SetInput( 1, this->GetBackwardUpdateBuffer() );

    bchfilter->GetOutput()->SetRequestedRegion( this->GetOutput()->GetRequestedRegion() );
    bchfilter->Update();
    VelocityFieldPointer Zb = bchfilter->GetOutput();
    Zb->DisconnectPipeline();

    // Finally get 0.5*( Z( v, K_fluid * u_forward ) - Z( -v, K_fluid * u_backward ) )
    typedef SubtractImageFilter<
       VelocityFieldType, VelocityFieldType, VelocityFieldType>  SubtracterType;

    typename SubtracterType::Pointer subtracter = SubtracterType::New();
    subtracter->SetInput( 0, Zf );
    subtracter->SetInput( 1, Zb );

    subtracter->GraftOutput( this->GetOutput() );
    
    m_Multiplier->SetConstant( 0.5 );
    m_Multiplier->SetInput( subtracter->GetOutput() );
    m_Multiplier->GraftOutput( this->GetOutput() );
    
    // Triggers in place update
    m_Multiplier->GetOutput()->SetRequestedRegion( this->GetOutput()->GetRequestedRegion() );
    m_Multiplier->Update();
    
    // Region passing stuff
    this->GraftOutput( m_Multiplier->GetOutput() );
    }


  // Smooth the velocity field
  if( this->GetSmoothVelocityField() )
    {
    this->SmoothVelocityField();
    }

}


template <class TFixedImage, class TMovingImage, class TField, class TFixedMaskImage, class TMovingMaskImage, class TWeightImage>
void
PolyaffineSymmetricLogDomainDemonsRegistrationFilter<TFixedImage,TMovingImage,TField,TFixedMaskImage,TMovingMaskImage,TWeightImage>
::PrintSelf(std::ostream& os, Indent indent) const
{ 
  Superclass::PrintSelf( os, indent );

  os << indent << "Intensity difference threshold: " << this->GetIntensityDifferenceThreshold() << std::endl;
  os << indent << "Multiplier: " << m_Multiplier << std::endl;
  os << indent << "Adder: " << m_Adder << std::endl;
  os << indent << "NumberOfBCHApproximationTerms: " << m_NumberOfBCHApproximationTerms << std::endl;
}


} // end namespace itk

#endif
