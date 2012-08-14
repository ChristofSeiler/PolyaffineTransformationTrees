#ifndef _itkManualPolyaffineStationaryVelocityFieldTransform_h_
#define _itkManualPolyaffineStationaryVelocityFieldTransform_h_

#include "itkPolyaffineStationaryVelocityFieldTransform.h"

#include <itkTransform.h>
#include <itkProcessObject.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkImage.h>
#include <itkTreeNode.h>
#include <itkMatrixOffsetTransformBase.h>

//#include <itkStationaryVelocityFieldExponential.h>
#include "itkExponentialDeformationFieldImageFilter2.h"

namespace itk
{


/**
 *
 * This class defines the non-linear stationary velocity field transformation.
 *
 * There are two templates for this class:
 *
 *   TScalarType  Type of the transformation parameters. Can be "float" or "double".
 *
 *   NDimensions  Dimension of the transformation space.
 *
 * The parameters of a stationary velocity field is a vector field. Since a vector field cannot
 * be simply represented by an simple linear array, we store this vector field into a VectorFieldType
 * which is simply an itk::Image where the type of voxels is itk::Vector. The methods
 * SetParametersAsVectorField and GetParametersAsVectorField respectively allow to set and get the
 * vector field.
 *
 * There is no I/O support as there is no parameters, the user might save the
 * vector field instead.
 *
 * @brief   Stationary velocity field transformation
 *
 * @class   PolyaffineStationaryVelocityFieldTransform (itk)
 *
 * @author  Vincent Garcia, Marco Lorenzi, Asclepios team, INRIA
 *
 */
template <class TScalarType=float, class TReferenceType=signed short, unsigned int NDimensions=3>
class ITK_EXPORT ManualPolyaffineStationaryVelocityFieldTransform : public PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>
{


public:

    typedef ManualPolyaffineStationaryVelocityFieldTransform                  Self;
    typedef PolyaffineStationaryVelocityFieldTransform<TScalarType, TReferenceType, NDimensions>  Superclass;
    typedef SmartPointer<Self>                                Pointer;
    typedef SmartPointer<const Self>                          ConstPointer;

    /**
     * Type of the input parameters
     */
    typedef typename Superclass::ParametersType               ParametersType;
    
    /**
     * Generic constructors.
     */
    itkNewMacro( Self );
    itkTypeMacro( ManualPolyaffineStationaryVelocityFieldTransform, PolyaffineStationaryVelocityFieldTransform );
    
    /**
     * Dimension of the domain space
     */
    itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
    itkStaticConstMacro(ParametersDimension, unsigned int, NDimensions);
    
    static const unsigned int NDimensionsHom = 4;
    typedef double MatrixElementType;
    typedef itk::MatrixOffsetTransformBase<MatrixElementType,NDimensions,NDimensions> MatrixTransformType;
    typedef vnl_vector< MatrixElementType > NodeValueType;
    typedef TreeNode< NodeValueType > GaussianTreeNodeType;
    typedef typename GaussianTreeNodeType::Pointer GaussianTreeNodePointerType;
    
    virtual void SetParametersToNode(const GaussianTreeNodePointerType node, const ParametersType& parameters, unsigned int& currentIndex, unsigned int currentLevel);
    
    virtual unsigned int GetNumberOfLevels() const;

protected:

    /**
     * Default constructor.
     */
    ManualPolyaffineStationaryVelocityFieldTransform(void);

    /**
     * Destructor.
     */
    virtual ~ManualPolyaffineStationaryVelocityFieldTransform(){};
    
};


} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkManualPolyaffineStationaryVelocityFieldTransform.txx"
#endif

#endif
