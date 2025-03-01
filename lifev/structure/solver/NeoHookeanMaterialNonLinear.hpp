//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief This file contains the definition for the St. Venant Kirchhoff linear material
 *
 *  @version 1.0
 *  @date 29-07-2010
 *  @author Paolo Tricerri,
 *  @author Gianmarco Mengaldo
 *
 *  @maintainer  Paolo Tricerri      <paolo.tricerri@epfl.ch>
 *  @contributor Gianmarco Mengaldo  <gianmarco.mengaldo@gmail.com>
 */

#ifndef _NEOHOOKEANMATERIAL_H_
#define _NEOHOOKEANMATERIAL_H_


#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>

namespace LifeV
{

template <typename MeshType>
class NeoHookeanMaterialNonLinear :
    public StructuralConstitutiveLaw<MeshType>
{
    //!@name Type definitions
    //@{

public:
    typedef StructuralConstitutiveLaw<MeshType>          super;

    typedef StructuralConstitutiveLawData            data_Type;

    typedef typename super::vector_Type              vector_Type;
    typedef typename super::matrix_Type              matrix_Type;

    typedef typename super::matrixPtr_Type           matrixPtr_Type;
    typedef typename super::vectorPtr_Type           vectorPtr_Type;
    typedef typename super::dataPtr_Type             dataPtr_Type;
    typedef typename super::displayerPtr_Type        displayerPtr_Type;

    typedef typename super::mapMarkerVolumesPtr_Type mapMarkerVolumesPtr_Type;
    typedef typename super::mapMarkerVolumes_Type mapMarkerVolumes_Type;
    typedef typename mapMarkerVolumes_Type::const_iterator mapIterator_Type;

    typedef typename super::vectorVolumes_Type       vectorVolumes_Type;
    typedef boost::shared_ptr<vectorVolumes_Type>    vectorVolumesPtr_Type;

    typedef std::vector<UInt>                        vectorIndexes_Type;
    typedef boost::shared_ptr<vectorIndexes_Type>    vectorIndexesPtr_Type;

    typedef typename super::mapMarkerIndexesPtr_Type mapMarkerIndexesPtr_Type;
    typedef typename super::mapMarkerIndexes_Type    mapMarkerIndexes_Type;
    typedef typename mapMarkerIndexes_Type::const_iterator mapIteratorIndex_Type;

    typedef typename super::FESpacePtr_Type          FESpacePtr_Type;
    typedef typename super::ETFESpacePtr_Type        ETFESpacePtr_Type;

    //Vector for vector parameters
    typedef typename super::vectorsParameters_Type       vectorsParameters_Type;
    typedef typename super::vectorsParametersPtr_Type    vectorsParametersPtr_Type;

    typedef MatrixSmall<3, 3>                          matrixSmall_Type;
    //@}



    //! @name Constructor &  Destructor
    //@{

    NeoHookeanMaterialNonLinear();

    virtual  ~NeoHookeanMaterialNonLinear();

    //@}

    //!@name Methods
    //@{

    //! Setup the created object of the class StructuralConstitutiveLaw
    /*!
      \param dFespace: the FiniteElement Space
      \param monolithicMap: the MapEpetra
      \param offset: the offset parameter used assembling the matrices
    */
    void setup ( const FESpacePtr_Type& dFESpace,
                 const ETFESpacePtr_Type& dETFESpace,
                 const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                 const UInt offset, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer );


    //! Compute the Stiffness matrix in StructuralSolver::buildSystem()
    /*!
      \param dataMaterial the class with Material properties data
    */
    void computeLinearStiff ( dataPtr_Type& /*dataMaterial*/, const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                              const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/ );


    //! Updates the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateJacobianMatrix ( const vector_Type& disp,
                                const dataPtr_Type& dataMaterial,
                                const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                const displayerPtr_Type& displayer);

   void updateJacobianMatrix ( const vector_Type& disp, const vector_Type& pressure, const dataPtr_Type& dataMaterial );



    //! Updates the nonlinear terms in the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param stiff: stiffness matrix provided from outside
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateNonLinearJacobianTerms ( matrixPtr_Type& jacobian,
                                        const vector_Type& disp,
                                        const dataPtr_Type& dataMaterial,
                                        const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                        const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/,
                                        const displayerPtr_Type& displayer);

    //! Interface method to compute the new Stiffness matrix in StructuralSolver::evalResidual and in
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void computeStiffness ( const vector_Type& disp, Real factor, const dataPtr_Type& dataMaterial,
                            const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                            const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/,
                            const displayerPtr_Type& displayer );

    void computeStiffness ( const vector_Type& sol, const vector_Type& pressure,
        		                         const dataPtr_Type& dataMaterial );

    //! Computes the new Stiffness vector for Neo-Hookean and Exponential materials in
    //! StructuralSolver given a certain displacement field.
    //! This function is used both in StructuralSolver::evalResidual and in StructuralSolver::updateSystem
    //! since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    // void computeVector( const vector_Type& sol,
    //                     Real factor,
    //                     const dataPtr_Type& dataMaterial,
    //                     const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
    //                     const displayerPtr_Type& displayer );


    //! Computes the deformation gradient F, the cofactor matrix Cof(F),
    //! the determinant of F (J = det(F)), the trace of right Cauchy-Green tensor tr(C)
    //! This function is used in StructuralConstitutiveLaw::computeStiffness
    /*!
      \param dk_loc: the elemental displacement
    */
    void computeKinematicsVariables ( const VectorElemental& dk_loc ) {}

    //! ShowMe method of the class (saved on a file the stiffness vector and the jacobian)
    void showMe ( std::string const& fileNameVectStiff,
                  std::string const& fileNameJacobain);


    //! Compute the First Piola Kirchhoff Tensor
    /*!
       \param firstPiola Epetra_SerialDenseMatrix that has to be filled
       \param tensorF Epetra_SerialDenseMatrix the deformation gradient
       \param cofactorF Epetra_SerialDenseMatrix cofactor of F
       \param invariants std::vector with the invariants of C and the detF
       \param material UInt number to get the material parameteres form the VenantElasticData class
    */
    void computeLocalFirstPiolaKirchhoffTensor ( Epetra_SerialDenseMatrix& firstPiola,
                                                 const Epetra_SerialDenseMatrix& tensorF,
                                                 const Epetra_SerialDenseMatrix& cofactorF,
                                                 const std::vector<Real>& invariants,
                                                 const UInt marker);


    //@}

    //! @name Get Methods
    //@{

    //! Get the Stiffness matrix
    matrixPtr_Type const stiffMatrix() const
    {
        return super::M_jacobian;
    }

    //! Get the stiffness vector
    vectorPtr_Type const stiffVector() const
    {
        return M_stiff;
    }

    void apply ( const vector_Type& sol, vector_Type& res,
                 const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                 const mapMarkerIndexesPtr_Type mapsMarkerIndexes);

    //@}



protected:

    //! construct the vectors for the parameters
    /*!
      \param VOID
      \return VOID
    */
    void setupVectorsParameters ( void );

    //! Vector: stiffness non-linear
    vectorPtr_Type                      M_stiff;

    //Create the indentity for F
    matrixSmall_Type                      M_identity;

};

template <typename MeshType>
NeoHookeanMaterialNonLinear<MeshType>::NeoHookeanMaterialNonLinear() :
    super           ( ),
    M_stiff         ( ),
    M_identity        ( )
{
}





template <typename MeshType>
NeoHookeanMaterialNonLinear<MeshType>::~NeoHookeanMaterialNonLinear()
{}





template <typename MeshType>
void
NeoHookeanMaterialNonLinear<MeshType>::setup ( const FESpacePtr_Type&                      dFESpace,
                                               const ETFESpacePtr_Type&                    dETFESpace,
                                               const boost::shared_ptr<const MapEpetra>&   monolithicMap,
                                               const UInt                                  offset,
                                               const dataPtr_Type&                         dataMaterial,
                                               const displayerPtr_Type&                    displayer)
{
    this->M_displayer = displayer;
    this->M_dataMaterial  = dataMaterial;

    //    std::cout<<"I am setting up the Material"<<std::endl;

    this->M_dispFESpace                     = dFESpace;
    this->M_dispETFESpace                     = dETFESpace;
    this->M_localMap                    = monolithicMap;
    this->M_offset                      = offset;
    this->M_dataMaterial                = dataMaterial;
    this->M_displayer                   = displayer;
    M_stiff.reset                   ( new vector_Type (*this->M_localMap) );

    M_identity (0, 0) = 1.0;
    M_identity (0, 1) = 0.0;
    M_identity (0, 2) = 0.0;
    M_identity (1, 0) = 0.0;
    M_identity (1, 1) = 1.0;
    M_identity (1, 2) = 0.0;
    M_identity (2, 0) = 0.0;
    M_identity (2, 1) = 0.0;
    M_identity (2, 2) = 1.0;

    // The 2 is because the law uses three parameters (mu, bulk).
    // another way would be to set up the number of constitutive parameters of the law
    // in the data file to get the right size. Note the comment below.
    this->M_vectorsParameters.reset ( new vectorsParameters_Type ( 2 ) );

    this->setupVectorsParameters();
}

template <typename MeshType>
void
NeoHookeanMaterialNonLinear<MeshType>::setupVectorsParameters ( void )
{
    // Paolo Tricerri: February, 20th
    // In each class, the name of the parameters has to inserted in the law
    // TODO: move the saving of the material parameters to more abstract objects
    //       such that in the class of the material we do not need to call explicitly
    //       the name of the parameter.

    // Number of volume on the local part of the mesh
    UInt nbElements = this->M_dispFESpace->mesh()->numVolumes();

    // Parameter mu
    // 1. resize the vector in the first element of the vector.
    (* (this->M_vectorsParameters) ) [0].resize ( nbElements );

    // Parameter bulk
    (* (this->M_vectorsParameters) ) [1].resize ( nbElements );

    for (UInt i (0); i < nbElements; i++ )
    {
        // Extracting the marker
        UInt markerID = this->M_dispFESpace->mesh()->element ( i ).markerID();

        Real mu = this->M_dataMaterial->mu ( markerID );
        Real bulk = this->M_dataMaterial->bulk ( markerID );

        ( (* (this->M_vectorsParameters) ) [0]) [ i ] = mu;
        ( (* (this->M_vectorsParameters) ) [1]) [ i ] = bulk;
    }
}


template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::computeLinearStiff (dataPtr_Type& /*dataMaterial*/,
                                                                const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                                                const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/)
{
    //! Empty method for neo-hookean material
}


template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::updateJacobianMatrix ( const vector_Type&       disp,
                                                                   const dataPtr_Type&      dataMaterial,
                                                                   const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                   const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                                   const displayerPtr_Type& displayer )
{
    this->M_jacobian.reset (new matrix_Type (*this->M_localMap) );

    displayer->leaderPrint (" \n*********************************\n  ");
    updateNonLinearJacobianTerms (this->M_jacobian, disp, dataMaterial, mapsMarkerVolumes, mapsMarkerIndexes, displayer);
    displayer->leaderPrint (" \n*********************************\n  ");
}


template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::updateNonLinearJacobianTerms ( matrixPtr_Type&       jacobian,
                                                                           const vector_Type&    disp,
                                                                           const dataPtr_Type&   dataMaterial,
                                                                           const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                           const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                                           const displayerPtr_Type&  displayer )
{

    using namespace ExpressionAssembly;

    displayer->leaderPrint ("   Non-Linear S-  updating non linear terms in the Jacobian Matrix (Neo-Hookean)");

    * (jacobian) *= 0.0;

    // mapIterator_Type it;
    // //mapIteratorIndex_Type itIndex;

    // vectorVolumesPtr_Type pointerListOfVolumes;
    // vectorIndexesPtr_Type pointerListOfIndexes;

    // for( it = (*mapsMarkerVolumes).begin(); it != (*mapsMarkerVolumes).end(); it++)
    //  {

    //     //Given the marker pointed by the iterator, let's extract the material parameters
    //     UInt marker = it->first;

    //     // Debug
    //     // UInt markerIndex = itIndex->first;
    //     // ASSERT( marker == markerIndex, "The list of volumes is referring to a marker that is not the same as the marker of index!!!");

    //     pointerListOfVolumes.reset( new vectorVolumes_Type(it->second) );
    //     pointerListOfIndexes.reset( new vectorIndexes_Type( (*mapsMarkerIndexes)[marker]) );

    //     Real mu     = dataMaterial->mu(marker);
    //     Real bulk   = dataMaterial->bulk(marker);

    //Macros to make the assembly more readable
#define deformationGradientTensor ( grad( this->M_dispETFESpace,  disp, this->M_offset) + value(this->M_identity) )
#define detDeformationGradientTensor det( deformationGradientTensor )
#define deformationGradientTensor_T  minusT(deformationGradientTensor)
#define RIGHTCAUCHYGREEN transpose(deformationGradientTensor) * deformationGradientTensor
#define firstInvariantC trace( RIGHTCAUCHYGREEN )
#define firstInvariantCbar pow( detDeformationGradientTensor, (-2.0/3.0) ) * firstInvariantC

    //Assembling Volumetric Part
    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value ( 1.0 / 2.0 ) * parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( value (2.0) *pow (detDeformationGradientTensor, 2.0) - detDeformationGradientTensor + value (1.0) ) * dot ( deformationGradientTensor_T, grad (phi_j) ) * dot ( deformationGradientTensor_T, grad (phi_i) )
              ) >> jacobian;

    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value ( - 1.0 / 2.0 ) * parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow (detDeformationGradientTensor, 2.0) - detDeformationGradientTensor + log (detDeformationGradientTensor) ) * dot ( deformationGradientTensor_T * transpose (grad (phi_j) ) * deformationGradientTensor_T,  grad (phi_i) )
              ) >> jacobian;


    //! ISOCHORIC PART
    //! 1. Stiffness matrix : int { -2/3 * mu * J^(-2/3) *( F^-T : \nabla \delta ) ( F : \nabla \v ) }
    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value (-2.0 / 3.0) * parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow (detDeformationGradientTensor, - (2.0 / 3.0) )  * dot ( deformationGradientTensor_T , grad (phi_j) ) * dot ( deformationGradientTensor , grad (phi_i) )
              ) >> jacobian;


    //! 2. Stiffness matrix : int { 2/9 * mu * ( Ic_iso )( F^-T : \nabla \delta ) ( F^-T : \nabla \v ) }
    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value (2.0 / 9.0) * parameter ( (* (this->M_vectorsParameters) ) [0] ) * firstInvariantCbar  * dot ( deformationGradientTensor_T , grad (phi_j) ) * dot ( deformationGradientTensor_T , grad (phi_i) )
              ) >> jacobian;

    //! 3. Stiffness matrix : int { mu * J^(-2/3) (\nabla \delta : \nabla \v)}
    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow (detDeformationGradientTensor, - (2.0 / 3.0) )  * dot ( grad (phi_j), grad (phi_i) )
              ) >> jacobian;

    //! 4. Stiffness matrix : int { -2/3 * mu * J^(-2/3) ( F : \nabla \delta ) ( F^-T : \nabla \v ) }
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value (-2.0 / 3.0) * parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow (detDeformationGradientTensor, - (2.0 / 3.0) )  * dot ( deformationGradientTensor , grad (phi_j) ) * dot ( deformationGradientTensor_T , grad (phi_i) )
              ) >> jacobian;



    //! 5. Stiffness matrix : int { 1/3 * mu * Ic_iso * (F^-T [\nabla \delta]^t F^-T ) : \nabla \v }
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value (1.0 / 3.0) * parameter ( (* (this->M_vectorsParameters) ) [0] ) * firstInvariantCbar  * dot ( ( deformationGradientTensor_T * transpose (grad (phi_j) ) *deformationGradientTensor_T ), grad (phi_i) )
              ) >> jacobian;

    //    }

    jacobian->globalAssemble();
}


template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::updateJacobianMatrix ( const vector_Type& disp,
		                                                           const vector_Type& pressure,
		                                                           const dataPtr_Type& dataMaterial )
{
    this->M_jacobian.reset (new matrix_Type (*this->M_localMap) );

    this->M_displayer->leaderPrint (" \n*********************************\n  ");
    this->M_displayer->leaderPrint ("   Non-Linear S-  updating non linear terms in the Jacobian Matrix (Mixed Neo-Hookean)");
    this->M_displayer->leaderPrint (" \n*********************************\n  ");


    *(this->M_jacobian) *= 0.0;

    using namespace ExpressionAssembly;

    auto F = grad( this->M_dispETFESpace,  disp, this->M_offset) + value(this->M_identity);
    auto FmT = minusT(F);
    auto C = transpose(F)*F;
    auto p = value(this->M_pressureETFESpace, pressure);

    auto J = det(F);
    auto Jm23= pow(J, -2./3.);
    auto dF= grad (phi_j);
    auto dFmTdF= value (-1.0) * FmT * transpose (dF) * FmT;
    auto dJ = J * FmT;
    auto dJdF = dot(dJ , dF );
    auto d2JdF = dJdF * FmT + J * dFmTdF;
    auto dJm23 = value(-2./3.) * Jm23 * FmT;
    auto dJm23dF = dot( dJm23, dF);
    auto d2Jm23dF = value(-2./3.) * Jm23 * dFmTdF + dJm23dF * FmT;

    auto I1= trace(C);
    auto I1bar = Jm23 * I1;
    auto dI1 = value(2.) * F;
    auto dI1bar = Jm23 * dI1 + I1 *dJm23;
    auto dI1dF = dot( dI1, dF);
    auto dI1bardF = dot( dI1bar, dF);
    auto d2I1dF = value(2.0) * dF;
    auto d2I1bardF = dJm23dF * dI1 + Jm23 * d2I1dF + I1 * d2Jm23dF + dI1dF * dJm23;

    auto mu = parameter ( (* (this->M_vectorsParameters) ) [0] );
//    auto dP = value(0.5) * mu * ( d2I1bardF +   dI1bardF  * dI1bar ) + p * ( J * dFmTdF + dJdF * FmT);

    auto FdF = dot(F, dF);
    auto FmTdF=dot(FmT, dF);

    auto dP = value(0.5 / 3.0) * mu * Jm23 * ( value(3.0) * dF
    		                                 + I1 * FmT * transpose(dF) * FmT
    										 - value(2.0) * FdF * FmT
    										 - value(2.0) * FmTdF * ( F + value(-1.0/3.0) * I1 * FmT ) )
             + p * ( J * dFmTdF + dJdF * FmT);
//    auto dP = value(0.5) * mu * dF + p * ( J * dFmTdF + dJdF * FmT) ; //( d2I1bardF +   dI1bardF  * dI1bar );
    //    	auto dP = grad(phi_j);

    this->M_displayer->leaderPrint ("   Non-Linear Mixed S - updating non linear terms in the Jacobian Matrix (Neo-Hookean)");

    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                dot ( dP , grad (phi_i) )
              ) >> this->M_jacobian;

    this->M_jacobian->globalAssemble();
    this->M_displayer->leaderPrint (" \n*********************************\n  ");
}




template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::apply ( const vector_Type& sol, vector_Type& res,
                                                    const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                    const mapMarkerIndexesPtr_Type mapsMarkerIndexes)
{
    computeStiffness (sol, 0., this->M_dataMaterial, mapsMarkerVolumes, mapsMarkerIndexes, this->M_displayer);
    res += *M_stiff;
}


template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::computeStiffness ( const vector_Type&       disp,
                                                               Real                     /*factor*/,
                                                               const dataPtr_Type&      dataMaterial,
                                                               const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                               const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                               const displayerPtr_Type& displayer )
{
    using namespace ExpressionAssembly;

    this->M_stiff.reset (new vector_Type (*this->M_localMap) );

    displayer->leaderPrint (" \n******************************************************************\n  ");
    displayer->leaderPrint (" Non-Linear S-  Computing the Neo-Hookean nonlinear stiffness vector"     );
    displayer->leaderPrint (" \n******************************************************************\n  ");

    M_stiff.reset (new vector_Type (*this->M_localMap) );
    * (M_stiff) *= 0.0;

    // mapIterator_Type it;
    // //mapIteratorIndex_Type itIndex;

    // vectorVolumesPtr_Type pointerListOfVolumes;
    // vectorIndexesPtr_Type pointerListOfIndexes;

    // for( it = (*mapsMarkerVolumes).begin(); it != (*mapsMarkerVolumes).end(); it++ )
    // {

    //     //Given the marker pointed by the iterator, let's extract the material parameters
    //     UInt marker = it->first;

    //     // Debug
    //     // UInt markerIndex = itIndex->first;
    //     // ASSERT( marker == markerIndex, "The list of volumes is referring to a marker that is not the same as the marker of index!!!");

    //     pointerListOfVolumes.reset( new vectorVolumes_Type(it->second) );
    //     pointerListOfIndexes.reset( new vectorIndexes_Type( (*mapsMarkerIndexes)[marker]) );

    //     Real mu     = dataMaterial->mu(marker);
    //     Real bulk   = dataMaterial->bulk(marker);

    //Computation of the volumetric part
    //
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                value (1.0 / 2.0) * parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow ( detDeformationGradientTensor , 2.0) - detDeformationGradientTensor + log (detDeformationGradientTensor) ) * dot (  deformationGradientTensor_T, grad (phi_i) )
              ) >> M_stiff;

    //Computation of the isochoric part
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow (detDeformationGradientTensor, -2.0 / 3.0) * (dot ( deformationGradientTensor - value (1.0 / 3.0) * firstInvariantC * deformationGradientTensor_T, grad (phi_i) ) )
              ) >> M_stiff;

    //    }
    this->M_stiff->globalAssemble();
}



template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::computeStiffness ( 	const vector_Type& disp,
																const vector_Type& pressure,
																const dataPtr_Type& dataMaterial )
{
    this->M_displayer->leaderPrint (" \n******************************************************************\n  ");
    this->M_displayer->leaderPrint (" Non-Linear S-  Computing the Mixed Neo-Hookean nonlinear stiffness vector"     );
    this->M_displayer->leaderPrint (" \n******************************************************************\n  ");

    M_stiff.reset (new vector_Type (*this->M_localMap) );
    * (M_stiff) *= 0.0;

//    std::cout << "\nMax pressure: " << pressure.normInf();
    using namespace ExpressionAssembly;

    auto F = grad( this->M_dispETFESpace,  disp, this->M_offset) + value(this->M_identity);
    auto FmT = minusT(F);
    auto C = transpose(F)*F;
    auto p = value(this->M_pressureETFESpace, pressure);

    auto J = det(F);
    auto Jm23= pow(J, -2./3.);

    auto I1= trace(C);

    auto mu = parameter ( (* (this->M_vectorsParameters) ) [0] );
    auto P = value(0.5) *mu * Jm23 * ( F + value(-1./3.) * I1 * FmT  ) + J * p * FmT;
//    auto P =  value(0.5) *mu * F + J * p * FmT; //p  * value(this->M_identity);
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
               dot( P, grad(phi_i) )
              ) >> M_stiff;

    this->M_stiff->globalAssemble();
}


template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::showMe ( std::string const& fileNameStiff,
                                                     std::string const& fileNameJacobian)
{
    this->M_stiff->spy (fileNameStiff);
    this->M_jacobian->spy (fileNameJacobian);
}

template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::computeLocalFirstPiolaKirchhoffTensor ( Epetra_SerialDenseMatrix& firstPiola,
        const Epetra_SerialDenseMatrix& tensorF,
        const Epetra_SerialDenseMatrix& cofactorF,
        const std::vector<Real>& invariants,
        const UInt marker)
{

    //Get the material parameters
    Real mu       = this->M_dataMaterial->mu (marker);
    Real bulk     = this->M_dataMaterial->bulk (marker);

    //Computing the first term \muJ^{-2/3}[F-(1/3)tr(C)F^{-T}]
    Epetra_SerialDenseMatrix firstTerm (tensorF);
    Epetra_SerialDenseMatrix copyCofactorF (cofactorF);
    Real scale ( 0.0 );
    scale = -1 * (1.0 / 3.0) * invariants[0];
    copyCofactorF.Scale ( scale );
    firstTerm += copyCofactorF;

    Real coef ( 0.0 );
    coef = mu * std::pow (invariants[3], -2.0 / 3.0);
    firstTerm.Scale ( coef );

    //Computing the second term (volumetric part) J*(bulk/2)(J-1+(1/J)*ln(J))F^{-T}
    Epetra_SerialDenseMatrix secondTerm (cofactorF);
    Real sCoef (0);
    sCoef = invariants[3] * (bulk / 2.0) * (invariants[3] - 1 + (1 / invariants[3]) * std::log (invariants[3]) );
    secondTerm.Scale ( sCoef );

    firstPiola += firstTerm;
    firstPiola += secondTerm;
}



template <typename MeshType>
inline StructuralConstitutiveLaw<MeshType>* createNeoHookeanMaterialNonLinear()
{
    return new NeoHookeanMaterialNonLinear<MeshType >();
}
namespace
{
static bool registerNH = StructuralConstitutiveLaw<LifeV::RegionMesh<LinearTetra> >::StructureMaterialFactory::instance().registerProduct ( "neoHookean", &createNeoHookeanMaterialNonLinear<LifeV::RegionMesh<LinearTetra> > );
}

#undef deformationGradientTensor
#undef detDeformationGradientTensor
#undef deformationGradientTensor_T
#undef RIGHTCAUCHYGREEN
#undef firstInvariantC
#undef firstInvariantCbar

} //Namespace LifeV

#endif /* __NEOHOOKENANMATERIAL_H */
