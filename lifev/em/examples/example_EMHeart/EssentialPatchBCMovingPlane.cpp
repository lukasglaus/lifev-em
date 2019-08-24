/*

 * EssentialPatchBCMovingPlane.cpp
 *
 *  Created on: Apr 5, 2019
 *      Author: pamstad
 */

#include "EssentialPatchBCMovingPlane.h"
#include <stdio.h>
#include <lifev/em/examples/example_EMHeart/EssentialPatchBC.hpp>
#include <lifev/em/examples/example_EMHeart/GenericFactory.hpp>
#include <cmath>
#include <vector>
#define PI 3.14159265359

namespace LifeV
{

void EssentialPatchBCMovingPlane::setup(const GetPot& dataFile, const std::string& name)
{
	super::setup(dataFile, name);
	
	for (UInt i (0); i < 3 ;++i )
	{

		starting_point[i] = dataFile( ("solid/boundary_conditions/" + m_Name + "/startingpoint").c_str(), 0.0, i );
		
	}

	for (UInt j (0); j < 3; ++j)
	{

		normal_vector[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/direction").c_str() , 0.0 , j );

	}

	m_maxDisplacement = dataFile ( ("solid/boundary_conditions/" + m_Name + "/displacement").c_str(), 1.0 );

}



const bool EssentialPatchBCMovingPlane::nodeOnPatch(const Vector3D& coord, const Real& time)
{

	bool nodeInArea = false;

	if((normal_vector[0]*coord[0] + normal_vector[1]*coord[1] + normal_vector[2]*coord[2] - normal_vector[0]*starting_point[0]- normal_vector[1]*starting_point[1]-normal_vector[2]*starting_point[2] - m_maxDisplacement) <= 0)
	{
		nodeInArea = true;
	}
	else
	{
		nodeInArea = false;
	}


	return nodeInArea;

}

const bool EssentialPatchBCMovingPlane::nodeOnPatchCurrent(const Vector3D& coord, const Real& time)
{
	bool nodeInArea = 0;
	
	if((normal_vector[0]*coord[0] + normal_vector[1]*coord[1] + normal_vector[2]*coord[2] - normal_vector[0]*starting_point[0]- normal_vector[1]*starting_point[1]-normal_vector[2]*starting_point[2] - activationFunction(time)) <= 0)
	{
		nodeInArea = true;
	}
	else
	{
		nodeInArea = false;
	}

	return nodeInArea;

}


void EssentialPatchBCMovingPlane::modifyPatchArea(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const int& newFlag, const Real& time)
{
	if ( solver.comm()->MyPID() == 0 ) std::cout << "WE ARE IN MODIFY PATCH AREA " << std::endl;

	auto p2FeSpace = solver.electroSolverPtr()->feSpacePtr();
    	auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
    	FESpace<RegionMesh<LinearTetra>, MapEpetra > p1FESpace (p2FeSpace->mesh(), "P1", 1, p2FeSpace->mesh()->comm());

    	//create an epetra vector to set it equal to one where it is part of patch
    	VectorEpetra p1ScalarFieldFaces (p1FESpace.map());

    	p1ScalarFieldFaces *= 0.0;

    	Int p1ScalarFieldFacesDof = p1ScalarFieldFaces.epetraVector().MyLength();
    	int globalIdArray[p1ScalarFieldFacesDof];
    	p1ScalarFieldFaces.blockMap().MyGlobalElements(globalIdArray);

        m_patchFlag = newFlag;
        const auto& mesh = solver.localMeshPtr(); 

	const auto& meshfull = solver.fullMeshPtr();
        unsigned int numNodesOnPatch(0); //here we just initalise an unsigned integer variable
        for (int j(0); j < mesh->numBoundaryFacets(); j++) //returns number of boundary facets
        {
              auto& face = mesh->boundaryFacet(j);
              auto faceFlag = face.markerID();
	                       
              int numPointsOnFace(0);

              for (int k(0); k < 3; ++k) 
              {
	                                   
	              	ID pointGlobalId = face.point(k).id();
                      	auto coord = face.point(k).coordinates();
                      	auto pointInPatch = nodeOnPatchCurrent(coord, time);

                      	if(pointInPatch == true)
                      	{
                      		++numPointsOnFace;
                       		for(int n = 0; n < p1ScalarFieldFacesDof; n++)
                      		{
                       			if(pointGlobalId == globalIdArray[n])
                       			{
                            			p1ScalarFieldFaces[pointGlobalId] = 1.0;

                       			}
                       		}
                       	}

              }
			
        }

            m_patchFacesLocationPtr.reset (new vector_Type (p2FeSpace->map() ));
            *m_patchFacesLocationPtr = p2FeSpace->feToFEInterpolate(p1FESpace, p1ScalarFieldFaces);
 }

vectorPtr_Type EssentialPatchBCMovingPlane::directionalVectorField(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time)
{
	Vector3D current_point_on_plane;
	Real distance;
	auto p2PositionVector = p2PositionVectorInitial(dFeSpace, solver);

        vectorPtr_Type p2PatchDisplacement (new VectorEpetra( dFeSpace->map(), Repeated ));
        auto nCompLocalDof = p2PatchDisplacement->epetraVector().MyLength() / 3;

        direction.normalize(); 

        if(normal_vector[2] != 0.0)
        {
	       	current_point_on_plane[2] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2] +activationFunction(time))/normal_vector[2];
	       	current_point_on_plane[1] = 0;
	       	current_point_on_plane[0] = 0;

	}
	else if (normal_vector[2] == 0.0 && normal_vector[1] != 0.0)
	{
	       	current_point_on_plane[1] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2]  + activationFunction(time))/normal_vector[1];
	       	current_point_on_plane[2] = 0;
	       	current_point_on_plane[0] = 0;
        }
        else if (normal_vector[2] == 0.0 && normal_vector[1] == 0.0 && normal_vector[0] != 0.0)
        {
        	current_point_on_plane[0] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2]  +activationFunction(time))/normal_vector[0];
        	current_point_on_plane[1] = 0;
        	current_point_on_plane[2] = 0;
        }
        else
        {	
        	std::cout << "A normal  vector in the data file of (0, 0 , 0) doesn't make sense" << std::endl;
        }

        for (int j (0); j < nCompLocalDof; ++j)
        {
                UInt iGID = p2PatchDisplacement->blockMap().GID (j);
                UInt jGID = p2PatchDisplacement->blockMap().GID (j + nCompLocalDof);
                UInt kGID = p2PatchDisplacement->blockMap().GID (j + 2 * nCompLocalDof);
      
                Vector3D coordinates;

                coordinates(0) = p2PositionVector[iGID];
                coordinates(1) = p2PositionVector[jGID];
                coordinates(2) = p2PositionVector[kGID];
		
                Vector3D QP; //define here the vector that goes from Q (point on plane) to point P

                QP = coordinates - current_point_on_plane;

                if(QP.dot(direction) <= 0) 
		{
                 	distance = abs(QP.dot(direction));
                }
                else
                {
                   	distance = 0.0;
                }

                Vector3D displacement_vector;
                displacement_vector[0] = distance*direction[0];
                displacement_vector[1] = distance*direction[1];
                displacement_vector[2] = distance*direction[2];

				

                (*p2PatchDisplacement)[iGID] = displacement_vector[0];
                (*p2PatchDisplacement)[jGID] = displacement_vector[1];
                (*p2PatchDisplacement)[kGID] = displacement_vector[2];

	        }
	
        return p2PatchDisplacement;

}



//REGISTER(EssentialPatchBC, EssentialPatchBCMovingPlane);

}//this is Klammer von LifeV namespace
