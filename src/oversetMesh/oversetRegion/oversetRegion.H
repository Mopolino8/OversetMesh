/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::oversetRegion

Description
    Overset mesh region

SourceFiles
    oversetRegion.C

\*---------------------------------------------------------------------------*/

#ifndef oversetRegion_H
#define oversetRegion_H

#include "fvMesh.H"
#include "triSurfaceMesh.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "indexedOctree.H"
#include "treeDataCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class oversetMesh;
class oversetFringe;

/*---------------------------------------------------------------------------*\
                      Class oversetRegion Declaration
\*---------------------------------------------------------------------------*/

class oversetRegion
{
    // Private data

        //- Mesh reference
        const fvMesh& mesh_;

        //- Overset mesh reference
        const oversetMesh& oversetMesh_;

        //- Value of region split index for this region
        label index_;

        //- List of region indices that provide donors
        //  Note: ordering is important as it defines precedence
        labelList donorRegions_;

        //- Hole patch names
        wordList holePatchNames_;

        //- Fringe selection algorithm
        autoPtr<oversetFringe> fringePtr_;


    // Demand-driven data

        //- Region cell labels
        mutable labelList* regionCellsPtr_;

        //- Acceptor region labels
        mutable labelList* acceptorRegionsPtr_;

        //- Acceptor cell labels
        mutable labelList* acceptorCellsPtr_;

        //- Donor cell labels.  Donor cell is defined as the nearest cell in
        //  donor regions to the acceptor cell
        mutable labelList* donorCellsPtr_;

        //- Hole cell labels
        mutable labelList* holeCellsPtr_;


        // Search tools

            //- Hole patches as a triangular surface
            mutable triSurfaceMesh* holeTriMeshPtr_;

            //- Triangular surface search algorithm
            mutable triSurfaceSearch* triSurfSearchPtr_;

            //- Cell octree for search in local region
            mutable indexedOctree<treeDataCell>* cellTreePtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        oversetRegion(const oversetRegion&);

        //- Disallow default bitwise assignment
        void operator=(const oversetRegion&);


        // Calculate region cells
        void calcRegionCells() const;

        // Calculate acceptor regions
        void calcAcceptorRegions() const;

        // Calculate acceptor cells
        void calcAcceptorCells() const;

        // Calculate donor cells
        void calcDonorCells() const;

        // Calculate hole cells
        void calcHoleCells() const;

        // Calculate cell octree
        void calcCellTree() const;

        // Clear all data
        void clearOut();


public:

    // Constructors

        //- Construct from dictionary
        oversetRegion
        (
            const fvMesh& mesh,
            const oversetMesh& oversetMesh,
            const dictionary& dict
        );


    //- Destructor
    ~oversetRegion();


    // Member Functions

        // Access

            //- Return mesh reference
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return region index as given by regionSplit
            label index() const
            {
                return index_;
            }

            //- Return list of donor region indices
            const labelList& donorRegions() const
            {
                return donorRegions_;
            }

            //- Return list of acceptor region indices
            const labelList& acceptorRegions() const;

            //- Return region cell labels
            const labelList& regionCells() const;


        // Donor-acceptor addressing

            //- Return acceptor cell list
            const labelList& acceptors() const;

            //- Return  donor cell list
            const labelList& donors() const;

            //- Return hole cell list
            const labelList& holes() const;


        // Search functionality

            //- Are hole patches present?
            bool holePatches() const
            {
                return !holePatchNames_.empty();
            }

            //- Hole triangular surface
            const triSurfaceMesh& holeTriMesh() const;

            //- Return triangular surface search object for hole patches
            const triSurfaceSearch& holeTriSurfSearch() const;

            //- Return cell octree for search in local region
            const indexedOctree<treeDataCell>& cellTree() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
