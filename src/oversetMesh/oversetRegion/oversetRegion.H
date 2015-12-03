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

    Overset region:
    - collects its cells by cell zone
    - contains list of donor regions
    - contains list of patch names used to define a hole (walls)
    - contains the fringe selection algorithm

    Overset region provides data on
    - cells within region (local processor)
    - local acceptor and hole cells
    - eligible local donors: cells that are neither acceptors nor holes

SourceFiles
    oversetRegion.C

\*---------------------------------------------------------------------------*/

#ifndef oversetRegion_H
#define oversetRegion_H

#include "fvMesh.H"
#include "donorAcceptorList.H"
#include "boundBox.H"
#include "triSurfaceMesh.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "optimisationSwitch.H"

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

        //- Overset region name and cell zone name
        const word name_;

        //- Overset region index
        const label index_;

        //- Mesh reference
        const fvMesh& mesh_;

        //- Overset mesh reference
        const oversetMesh& oversetMesh_;

        //- Cell zone index
        const label zoneIndex_;

        //- List of region names that provide donors
        //  Note: ordering is important as it defines precedence
        wordList donorRegionNames_;

        //- Fringe selection algorithm
        autoPtr<oversetFringe> fringePtr_;


    // Demand-driven data

        //- List of region indices that provide donors
        //  Note: ordering is important as it defines precedence
        mutable labelList* donorRegionsPtr_;

        //- Acceptor region labels
        mutable labelList* acceptorRegionsPtr_;

        //- Acceptor cell data
        mutable donorAcceptorList* acceptorCellsPtr_;

        //- Donor cell data.  Donor cell is defined as the nearest cell in
        //  donor regions to the acceptor cell
        mutable donorAcceptorList* donorCellsPtr_;

        //- Hole cell labels
        mutable labelList* holeCellsPtr_;

        //- Eligible donor cells from current region
        //  Eligible donors are region cells that are not holes or acceptors
        mutable labelList* eligibleDonorCellsPtr_;


        // Hole cutting tools

            //- Hole patches as a triangular surface
            mutable triSurface* holeTriMeshPtr_;

            //- Triangular hole surface search algorithm.  Contains complete
            //  cutting surface
            //  Used in hole cutting
            mutable triSurfaceSearch* holeSearchPtr_;


        // Donor search tools

            //- Overset region local (processor) bounding box
            mutable boundBox* localBoundsPtr_;

            //- Overset region global (region) bounding box
            mutable boundBox* globalBoundsPtr_;

            //- Cell octree for search in local region
            //  Used in donor cell search
            mutable indexedOctree<treeDataCell>* cellSearchPtr_;


        //- Donor fill-in.  One in donorFraction cells is expected to be
        //  a donor
        static const debug::optimisationSwitch donorFraction;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        oversetRegion(const oversetRegion&);

        //- Disallow default bitwise assignment
        void operator=(const oversetRegion&);


        //- Prepare hole triangular mesh
        void calcHoleTriMesh() const;

        // Calculate donor regions
        void calcDonorRegions() const;

        // Calculate acceptor regions
        void calcAcceptorRegions() const;

        // Calculate donor and acceptor cells
        void calcDonorAcceptorCells() const;

        // Calculate hole cells
        void calcHoleCells() const;

        // Calculate eligible donor cells
        void calcEligibleDonorCells() const;

        // Calculate bounds
        void calcBounds() const;

        // Calculate cell octree
        void calcCellSearch() const;

        // Clear all data
        void clearOut() const;


public:

    // Constructors

        //- Construct from dictionary
        oversetRegion
        (
            const word& name,
            const label index,
            const fvMesh& mesh,
            const oversetMesh& oversetMesh,
            const dictionary& dict
        );


    //- Destructor
    ~oversetRegion();


    // Member Functions

        // Access

            //- Return region name and cell zone name
            const word& name() const
            {
                return name_;
            }

            //- Return region index
            label index() const
            {
                return index_;
            }

            //- Return mesh reference
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return overset mesh reference
            const oversetMesh& overset() const
            {
                return oversetMesh_;
            }

            //- Return region cell zone index
            label zoneIndex() const
            {
                return zoneIndex_;
            }

            //- Return region cell zone
            const cellZone& zone() const
            {
                return mesh_.cellZones()[zoneIndex()];
            }

            //- Return list of donor region indices
            const labelList& donorRegions() const;

            //- Return list of acceptor region indices
            const labelList& acceptorRegions() const;

            //- Return region cell labels
            const labelList& regionCells() const;


        // Donor-acceptor addressing

            //- Return acceptor cell list
            const donorAcceptorList& acceptors() const;

            //- Return donor cell list
            const donorAcceptorList& donors() const;

            //- Return hole cell list
            const labelList& holes() const;

            //- Return eligible donor cell list
            const labelList& eligibleDonors() const;


        // Hole cutting search functionality

            //- Are hole patches present?
            bool holePatchesPresent() const;

            //- Hole triangular surface
            const triSurface& holeTriMesh() const;

            //- Return triangular surface search object for hole patches
            const triSurfaceSearch& holeSearch() const;


        // Donor cell search functionality

            //- Return local (processor) bounding box for region
            const boundBox& localBounds() const;

            //- Return global (region) bounding box for region
            const boundBox& globalBounds() const;

            //- Return cell octree for search in local region
            //  Only eligible donors are added to the list
            const indexedOctree<treeDataCell>& cellSearch() const;


        //- Update for topo change or mesh motion
        void update() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
