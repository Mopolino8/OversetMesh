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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "overlapFringe.H"
#include "oversetRegion.H"
#include "oversetMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(overlapFringe, 0);
    addToRunTimeSelectionTable(oversetFringe, overlapFringe, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::overlapFringe::calcAddressing() const
{
    if (holesPtr_ || acceptorsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::overlapFringe::calcAddressing() const"
        )   << "Hole-acceptor addressing already calculated"
            << abort(FatalError);
    }

    Info<< "overlapFringe::calcAddressing" << endl;

    // Algorithm
    // - visit all donor regions of the master (= current) region
    // - for each donor region, mark live cells that can be donors,
    //   excluding hole and acceptor cells
    // - for all cells of master region find whether they overlap with the
    //   cells on the other sides and of which type:
    //       - if hole, master cell is also a hole
    //       - if live, master cell is acceptor

    // Get addressing from master region
    const labelList& masterCells = region().regionCells();

    const labelList& donorRegions = region().donorRegions();

    const fvMesh& mesh = region().mesh();

    // Get cell centres
    const vectorField& cc = mesh.cellCentres();

    // Make a hash set to collect acceptor points
    labelHashSet acceptorSet;

    // Make a hash set to collect hole points
    labelHashSet holeSet;

    forAll (donorRegions, regionI)
    {
        // Get reference to current donor region
        const oversetRegion& curDonorRegion =
            region().overset().regions()[donorRegions[regionI]];

            const labelList& curDonors = curDonorRegion.eligibleDonors();

        // Get donor tree
        const indexedOctree<treeDataCell>& tree =
            curDonorRegion.cellSearch();

        // NOTE: PARALLELISATION MISSING
        // Overlap search should be performed on all processors
        // HJ: ToDo

        // If the tree is empty on local processor, do not search
        if (tree.nodes().empty())
        {
            continue;
        }

        scalar span = tree.bb().mag();

        // Go through all master cells and find the nearest cell in the current
        // donor region
        forAll (masterCells, mcI)
        {
            const label& curCell = masterCells[mcI];
            const vector& curCentre = cc[curCell];

            // Find nearest cell in the current master region
            pointIndexHit pih = tree.findNearest(curCentre, span);

            if (pih.hit())
            {
                if (mesh.pointInCell(curCentre, curDonors[pih.index()]))
                {
                    // Found nearest cell in donor region.  This is a
                    // potential acceptor cell
                    acceptorSet.insert(curCell);
                }
            }
        }
    }

    // Collected all overlap donor cells

    // Note: add fringe minimisation here
    // HJ, 20/Jun/2015

    // Collect holes
    holesPtr_ = new labelList(holeSet.sortedToc());

    // Collect acceptors
    acceptorsPtr_ = new labelList(acceptorSet.sortedToc());
}


void Foam::overlapFringe::clearAddressing()
{
    deleteDemandDrivenData(holesPtr_);
    deleteDemandDrivenData(acceptorsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::overlapFringe::overlapFringe
(
    const fvMesh& mesh,
    const oversetRegion& region,
    const dictionary& dict
)
:
    oversetFringe(mesh, region, dict),
    holesPtr_(NULL),
    acceptorsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::overlapFringe::~overlapFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::overlapFringe::holes() const
{
    if (!holesPtr_)
    {
        calcAddressing();
    }

    return *holesPtr_;
}


const Foam::labelList& Foam::overlapFringe::acceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }

    return *acceptorsPtr_;
}


void Foam::overlapFringe::update()
{
    Info<< "overlapFringe::update()" << endl;
    clearAddressing();
}


// ************************************************************************* //
