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

    // Allocate storage for holes and acceptors
    holesPtr_ = new labelList(masterCells.size());
    labelList& masterHoles = *holesPtr_;
    label nMasterHoles = 0;

    acceptorsPtr_ = new labelList(masterCells.size());
    labelList& masterAcceptors = *acceptorsPtr_;
    label nMasterAcceptors = 0;

    forAll (donorRegions, regionI)
    {
        // Get reference to current donor region
        const oversetRegion& curDonorRegion =
            region().overset().region(donorRegions[regionI]);

        // Mark eligible cells from current donor region
        const labelList& rc = curDonorRegion.regionCells();

        // Prepare donor mask for donor region cells
        boolList donorMask(mesh.nCells(), false);

        // Prepare hole mask for donor region cells
        boolList holeMask(mesh.nCells(), false);

        // Mark region cells as eligible
        forAll (rc, rcI)
        {
            donorMask[rc[rcI]] = true;
        }

        // Remove all hole cells from region
        const labelList& hc = curDonorRegion.holes();

        forAll (hc, hcI)
        {
            donorMask[hc[hcI]] = false;
            holeMask[hc[hcI]] = true;
        }

        // Remove all acceptor cells from region
        const labelList& ac = curDonorRegion.acceptors();

        forAll (ac, acI)
        {
            donorMask[ac[acI]] = false;
        }

        // End of preparation of current donor region

        // Get donor search
//         const holeTriSurfSearch& donorTree = curDonorRegion.cellTree();

        // Go through all master cells and find the nearest cell in the current
        // donor region
        forAll (masterCells, mcI)
        {
            const label& curCell = masterCells[mcI];
            const vector& curCentre = cc[curCell];

            // Find nearest hole cell in the current master region

                // N-squared search - testing

            scalar minDistance = GREAT;
            label nearest = -1;
            scalar magSqrDist;

            forAll (rc, rcI)
            {
                magSqrDist = magSqr(cc[rc[rcI]] - curCentre);

                if (magSqrDist < minDistance)
                {
                    nearest = rc[rcI];
                    minDistance = magSqrDist;
                }
            }

            // Check for point in cell
            if (nearest > -1)
            {
//                 if (mesh.pointInCell(curCentre, nearest))
                if (mesh.pointInCellBB(curCentre, nearest))
                {
                    // Found nearest cell in donor region

                    if (holeMask[nearest])
                    {
                        // If nearest cell is a hole, master cell is also
                        // a hole
                        masterHoles[nMasterHoles] = curCell;
                        nMasterHoles++;
                    }
                    else if (donorMask[nearest])
                    {
                        // If nearest cell is a donor, master cell is
                        // acceptor
                        masterAcceptors[nMasterAcceptors] = curCell;
                        nMasterAcceptors++;
                    }
                }
            }
        }
    }

    // Resize hole and acceptor lists
    masterAcceptors.setSize(nMasterAcceptors);
    masterHoles.setSize(nMasterHoles);
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
}


// ************************************************************************* //
