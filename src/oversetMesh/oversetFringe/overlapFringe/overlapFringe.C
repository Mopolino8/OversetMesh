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
    if (acceptorsPtr_)
    {
        FatalErrorIn("void overlapFringe::calcAddressing() const")
            << "Acceptor addressing already calculated"
            << abort(FatalError);
    }

    Info<< "overlapFringe::calcAddressing" << endl;

    // Algorithm
    // - visit all donor regions of the master (= current) region
    // - for each donor region, mark live cells that can be donors,
    //   excluding hole cells
    // - for all cells of master region find whether they overlap with the
    //   cells on the other sides and of which type:
    //       - if hole, master cell is also a hole
    //       - if live, master cell is acceptor

    // Get mesh
    const fvMesh& mesh = region().mesh();

    // Get cell centres
    const vectorField& cc = mesh.cellCentres();

    // Prepare list of potential acceptor cells

    // Get addressing from master region
    const labelList& masterCells = region().regionCells();

    // Get holes in form of mask
    boolList holeMask(region().mesh().nCells(), false);

    // Get and mark holes
    const labelList& masterHoles = region().holes();

    forAll (masterHoles, mhI)
    {
        holeMask[masterHoles[mhI]] = true;
    }

    // Collect all cells eligible to be acceptors
    // Note: cells that act as donors for other regions should be excluded
    // Not sure how to do this now.
    // HJ, 15/Sep/2015

    donorAcceptorList acCand(masterCells.size() - masterHoles.size());

    label nAcCand = 0;

    forAll (masterCells, mcI)
    {
        // If cell is not a hole, add it to candidates
        if (!holeMask[masterCells[mcI]])
        {
            acCand[nAcCand] = donorAcceptor
            (
                masterCells[mcI],
                Pstream::myProcNo(),
                cc[masterCells[mcI]]
            );

            nAcCand++;
        }
    }

    // Perform search to see which cells can find valid donors
    const labelList& dr = region().donorRegions();

    if (Pstream::parRun())
    {
        // Make a global list of all acceptors
        donorAcceptorListList globalDonorAcceptor(Pstream::nProcs());

        // Copy local acceptor list into processor slot
        globalDonorAcceptor[Pstream::myProcNo()] = acCand;

        // Gather-scatter acceptor data before donor indentification
        Pstream::gatherList(globalDonorAcceptor);
        Pstream::scatterList(globalDonorAcceptor);

        // Donor identification: search for donors for all processors
        // using local donor regions
        
        // Go through all donor regions and identify donor cells
        forAll (dr, drI)
        {
            if (dr[drI] == region().index())
            {
                FatalErrorIn("void overlapFringe::calcAddressing() const")
                    << "Region " << region().name()
                    << " specified as the donor of itself.  "
                    << "List of donors: " << dr << nl
                    << "This is not allowed: check oversetMesh definition"
                    << abort(FatalError);
            }

            // Get reference to current donor region
            const oversetRegion& curDonorRegion =
                region().overset().regions()[dr[drI]];

            const labelList& curDonors = curDonorRegion.eligibleDonors();

            // Get donor tree
            const indexedOctree<treeDataCell>& tree =
                curDonorRegion.cellSearch();

            // If the tree is empty on local processor, do not search
            if (tree.nodes().empty())
            {
                continue;
            }

            scalar span = tree.bb().mag();

            // Go through all processors and see if local donor can be found

            // For all acceptors, perform donor search
            // Searching for donor cells on local processors using the
            // requested acceptor data from all processors
            forAll (globalDonorAcceptor, procI)
            {
                List<donorAcceptor>& curDA = globalDonorAcceptor[procI];

                forAll (curDA, daI)
                {
                    if (!curDA[daI].donorFound())
                    {
                        const vector& curP = curDA[daI].acceptorPoint();

                        // Find nearest cell with octree
                        // Note: octree only contains eligible cells
                        // HJ, 10/Jan/2015
                        pointIndexHit pih = tree.findNearest(curP, span);

                        if (pih.hit())
                        {
                            // Note: this check needs to be improved because
                            // pointInCell is not sufficiently robust.
                            // HJ, 16/Sep/2015
                            if
                            (
                                mesh.pointInCell
                                (
                                    curP,
                                    curDonors[pih.index()]
                                )
                            )
                            {
                                // Found nearest cell in donor region.
                                // This cell is a potential acceptor
                                curDA[daI].setDonor
                                (
                                    curDonors[pih.index()],
                                    Pstream::myProcNo(),
                                    cc[curDonors[pih.index()]]
                                );
                            }
                        }
                    }
                }
            } // end for all processor lists
        } // end for all dr

        // Gather-scatter acceptor data after donor search

        // At this point, each processor has filled parts of every other
        // processors's list.  Therefore, a simple gather-scatter will not do
        // Algorithm:
        // - send all processor data to master
        // - master recombines
        // - scatter to all processors
        if (Pstream::master())
        {
            // Receive data from all processors and recombine
            for (label procI = 1; procI < Pstream::nProcs(); procI++)
            {
                // Receive list from slave
                IPstream fromSlave
                (
                    Pstream::blocking,
                    procI
                );

                donorAcceptorListList otherDonorAcceptor(fromSlave);

                // Perform recombination
                // If slave has found the acceptor and recombined list
                // did not, copy the data from the slave into the recombined
                // list
                // If two processors have found the acceptor, report error
                forAll (globalDonorAcceptor, pI)
                {
                    // Get reference to recombined list
                    donorAcceptorList& recombined =
                        globalDonorAcceptor[pI];

                    // Get reference to candidate list
                    const donorAcceptorList& candidate =
                        otherDonorAcceptor[pI];

                    // Compare candidate with recombined list
                    // If candidate has found the donor, record it in the
                    // recombined list
                    forAll (candidate, cI)
                    {
                        if (candidate[cI].donorFound())
                        {
                            if (!recombined[cI].donorFound())
                            {
                                // Candidate has found the donor
                                // Record donor and donor processor
                                recombined[cI].setDonor
                                (
                                    candidate[cI].donorCell(),
                                    candidate[cI].donorProcNo(),
                                    candidate[cI].donorPoint()
                                );
                            }
                        }
                    }
                }
            }
        }
        else
        {
            // Slave processor: send global list to master
            OPstream toMaster
            (
                Pstream::nonBlocking,
                Pstream::masterNo()
            );

            toMaster << globalDonorAcceptor;
        }

        // Scatter recombined list to all processors
        Pstream::scatter(globalDonorAcceptor);

        // Copy local acceptor list from processor slot
        acCand = globalDonorAcceptor[Pstream::myProcNo()];
    }
    else
    {
        // Serial run
        forAll (dr, drI)
        {
            if (dr[drI] == region().index())
            {
                FatalErrorIn("void overlapFringe::calcAddressing() const")
                    << "Region " << region().name()
                    << " specified as the donor of itself.  "
                    << "List of donors: " << dr << nl
                    << "This is not allowed: check oversetMesh definition"
                    << abort(FatalError);
            }

            // Get reference to current donor region
            const oversetRegion& curDonorRegion =
                region().overset().regions()[dr[drI]];

            const labelList& curDonors = curDonorRegion.eligibleDonors();

            // Get donor tree
            const indexedOctree<treeDataCell>& tree =
                curDonorRegion.cellSearch();

            // If the tree is empty on local processor, do not search
            if (tree.nodes().empty())
            {
                continue;
            }

            scalar span = tree.bb().mag();

            // Go through all candidates find the nearest cell in the current
            // donor region
            forAll (acCand, acI)
            {
                if (!acCand[acI].donorFound())
                {
                    const label& curCell = acCand[acI].acceptorCell();
                    const vector& curCentre = cc[curCell];

                    // Find nearest cell in the current master region
                    pointIndexHit pih = tree.findNearest(curCentre, span);

                    if (pih.hit())
                    {
                        // Note: this check needs to be improved because
                        // pointInCell is not sufficiently robust.
                        // HJ, 16/Sep/2015
                        if
                        (
                            mesh.pointInCell(curCentre, curDonors[pih.index()])
                        )
                        {
                            // Found nearest cell in donor region.
                            // This cell is a potential acceptor
                            acCand[acI].setDonor
                            (
                                curDonors[pih.index()],
                                Pstream::myProcNo(),
                                cc[curDonors[pih.index()]]
                            );
                        }
                    }
                }
            }
        } // end for all dr
    }

    // Collected all overlap donor cells

    // Note: add fringe minimisation here
    // HJ, 20/Jun/2015

    // Collect acceptors
    acceptorsPtr_ = new labelList(acCand.size());
    labelList& acc = *acceptorsPtr_;

    label nAcc = 0;

    forAll (acCand, acI)
    {
        if (acCand[acI].donorFound())
        {
            acc[nAcc] = acCand[acI].acceptorCell();

            nAcc++;
        }
    }

    acc.setSize(nAcc);

    Pout<< "Region " <<  region().name() << " found " << nAcc
        << " overlap acceptors"
        << endl;
}


void Foam::overlapFringe::clearAddressing() const
{
    deleteDemandDrivenData(fringeHolesPtr_);
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
    fringeHolesPtr_(NULL),
    acceptorsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::overlapFringe::~overlapFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::overlapFringe::fringeHoles() const
{
    if (!fringeHolesPtr_)
    {
        // fringeHoles currently empty.  HJ, 10/Sep/2015
        fringeHolesPtr_ = new labelList();
    }

    return *fringeHolesPtr_;
}


const Foam::labelList& Foam::overlapFringe::acceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }

    return *acceptorsPtr_;
}


void Foam::overlapFringe::update() const
{
    Info<< "overlapFringe::update() const" << endl;
    clearAddressing();
}


// ************************************************************************* //
