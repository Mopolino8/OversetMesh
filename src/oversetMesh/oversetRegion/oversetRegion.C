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

\*---------------------------------------------------------------------------*/

#include "oversetRegion.H"
#include "oversetMesh.H"
#include "oversetFringe.H"
#include "polyPatchID.H"
#include "triSurfaceTools.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::debug::optimisationSwitch
Foam::oversetRegion::donorFraction
(
    "donorFraction",
    30
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::oversetRegion::calcDonorRegions() const
{
    if (donorRegionsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcDonorRegions() const")
            << "Donor regions already calculated"
            << abort(FatalError);
    }

    // Get regions
    const PtrList<oversetRegion>& regions = oversetMesh_.regions();

    donorRegionsPtr_ = new labelList(donorRegionNames_.size());
    labelList& donRegions = *donorRegionsPtr_;

    forAll (donorRegionNames_, drI)
    {
        bool found = false;

        // Get name to search for
        const word& curName = donorRegionNames_[drI];

        // Find donor region name
        forAll (regions, orI)
        {
            if (regions[orI].name() == curName)
            {
                // Found donor region name in the list
                found = true;

                donRegions[drI] = orI;
                break;
            }
        }

        if (!found)
        {
            FatalErrorIn("void oversetRegion::calcDonorRegions() const")
                << "For region " << name() << " cannot find donor region "
                << curName << ".  Please check overset definition"
                << abort(FatalError);
        }
    }
}


void Foam::oversetRegion::calcAcceptorRegions() const
{
    if (acceptorRegionsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcAcceptorRegions() const")
            << "Acceptor regions already calculated"
            << abort(FatalError);
    }

    const PtrList<oversetRegion>& regions = oversetMesh_.regions();

    acceptorRegionsPtr_ = new labelList(regions.size());
    labelList& accRegions = *acceptorRegionsPtr_;

    label nAccRegions = 0;

    // Go through all regions apart from the current and check if
    // this region appears in the donor list
    forAll (regions, regionI)
    {
        // Skip current region
        if (regionI == index())
        {
            continue;
        }

        const labelList& remoteDonors = regions[regionI].donorRegions();

        forAll (remoteDonors, rdI)
        {
            if (remoteDonors[rdI] == index())
            {
                // Remote region lists current region as donor
                // Therefore, it is acceptor of this region
                accRegions[nAccRegions] = regionI;
                nAccRegions++;
            }
        }
    }

    // Reset size of acceptor regions
    accRegions.setSize(nAccRegions);
}


void Foam::oversetRegion::calcDonorAcceptorCells() const
{
    if (donorCellsPtr_ || acceptorCellsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcDonorAcceptorCells() const")
            << "Donor cells already calculated"
            << abort(FatalError);
    }

    // Get regions
    const PtrList<oversetRegion>& regions = oversetMesh_.regions();

    // Behave as acceptor: prepare acceptor data for search
    // Operating on acceptor cells for current region on local processor

    // Prepare and broadcast all acceptors from this region

    // Get cell centres
    const vectorField& cc = mesh_.cellCentres();

    // Get list of donor regions
    const labelList& dr = donorRegions();

    // Get list of local acceptor cell labels from fringe
    const labelList& a = fringePtr_().acceptors();

    if (Pstream::parRun())
    {
        // Make a global list of all acceptors
        donorAcceptorListList globalDonorAcceptor(Pstream::nProcs());

        // Memory management
        {
            // Insert local acceptors into a global list
            donorAcceptorList& curDA =
                globalDonorAcceptor[Pstream::myProcNo()];

            curDA.setSize(a.size());

            forAll (a, aI)
            {
                curDA[aI] = donorAcceptor
                (
                    a[aI],
                    Pstream::myProcNo(),
                    cc[a[aI]]
                );
            }
        }

        // Gather-scatter acceptor data before donor indentification
        Pstream::gatherList(globalDonorAcceptor);
        Pstream::scatterList(globalDonorAcceptor);

        // Donor identification: search for donors for all processors
        // using local donor regions

        // Create a list to record local donors.  Guess the size as
        //  donorFraction of number of cells
        DynamicList<donorAcceptor> localDonors
        (
            Foam::max(mesh_.nCells()/donorFraction(), 100)
        );

        // Go through all donor regions and identify donor cells
        forAll (dr, drI)
        {
            if (dr[drI] == index())
            {
                FatalErrorIn
                (
                    "void oversetRegion::calcDonorAcceptorCells() const"
                )   << "Region " << name() << " specified as the donor "
                    << "of itself.  List of donors: " << dr << nl
                    << "This is not allowed: check oversetMesh definition"
                    << abort(FatalError);
            }

            // Get current donor region
            const oversetRegion& curDonorRegion = regions[dr[drI]];

            const labelList& curDonors = regions[dr[drI]].eligibleDonors();

            // Get donor region tree
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
                            // Found a hit.  Additional check for point in cell
                            if
                            (
                                mesh_.pointInCell
//                                 mesh_.pointInCellBB
                                (
                                    curP,
                                    curDonors[pih.index()]
                                )
                            )
                            {
                                // Identified local donor.  Set donor data
                                // to let acceptor know about the match
                                curDA[daI].setDonor
                                (
                                    curDonors[pih.index()],
                                    Pstream::myProcNo()
                                );

                                // Record local donor data.  Note: acceptor
                                // processor may be remote
                                localDonors.append(curDA[daI]);

                                // Note:
                                // If the interpolation stencil requires
                                // larger neighbourhood, additional donor cells
                                // can be located and gathered here.
                                // Issues arise if parts of the interpolation
                                // stencil straddle the processor boundary
                                // HJ, 19/Jun/2015
                            }
                        }
                    }
                }
            }
        }

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
                        if
                        (
                            !recombined[cI].donorFound() 
                         && candidate[cI].donorFound()
                        )
                        {
                            // Candidate has found the donor
                            // Record donor and donor processor
                            recombined[cI].setDonor
                            (
                                candidate[cI].donorCell(),
                                candidate[cI].donorProcNo()
                            );
                        }
                        // Checking
                        else if
                        (
                            recombined[cI].donorFound() 
                         && candidate[cI].donorFound()
                        )
                        {
                            if
                            (
                                (
                                    recombined[cI].donorCell()
                                 != candidate[cI].donorCell()
                                )
                             || (
                                    recombined[cI].donorCell()
                                 != candidate[cI].donorCell()
                                )
                            )
                            {
                                FatalErrorIn
                                (
                                    "void oversetRegion::"
                                    "calcDonorAcceptorCells() const"
                                )   << "Region " << name()
                                    << ": Multiple parallel donor found: "
                                    << recombined[cI] << nl
                                    << candidate[cI]
                                    << abort(FatalError);
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

        // Check if donors have been found for all local acceptors
        {
            // Grab local acceptors
            acceptorCellsPtr_ = new donorAcceptorList
            (
                globalDonorAcceptor[Pstream::myProcNo()]
            );

            donorAcceptorList& curDA = *acceptorCellsPtr_;

            labelList nDonorsFromProc(Pstream::nProcs(), 0);

            label nUncoveredAcceptors = 0;

            forAll (curDA, daI)
            {
                if (!curDA[daI].donorFound())
                {
                    // Donor not found globally
                    WarningIn
                    (
                        "void oversetRegion::calcDonorAcceptorCells() const"
                    )   << "Donor not found for cell "
                        << curDA[daI].acceptorCell() << " on processor "
                        << curDA[daI].acceptorProcNo()
                        << endl;

                    nUncoveredAcceptors++;
                }
                else
                {
                    // Assemble donor statistics
                    nDonorsFromProc[curDA[daI].donorProcNo()]++;
                }
            }

//             Pout<< "Region " << index()
//                 << " number of processor donors for " << curDA.size()
//                 << " local acceptors per processor: " << nDonorsFromProc
//                 << endl;

            // Check for uncovered acceptors
            if (nUncoveredAcceptors > 0)
            {
                FatalErrorIn
                (
                    "void oversetRegion::calcDonorAcceptorCells() const"
                )   << "Inconsistency in donor cell data assembly"
                    << abort(FatalError);
            }

            // Sanity check
            if
            (
                nUncoveredAcceptors == 0
             && sum(nDonorsFromProc) != curDA.size()
            )
            {
                FatalErrorIn
                (
                    "void oversetRegion::calcDonorAcceptorCells() const"
                )   << "Inconsistency in donor cell data assembly"
                    << abort(FatalError);
            }
        }

        donorCellsPtr_ = new donorAcceptorList();
        donorAcceptorList& d = *donorCellsPtr_;

        d.transfer(localDonors.shrink());

        // Sanity check local donors
        labelList nDonorsToProc(Pstream::nProcs(), 0);

        forAll (d, dI)
        {
            // Assemble donor statistics
            nDonorsToProc[d[dI].acceptorProcNo()]++;
        }

//         Pout<< "Region " << index()
//             << " number of local donors = " << d.size()
//             << " per processor: " << nDonorsToProc
//             << endl;
    }
    else
    {
        // Serial run.  All donors and acceptors are local

        // Note: in a serial run, donorCells and acceptorCells
        // contain identical data.  This is duplicated for the moment
        // HJ, 1/May/2015

        acceptorCellsPtr_ = new donorAcceptorList(a.size());
        donorAcceptorList& DA = *acceptorCellsPtr_;

        // Insert local acceptors into the list
        forAll (a, aI)
        {
            DA[aI] = donorAcceptor
            (
                a[aI],
                Pstream::myProcNo(),
                cc[a[aI]]
            );
        }

        // Go through all donor regions and identify donor cells
        forAll (dr, drI)
        {
            if (dr[drI] == index())
            {
                FatalErrorIn
                (
                    "void oversetRegion::calcDonorAcceptorCells() const"
                )   << "Region " << index() << " specified as the donor "
                    << "of itself.  List of donors: " << dr << nl
                    << "This is not allowed: check oversetMesh definition"
                    << abort(FatalError);
            }

            // Get current donor region
            const oversetRegion& curDonorRegion = regions[dr[drI]];

            const labelList& curDonors = regions[dr[drI]].eligibleDonors();

            // Get donor region tree
            const indexedOctree<treeDataCell>& tree =
                curDonorRegion.cellSearch();

            scalar span = tree.bb().mag();

            forAll (DA, daI)
            {
                const vector& curP = DA[daI].acceptorPoint();

                // Find nearest cell with octree
                // Note: octree only contains eligible cells
                // HJ, 10/Jan/2015
                pointIndexHit pih = tree.findNearest(curP, span);

                if (pih.hit())
                {
                    // Found a hit.  Additional check for point in cell
                    if
                    (
                        mesh_.pointInCell
//                         mesh_.pointInCellBB
                        (
                            curP,
                            curDonors[pih.index()]
                        )
                    )
                    {
                        // Identified donor.  Set donor data
                        // to let acceptor know about the match
                        DA[daI].setDonor
                        (
                            curDonors[pih.index()],
                            Pstream::myProcNo()
                        );
                    }
                }
            }
        } // End of all donor regions

        // Check if donors have been found for all local acceptors
        {
            label nUncoveredAcceptors = 0;

            forAll (DA, daI)
            {
                if (!DA[daI].donorFound())
                {
                    // Donor not found =
                    Info<< "Donor not found for cell "
                        << DA[daI].acceptorCell()
                        << endl;

                    nUncoveredAcceptors++;
                }
            }

            // Check for uncovered acceptors
            if (nUncoveredAcceptors > 0)
            {
                FatalErrorIn
                (
                    "void oversetRegion::calcDonorAcceptorCells() const"
                )   << "Inconsistency in donor cell data assembly"
                    << abort(FatalError);
            }
        }

        // Since in serial execution donor and acceptor data is identical
        // copy the acceptor list into donor list after the search and check
        // have been performed
        donorCellsPtr_ = new donorAcceptorList(DA);
    }
}


void Foam::oversetRegion::calcHoleCells() const
{
    if (holeCellsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcHoleCells() const")
            << "Hole cells already calculated"
            << abort(FatalError);
    }

    // Algorithm
    // - go through all regions apart from the current region
    // - get access to hole surface search
    // - mark as hole all cells that fall outside of hole search
    // - combine all searches into a single inside-outside list

    // Get local cell indices
    const labelList& rc = regionCells();

    // Prepare local cell centres for inside-outside search on all regions
    vectorField localC(rc.size());

    const vectorField& c = mesh().cellCentres();

    forAll (localC, i)
    {
        localC[i] = c[rc[i]];
    }

    // Prepare hole mask
    boolList holeMask(mesh().nCells(), false);

    // Mark all cells indicated by fringe
    // Mask additional hole cells as defined by the fringe
    const labelList& fringeHoles = fringePtr_().holes();

    forAll (fringeHoles, i)
    {
        holeMask[fringeHoles[i]] = true;
    }

    // Mark all hole cells using their hole boundary patch inside search

    if (overset().holePatchesPresent())
    {
        // Get reference to hole search
        const triSurfaceSearch& holeSearch = overset().holeSearch();

        boolList regionInside = holeSearch.calcInside(localC);

        // Note: hole mask has the size of all mesh cells and regionInside
        // only of the size of local region
        forAll (regionInside, i)
        {
            holeMask[rc[i]] |= regionInside[i];
        }

        // Count hole cells
        label nHoleCells = 0;

        forAll (rc, i)
        {
            if (holeMask[rc[i]])
            {
                nHoleCells++;
            }
        }

        // Allocate hole cells storage
        holeCellsPtr_ = new labelList(nHoleCells);
        labelList& ch = *holeCellsPtr_;

        // Reset counter and collect hole cells
        nHoleCells = 0;

        forAll (rc, i)
        {
            if (holeMask[rc[i]])
            {
                ch[nHoleCells] = rc[i];
                nHoleCells++;
            }
        }
    }
    else
    {
        // No hole search = no holes
        holeCellsPtr_ = new labelList();
    }
}


void Foam::oversetRegion::calcEligibleDonorCells() const
{
    if (eligibleDonorCellsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcEligibleDonorCells() const")
            << "Hole cells already calculated"
            << abort(FatalError);
    }

    const labelList& rc = regionCells();

    // Prepare donor mask for region cells
    boolList donorMask(mesh().nCells(), false);

    // Mark region cells as eligible
    forAll (rc, rcI)
    {
        donorMask[rc[rcI]] = true;
    }

    // Remove all hole cells from region
    const labelList& hc = holes();

    forAll (hc, hcI)
    {
        donorMask[hc[hcI]] = false;
    }

    // Remove all acceptor cells from fringe
    const labelList& ac = fringePtr_().acceptors();

    forAll (ac, acI)
    {
        donorMask[ac[acI]] = false;
    }

    // Count and collect eligible donors
    label nEligibleDonors = 0;

    forAll (donorMask, cellI)
    {
        if (donorMask[cellI])
        {
            nEligibleDonors++;
        }
    }

    eligibleDonorCellsPtr_ = new labelList(nEligibleDonors);
    labelList& ed = *eligibleDonorCellsPtr_;

    // Reset counter and collect cells
    nEligibleDonors = 0;

    forAll (donorMask, cellI)
    {
        if (donorMask[cellI])
        {
            ed[nEligibleDonors] = cellI;
            nEligibleDonors++;
        }
    }
}


void Foam::oversetRegion::calcCellSearch() const
{
    if (cellSearchPtr_)
    {
        FatalErrorIn("void oversetRegion::calcCellSearch() const")
            << "Cell tree already calculated"
            << abort(FatalError);
    }

    // Create the octree search for this region.  It will be used by other
    // regions when searching for donor cells

    // Reconsider search boxes: only capture local cells
    treeBoundBox overallBb(mesh_.points());
    Random rndGen(123456);
    overallBb = overallBb.extend(rndGen, 1E-4);
    overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    // Search
    cellSearchPtr_ = new indexedOctree<treeDataCell>
    (
        treeDataCell
        (
            false,  //  Cache bb.  Reconsider for moving mesh cases
            mesh_,
            eligibleDonors()
        ),
        overallBb,  // overall search domain
        8,          // maxLevel
        10,         // leafsize
        3.0         // duplicity
    );
}


void Foam::oversetRegion::clearOut()
{
    deleteDemandDrivenData(donorRegionsPtr_);
    deleteDemandDrivenData(acceptorRegionsPtr_);

    deleteDemandDrivenData(acceptorCellsPtr_);
    deleteDemandDrivenData(donorCellsPtr_);
    deleteDemandDrivenData(holeCellsPtr_);
    deleteDemandDrivenData(eligibleDonorCellsPtr_);

    deleteDemandDrivenData(cellSearchPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetRegion::oversetRegion
(
    const word& name,
    const label index,
    const fvMesh& mesh,
    const oversetMesh& oversetMesh,
    const dictionary& dict
)
:
    name_(name),
    index_(index),
    mesh_(mesh),
    oversetMesh_(oversetMesh),
    zoneIndex_(mesh_.cellZones().findZoneID(name_)),
    donorRegionNames_(dict.lookup("donorRegions")),
    fringePtr_(),
    donorRegionsPtr_(NULL),
    acceptorRegionsPtr_(NULL),

    acceptorCellsPtr_(NULL),
    donorCellsPtr_(NULL),
    holeCellsPtr_(NULL),
    eligibleDonorCellsPtr_(NULL),

    cellSearchPtr_(NULL)
{
    // Check zone index
    if (zoneIndex_ < 0)
    {
        FatalErrorIn
        (
            "oversetRegion::oversetRegion\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const oversetMesh& oversetMesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Cannot find cell zone for region " << name << nl
            << "Available cell zones: " << mesh_.cellZones().names()
            << abort(FatalError);
    }

    fringePtr_ = oversetFringe::New
    (
        mesh,
        *this,
        dict.subDict("fringe")
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::oversetRegion::~oversetRegion()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::labelList& Foam::oversetRegion::regionCells() const
{
    return mesh_.cellZones()[zoneIndex_];
}


const Foam::labelList& Foam::oversetRegion::donorRegions() const
{
    if (!donorRegionsPtr_)
    {
        calcDonorRegions();
    }

    return *donorRegionsPtr_;
}


const Foam::labelList& Foam::oversetRegion::acceptorRegions() const
{
    if (!acceptorRegionsPtr_)
    {
        calcAcceptorRegions();
    }

    return *acceptorRegionsPtr_;
}


const Foam::donorAcceptorList& Foam::oversetRegion::acceptors() const
{
    if (!acceptorCellsPtr_)
    {
        calcDonorAcceptorCells();
    }

    return *acceptorCellsPtr_;
}


const Foam::donorAcceptorList& Foam::oversetRegion::donors() const
{
    if (!donorCellsPtr_)
    {
        calcDonorAcceptorCells();
    }

    return *donorCellsPtr_;
}


const Foam::labelList& Foam::oversetRegion::holes() const
{
    if (!holeCellsPtr_)
    {
        calcHoleCells();
    }

    return *holeCellsPtr_;
}


const Foam::labelList& Foam::oversetRegion::eligibleDonors() const
{
    if (!eligibleDonorCellsPtr_)
    {
        calcEligibleDonorCells();
    }

    return *eligibleDonorCellsPtr_;
}

const Foam::indexedOctree<Foam::treeDataCell>&
Foam::oversetRegion::cellSearch() const
{
    if (!cellSearchPtr_)
    {
        calcCellSearch();
    }

    return *cellSearchPtr_;
}


// ************************************************************************* //
