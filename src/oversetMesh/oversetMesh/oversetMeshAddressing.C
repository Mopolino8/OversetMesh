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

#include "oversetMesh.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "polyPatchID.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::oversetMesh::calcCellClassification() const
{
    if (acceptorCellsPtr_ || donorCellsPtr_ || holeCellsPtr_)
    {
        FatalErrorIn("void oversetMesh::calcCellClassification() const")
            << "Cell classification already calculated"
            << abort(FatalError);
    }

    // Count acceptor and donor and hole cells
    label nAcceptorCells = 0;
    label nDonorCells = 0;
    label nHoleCells = 0;

    forAll (regions_, regionI)
    {
        nAcceptorCells += regions_[regionI].acceptors().size();
        nDonorCells += regions_[regionI].donors().size();
        nHoleCells += regions_[regionI].holes().size();
    }

    Pout<< "Number of acceptor cells: " << nAcceptorCells << endl;
    acceptorCellsPtr_ = new labelList(nAcceptorCells);
    labelList& acceptor = *acceptorCellsPtr_;

    Pout<< "Number of donor cells: " << nDonorCells << endl;
    donorCellsPtr_ = new labelList(nDonorCells);
    labelList& donor = *donorCellsPtr_;

    Pout<< "Number of hole cells: " << nHoleCells << endl;
    holeCellsPtr_ = new labelList(nHoleCells);
    labelList& hole = *holeCellsPtr_;

    // Reset counters
    nAcceptorCells = 0;
    nDonorCells = 0;
    nHoleCells = 0;

    forAll (regions_, regionI)
    {
        // Acceptors
        const donorAcceptorList& curAcceptors = regions_[regionI].acceptors();

        forAll (curAcceptors, aI)
        {
            acceptor[nAcceptorCells] = curAcceptors[aI].acceptorCell();
            nAcceptorCells++;
        }

        // Donors
        const donorAcceptorList& curDonors = regions_[regionI].donors();

        forAll (curDonors, dI)
        {
            donor[nDonorCells] = curDonors[dI].donorCell();
            nDonorCells++;
        }

        // Holes
        const labelList& curHoles = regions_[regionI].holes();

        forAll (curHoles, hI)
        {
            hole[nHoleCells] = curHoles[hI];
            nHoleCells++;
        }
    }

    // Check donor and acceptor assembly
    Info<< "Checking donor-acceptor assembly" << endl;

    // Check for multiple acceptors

    // Check for donors that are holes

    // Check for donors that are acceptors
}


void Foam::oversetMesh::calcDomainMarkup() const
{
    if (oversetTypesPtr_ || regionIDPtr_)
    {
        FatalErrorIn("void oversetMesh::calcDomainMarkup() const")
            << "Domain markup already calculated"
            << abort(FatalError);
    }

    // Overset types
    oversetTypesPtr_ = new volScalarField
    (
        IOobject
        (
            "oversetTypes",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, ACTIVE)
    );
    volScalarField& overTypes = *oversetTypesPtr_;

    scalarField& overTypesIn = overTypes.internalField();

    // Mark donor cells
    const labelList& d = donorCells();

    forAll (d, cellI)
    {
        overTypesIn[d[cellI]] = DONOR;
    }

    // Mark acceptor cells
    const labelList& a = acceptorCells();

    forAll (a, cellI)
    {
        overTypesIn[a[cellI]] = ACCEPTOR;
    }

    // Mark all hole cells
    const labelList& h = holeCells();

    forAll (h, cellI)
    {
        overTypesIn[h[cellI]] = HOLE;
    }


    // Update boundary values
    forAll (overTypes.boundaryField(), patchI)
    {
        if (!overTypes.boundaryField()[patchI].coupled())
        {
            overTypes.boundaryField()[patchI] =
                overTypes.boundaryField()[patchI].patchInternalField();
        }
    }

    // Region ID
    regionIDPtr_ = new labelList(mesh().nCells(), -1);
    labelList& rID = *regionIDPtr_;

    // Mark regions

    forAll (regions_, regionI)
    {
        const labelList& curCells = regions_[regionI].zone();

        forAll (curCells, curCellI)
        {
            rID[curCells[curCellI]] = regionI;
        }
    }

    // Check regions
    if (min(rID) < 0)
    {
        FatalErrorIn("void oversetMesh::calcDomainMarkup() const")
            << "Found cells without region ID.  Please check overset setup"
            << abort(FatalError);
    }
}


void Foam::oversetMesh::calcGamma() const
{
    if (gammaPtr_ || gammaExtPtr_ || sGammaPtr_)
    {
        FatalErrorIn("void oversetMesh::calcGamma() const")
            << "Markup fields already calculated"
            << abort(FatalError);
    }

    // Trigger the calculation of cell centres to avoid tangled parallel
    // communications in case the lazy evaluation mechanism is invoked in
    // GeometricField::evaluate() member function. Temporary solution.
    // To-do: examine stack trace for fvMesh::makeC(), VV, 8/Feb/2015.
    mesh().C();

    // Fluid cells indicator, marking only live cells
    gammaPtr_ = new volScalarField
    (
        IOobject
        (
            "oversetGamma",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("one", dimless, 1)
    );
    volScalarField& g = *gammaPtr_;
    scalarField& gIn = g.internalField();

    // Non-hole cells indicator, marking live and acceptor cells
    gammaExtPtr_ = new volScalarField
    (
        IOobject
        (
            "oversetGammaExt",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("one", dimless, 1)
    );
    volScalarField& gExt = *gammaExtPtr_;
    scalarField& gExtIn = gExt.internalField();

    // Mark hole cells
    const labelList& holes = holeCells();

    forAll (holes, holeI)
    {
        gIn[holes[holeI]] = 0;
        gExtIn[holes[holeI]] = 0;
    }

    // Mark acceptor cells: they are marked with 1 in gammaExt, but with
    // 0 in gamma
    const labelList& acceptors = acceptorCells();

    forAll (acceptors, accI)
    {
        gIn[acceptors[accI]] = 0;
    }

    // Before collecting active faces, update coupled boundaries
    // to make sure patchNeighbourField is correct

    // Not allowed to call correctBoundaryConditions.  HJ, 16/Apr/2012
    // Evaluate coupled boundaries and copy out the uncoupled ones

    g.boundaryField().evaluateCoupled();

    forAll (g.boundaryField(), patchI)
    {
        if (!g.boundaryField()[patchI].coupled())
        {
            g.boundaryField()[patchI] =
                g.boundaryField()[patchI].patchInternalField();
        }
    }

    gExt.boundaryField().evaluateCoupled();

    forAll (gExt.boundaryField(), patchI)
    {
        if (!gExt.boundaryField()[patchI].coupled())
        {
            gExt.boundaryField()[patchI] =
                gExt.boundaryField()[patchI].patchInternalField();
        }
    }

    // Calculate face mask
    sGammaPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "oversetSGamma",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("one", dimless, 0)
    );
    surfaceScalarField& sg = *sGammaPtr_;
    scalarField& sgIn = sg.internalField();

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    // Internal faces: flux is live between all active and acceptor cells
    forAll (sgIn, faceI)
    {
        // If both cells are live, flux is live
        if
        (
            gExtIn[owner[faceI]] > SMALL
         && gExtIn[neighbour[faceI]] > SMALL
        )
        {
            sgIn[faceI] = 1;
        }
    }

    // Remove faces where both owner and neighbour are acceptors
    volScalarField gAcc = gExt - g;
    gAcc.boundaryField().evaluateCoupled();

    const scalarField& gAccIn = gAcc.internalField();

    // Kill fluxes between two acceptor cells
    forAll (sgIn, faceI)
    {
        // If both cells are live, flux is live
        if
        (
            gAccIn[owner[faceI]] > SMALL
         && gAccIn[neighbour[faceI]] > SMALL
        )
        {
            sgIn[faceI] = 0;
        }
    }

    surfaceScalarField::GeometricBoundaryField& sgPatches =
        sg.boundaryField();

    const volScalarField::GeometricBoundaryField& gExtPatches =
        gExt.boundaryField();

    const volScalarField::GeometricBoundaryField& gAccPatches =
        gAcc.boundaryField();


    forAll (gExtPatches, patchI)
    {
        if (gExtPatches[patchI].coupled())
        {
            scalarField& gP = sgPatches[patchI];

            // For coupled patches, check gammaExt
            const scalarField gammaOwn =
                gExtPatches[patchI].patchInternalField();

            const scalarField gammaNei =
                gExtPatches[patchI].patchNeighbourField();

            forAll (gammaOwn, faceI)
            {
                if
                (
                    gammaOwn[faceI] > SMALL
                 && gammaNei[faceI] > SMALL
                )
                {
                    gP[faceI] = 1;
                }
            }

            const scalarField gammaAccOwn =
                gAccPatches[patchI].patchInternalField();

            const scalarField gammaAccNei =
                gAccPatches[patchI].patchNeighbourField();

            forAll (gammaOwn, faceI)
            {
                if
                (
                    gammaAccOwn[faceI] > SMALL
                 && gammaAccNei[faceI] > SMALL
                )
                {
                    gP[faceI] = 0;
                }
            }
        }
        else
        {
            // For regular patches, check live cells only to achieve
            // correct global mass adjustment.
            // HJ, 21/May/2012
            const scalarField gammaFc =
                g.boundaryField()[patchI].patchInternalField();

            scalarField& gP = sgPatches[patchI];

            forAll (gammaFc, faceI)
            {
                if (gammaFc[faceI] > SMALL)
                {
                   gP[faceI] = 1;
                }
            }
        }
    }
}


void Foam::oversetMesh::calcFringeFaces() const
{
    if
    (
        fringeFacesPtr_
     || fringeFaceCellsPtr_
     || fringeFaceFlipsPtr_
     || acceptorInternalFacesPtr_
    )
    {
        FatalErrorIn("void oversetMesh::calcFringeFaces() const")
            << "Face masks already calculated"
            << abort(FatalError);
    }

    labelList acCellIndicator(mesh().nCells(), -1);

    // Mark hole cells with -2
    const labelList& hc = holeCells();

    forAll (hc, hcI)
    {
        acCellIndicator[hc[hcI]] = -2;
    }

    // Mark acceptor cells with their index
    const labelList& acc = acceptorCells();

    forAll (acc, accI)
    {
        acCellIndicator[acc[accI]] = accI;
    }

    DynamicList<label> acF(2*acc.size());
    DynamicList<label> acFC(2*acc.size());
    DynamicList<bool> acFF(2*acc.size());

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    const volScalarField::GeometricBoundaryField& gammaPatches =
        gamma().boundaryField();

    // Fringe face collection

    forAll (neighbour, faceI)
    {
        // Old check done using mag gamma.  Changed for internal cells
        if
        (
            acCellIndicator[owner[faceI]] == -1
         && acCellIndicator[neighbour[faceI]] > -1
        )
        {
            // Owner is live, neighbour AC.  Its AC index is in
            // acCellIndicator
            acF.append(faceI);
            acFC.append(acCellIndicator[neighbour[faceI]]);
            acFF.append(false);
        }
        else if
        (
            acCellIndicator[owner[faceI]] > -1
         && acCellIndicator[neighbour[faceI]] == -1
        )
        {
            // Neighbour is live, owner AC.  Its AC index is in
            // acCellIndicator
            acF.append(faceI);
            acFC.append(acCellIndicator[owner[faceI]]);
            acFF.append(true);
        }
    }

    forAll (gammaPatches, patchI)
    {
        // Note: take faceCells from fvPatch (because of empty)
        const labelList& fc = mesh().boundary()[patchI].faceCells();
        const label start = mesh().boundaryMesh()[patchI].start();

        if (gammaPatches[patchI].coupled())
        {
            const scalarField gammaOwn =
                gammaPatches[patchI].patchInternalField();

            const scalarField gammaNei =
                gammaPatches[patchI].patchNeighbourField();

            forAll (gammaOwn, patchFaceI)
            {
                // Note: we assume that there are no live cells neighbouring
                // hole cells, only acceptors. This is forced in overset fringe
                // algorithms (overlapFringe specifically).
                if
                (
                    mag(gammaNei[patchFaceI] - gammaOwn[patchFaceI]) > SMALL
                )
                {
                    if (acCellIndicator[fc[patchFaceI]] > -1)
                    {
                        // Owner cell is AC
                        acF.append(start + patchFaceI);
                        acFC.append(acCellIndicator[fc[patchFaceI]]);
                        acFF.append(true);
                    }
                    else
                    {
                        // Neighbour cell is AC
                        acF.append(start + patchFaceI);
                        acFC.append(-1);
                        acFF.append(false);
                    }
                }
            }
        }
    }

    // Pack the data
    acF.shrink();
    acFC.shrink();
    acFF.shrink();

    fringeFacesPtr_ = new labelList(acF.xfer());

    fringeFaceCellsPtr_ = new labelList(acFC.xfer());
    fringeFaceFlipsPtr_ = new boolList(acFF.xfer());

    // Acceptor internal face collection

    DynamicList<label> acInternalF(2*acc.size());

    forAll (neighbour, faceI)
    {
        // Old check done using mag gamma.  Changed for internal cells
        if
        (
            acCellIndicator[owner[faceI]] > -1
         && acCellIndicator[neighbour[faceI]] > -1
        )
        {
            // Found face between two acceptor cells
            acInternalF.append(faceI);
        }
    }

    const volScalarField::GeometricBoundaryField& gammaExtPatches =
        gammaExt().boundaryField();

    forAll (gammaPatches, patchI)
    {
        // Note: take faceCells from fvPatch (because of empty)
        const label start = mesh().boundaryMesh()[patchI].start();

        if (gammaPatches[patchI].coupled())
        {
            // Create acceptor mask
            const scalarField gammaAccOwn =
                gammaExtPatches[patchI].patchInternalField()
              - gammaPatches[patchI].patchInternalField();

            const scalarField gammaAccNei =
                gammaExtPatches[patchI].patchNeighbourField()
              - gammaPatches[patchI].patchNeighbourField();

            forAll (gammaAccOwn, patchFaceI)
            {
                if
                (
                    mag(gammaAccNei[patchFaceI]) > SMALL
                 && mag(gammaAccOwn[patchFaceI]) > SMALL
                )
                {
                    acInternalF.append(start + patchFaceI);
                }
            }
        }
    }

    // Pack the data
    acInternalF.shrink();

    acceptorInternalFacesPtr_ = new labelList(acInternalF.xfer());
}


void Foam::oversetMesh::calcHoleFaces() const
{
    if
    (
        holeFacesPtr_
     || holeFaceCellsPtr_
     || holeFaceFlipsPtr_
     || holeInternalFacesPtr_
    )
    {
        FatalErrorIn("void oversetMesh::calcHoleFaces() const")
            << "Face masks already calculated"
            << abort(FatalError);
    }

    labelList hole(mesh().nCells(), -1);

    // Mark hole cells with their index
    const labelList& hc = holeCells();

    forAll (hc, hcI)
    {
        hole[hc[hcI]] = hcI;
    }

    DynamicList<label> hF(2*hc.size());
    DynamicList<label> hFC(2*hc.size());
    DynamicList<bool> hFF(2*hc.size());

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    forAll (neighbour, faceI)
    {
        // Old check done using mag gamma.  Changed for internal cells
        if
        (
            hole[owner[faceI]] == -1
         && hole[neighbour[faceI]] > -1
        )
        {
            // Owner is live (acceptor), neighbour H.  Its H index is in
            // hole
            hF.append(faceI);
            hFC.append(hole[neighbour[faceI]]);
            hFF.append(false);
        }
        else if
        (
            hole[owner[faceI]] > -1
         && hole[neighbour[faceI]] == -1
        )
        {
            // Neighbour is live (acceptor), owner H.  Its H index is in
            // hole
            hF.append(faceI);
            hFC.append(hole[owner[faceI]]);
            hFF.append(true);
        }
    }

    // Bug fix: need to use gammaExt instead of gamma to mark hole faces
    const volScalarField::GeometricBoundaryField& gammaExtPatches =
        gammaExt().boundaryField();

    forAll (gammaExtPatches, patchI)
    {
        // Note: take faceCells from fvPatch (because of empty)
        const labelList& fc = mesh().boundary()[patchI].faceCells();
        const label start = mesh().boundaryMesh()[patchI].start();

        if (gammaExtPatches[patchI].coupled())
        {
            const scalarField gammaExtOwn =
                gammaExtPatches[patchI].patchInternalField();

            const scalarField gammaExtNei =
                gammaExtPatches[patchI].patchNeighbourField();

            forAll (gammaExtOwn, patchFaceI)
            {
                // Note: we assume that there are no live cells neighbouring
                // hole cells, only acceptors. This is forced in overset fringe
                // algorithms (overlapFringe specifically).
                if
                (
                    mag(gammaExtNei[patchFaceI] - gammaExtOwn[patchFaceI])
                  > SMALL
                )
                {
                    if (hole[fc[patchFaceI]] > -1)
                    {
                        // Owner cell is H
                        hF.append(start + patchFaceI);
                        hFC.append(hole[fc[patchFaceI]]);
                        hFF.append(true);
                    }
                    else
                    {
                        // Neighbour cell is H
                        hF.append(start + patchFaceI);
                        hFC.append(-1);
                        hFF.append(false);
                    }
                }
            }
        }
    }

    // Pack the data
    hF.shrink();
    hFC.shrink();
    hFF.shrink();

    holeFacesPtr_ = new labelList(hF.xfer());

    holeFaceCellsPtr_ = new labelList(hFC.xfer());
    holeFaceFlipsPtr_ = new boolList(hFF.xfer());


    // Hole internal faces

    // Hole faces are the ones where both owner and neighbour is a hole

    DynamicList<label> hIF(2*hc.size());

    forAll (neighbour, faceI)
    {
        // Old check done using mag gamma.  Changed for internal cells
        if
        (
            hole[owner[faceI]] > -1
         && hole[neighbour[faceI]] > -1
        )
        {
            hIF.append(faceI);
        }
    }

    forAll (gammaExtPatches, patchI)
    {
        // Note: take faceCells from fvPatch (because of empty)
        const label start = mesh().boundaryMesh()[patchI].start();

        if (gammaExtPatches[patchI].coupled())
        {
            const scalarField gammaExtOwn =
                gammaExtPatches[patchI].patchInternalField();

            const scalarField gammaExtNei =
                gammaExtPatches[patchI].patchNeighbourField();

            forAll (gammaExtOwn, patchFaceI)
            {
                if
                (
                    gammaExtNei[patchFaceI] < SMALL
                 && gammaExtOwn[patchFaceI] < SMALL
                )
                {
                    hIF.append(start + patchFaceI);
                }
            }
        }
    }

    // Pack the data
    hIF.shrink();
    holeInternalFacesPtr_ = new labelList(hIF.xfer());
}


void Foam::oversetMesh::calcInterpolationAddressing() const
{
    if
    (
        mapPtr_
     || nDonorsToProcessorMapPtr_
     || donorsToAcceptorsMapPtr_
    )
    {
        FatalErrorIn("void oversetMesh::calcInterpolationAddressing() const")
            << "Fringe addressing already calculated"
            << abort(FatalError);
    }

    // Create donor mask excluding acceptor and hole cells from the extended
    // (neighbourhood) donor selection
    boolList donorMask(mesh().nCells(), true);

    // Mask hole cells
    const labelList& hc = holeCells();

    forAll (hc, hcI)
    {
        donorMask[hc[hcI]] = false;
    }

    // Mark the acceptors with acceptor cell index and exclude them from the
    // extended donor selection
    const labelList& ac = acceptorCells();

    labelList acIndex(mesh().nCells(), -1);

    forAll (ac, acI)
    {
        acIndex[ac[acI]] = acI;

        // Acceptor cells cannot be extended (neighbouring) donors
        donorMask[ac[acI]] = false;
    }

    // Mark the donors with donor cell index
    const labelList& dc = donorCells();

    labelList dcIndex(mesh().nCells(), -1);

    forAll (dc, dcI)
    {
        dcIndex[dc[dcI]] = dcI;
    }

    // Create list containing local donors
    localDonorsPtr_ = new labelList(ac.size(), -1);
    labelList& locDonors = *localDonorsPtr_;

    // Create local neighbouring donors
    localNeighbouringDonorsPtr_ = new labelListList(ac.size());
    labelListList& locNeiDonors = *localNeighbouringDonorsPtr_;

    // Get cell-cell addressing for neighbouring donor search
    const labelListList& cc = mesh().cellCells();

    // Create local donor to local acceptor addressing
    localDonorAddrPtr_ = new labelList(ac.size(), -1);
    labelList& locDonorAddr = *localDonorAddrPtr_;

    if (Pstream::parRun())
    {
        // Parallel run: create mapDistribute
        // Algorithm:
        // 1) First calculate sending map
        //    - Loop through all regions, marking all donors
        //    - For each donor, we must keep track of the processor we are
        //      sending donor data to. Note that a single donor data can go to
        //      multiple processors.
        //    * Note: We currently do not take into account the possibility of
        //      having a donor from the neighbourhood of master donor being on
        //      another processor, i.e. all donors for acceptor on processor I
        //      have to come from the same processor J (where I = J means that
        //      no communication is required). This will only very locally
        //      affect the accuracy of interpolation, which I think is not
        //      important at this stage. VV, 20/Oct/2016.

        // Create counter for all donors on this processor
        label nDonors = 0;

        // Loop through all regions
        forAll(regions_, regionI)
        {
            // Get the list of donors
            const donorAcceptorList& curDonors = regions_[regionI].donors();

            // Loop through all donors
            forAll (curDonors, dI)
            {





                if (curDonors[dI].acceptorProcNo() == Pstream::myProcNo())
                {
                    // Local donor and acceptor
                    // Data has already been recorded from the acceptor side
                    // Add debug check?
                }
                else
                {
                    // Get donor cell index
                    const label& donorCellI = curDonors[dI].donorCell();

                    // Record local donor for a remote acceptor:
                    // Data will be sent to the master
                    remDonors[nRemoteDonors] = donorCellI;

                    procRemoteDonors[nRemoteDonors] = curDonors[dI];

                    // Get neighbouring cells of this cells
                    const labelList& curNbrs = cc[donorCellI];

                    // Mark neighbouring donors based on extension level
                    markNeighbouringDonors
                    (
                        curNbrs,
                        cc,
                        donorMask,
                        extDonors,
                        interpolationPtr_->extensionLevel()
                    );

                    // Transfer the list from hash set into
                    // localNeighbouringDonors list for this master donor
                    // (acceptor)
                    remNeiDonors[nRemoteDonors] = extDonors.toc();

                    // Clear the hash set, keeping the allocated memory
                    extDonors.clear();

                    // Increment the counter
                    nRemoteDonors++;
                }
            }
        }






        // Parallel run: handle remote donors
        // Algorithm:
        // 1) Collect local donors and mark them directly in the list
        // 2) For remote donors providing data for the local processor,
        //    mark the processor they come from and their location
        //    in the local interpolation list
        // 3) Communicate the local donor and acceptor data to the master,
        //    which assembles the addressing

        // Acceptor list will be filled in two loops:
        // - first, local donors are inserted into the acceptor array
        // - second, remote donors are filled and sent to the master processor
        // - third, master processor reshuffles the processor donor data
        //   based on the acceptor processor index
        //   and sends the data to acceptor processors
        // - fourth, the remote donor data is unpacked into the acceptor array
        // - fifth, donor-to-acceptor interpolation is executed in a straight
        //   loop

        // Create remote donor addressing
        remoteDonorsPtr_ = new labelList(dc.size(), -1);
        labelList& remDonors = *remoteDonorsPtr_;

        // Create remote neighbouring donors addressing
        remoteNeighbouringDonorsPtr_ = new labelListList(dc.size());
        labelListList& remNeiDonors = *remoteNeighbouringDonorsPtr_;

        // Prepare remote donor list for master processor to
        // calculate addressing.  This list contains local
        // donor cells whose acceptor is on a different processor
        donorAcceptorListList globalRemoteDonors(Pstream::nProcs());
        donorAcceptorList& procRemoteDonors =
            globalRemoteDonors[Pstream::myProcNo()];

        // Size the remote donor list for the local processor
        procRemoteDonors.setSize(dc.size());

        // Create remote acceptor addressing
        remoteAcceptorAddrPtr_ = new labelList(ac.size(), -1);
        labelList& remAcceptorAddr = *remoteAcceptorAddrPtr_;

        // Prepare local acceptor list for master processor
        // to calculate addressing.  This list contains local
        // acceptor cells whose donor is on a different processor
        donorAcceptorListList globalRemoteAcceptors(Pstream::nProcs());
        donorAcceptorList& procRemoteAcceptors =
            globalRemoteAcceptors[Pstream::myProcNo()];

        // Size the remote acceptor list for local processor
        procRemoteAcceptors.setSize(ac.size());

        // Count local donors to local acceptors
        label nLocalAddr = 0;

        // Count remote donors
        label nRemoteDonors = 0;

        // Count remote acceptors
        label nRemoteAcceptors = 0;

        // Collect and mark acceptors that are local

        // Note:
        // In order to avoid creating a donorAcceptorList across all regions
        // use a counter which follows the number of analysed acceptors
        // for the local mesh across all regions
        // HJ, 1/May/2015
        label acceptorI = 0;

        forAll (regions_, regionI)
        {
            // Analyse the acceptor list for local acceptors
            const donorAcceptorList& curAcceptors =
                regions_[regionI].acceptors();

            // Create working hash set of extended neighbouring donors
            labelHashSet extDonors
            (
                polyMesh::facesPerCell_*
                Foam::pow(interpolationPtr_->extensionLevel() + 1, 3)
            );

            forAll (curAcceptors, aI)
            {
                if (curAcceptors[aI].donorProcNo() == Pstream::myProcNo())
                {
                    // Get donor cell index
                    const label& donorCellI = curAcceptors[aI].donorCell();

                    // Local donor and acceptor
                    locDonors[nLocalAddr] = donorCellI;
                    locDonorAddr[nLocalAddr] = acceptorI;

                    // Get neighbouring cells of this cells
                    const labelList& curNbrs = cc[donorCellI];

                    // Mark neighbouring donors based on extension level
                    markNeighbouringDonors
                    (
                        curNbrs,
                        cc,
                        donorMask,
                        extDonors,
                        interpolationPtr_->extensionLevel()
                    );

                    // Transfer the list from hash set into
                    // localNeighbouringDonors list for this master donor
                    // (acceptor)
                    locNeiDonors[nLocalAddr] = extDonors.toc();

                    // Clear the hash set, keeping the allocated memory
                    extDonors.clear();

                    // Increment the counter
                    nLocalAddr++;
                }
                else
                {
                    // Record local acceptor with a remote donor:
                    // Data will be received from master
                    remAcceptorAddr[nRemoteAcceptors] = acceptorI;

                    procRemoteAcceptors[nRemoteAcceptors] = curAcceptors[aI];

                    nRemoteAcceptors++;
                }

                // Increment acceptor counter
                acceptorI++;
            }

            // Analyse the donor list
            const donorAcceptorList& curDonors = regions_[regionI].donors();

            forAll (curDonors, dI)
            {
                if (curDonors[dI].acceptorProcNo() == Pstream::myProcNo())
                {
                    // Local donor and acceptor
                    // Data has already been recorded from the acceptor side
                    // Add debug check?
                }
                else
                {
                    // Get donor cell index
                    const label& donorCellI = curDonors[dI].donorCell();

                    // Record local donor for a remote acceptor:
                    // Data will be sent to the master
                    remDonors[nRemoteDonors] = donorCellI;

                    procRemoteDonors[nRemoteDonors] = curDonors[dI];

                    // Get neighbouring cells of this cells
                    const labelList& curNbrs = cc[donorCellI];

                    // Mark neighbouring donors based on extension level
                    markNeighbouringDonors
                    (
                        curNbrs,
                        cc,
                        donorMask,
                        extDonors,
                        interpolationPtr_->extensionLevel()
                    );

                    // Transfer the list from hash set into
                    // localNeighbouringDonors list for this master donor
                    // (acceptor)
                    remNeiDonors[nRemoteDonors] = extDonors.toc();

                    // Clear the hash set, keeping the allocated memory
                    extDonors.clear();

                    // Increment the counter
                    nRemoteDonors++;
                }
            }
        } // End of for all regions

        // Reset the size of lists
        locDonors.setSize(nLocalAddr);
        locNeiDonors.setSize(nLocalAddr);
        locDonorAddr.setSize(nLocalAddr);

        remAcceptorAddr.setSize(nRemoteAcceptors);
        procRemoteAcceptors.setSize(nRemoteAcceptors);

        // Reset the size of lists
        remDonors.setSize(nRemoteDonors);
        remNeiDonors.setSize(nRemoteDonors);
        procRemoteDonors.setSize(nRemoteDonors);

        // Gather remote donor and acceptor data before identification
        Pstream::gatherList(globalRemoteDonors);
        Pstream::gatherList(globalRemoteAcceptors);

        // Now all data is available on master.  Perform analysis
        if (Pstream::master())
        {
            // Make a map of donor cells on all processors
            // May key is the donor cell index and its value is the
            // index in donor list
            List<Map<label> > globalDonorLookup(Pstream::nProcs());

            forAll (globalDonorLookup, procI)
            {
                Map<label>& procDonorLookup = globalDonorLookup[procI];

                // Get processor donor data
                const donorAcceptorList& procDonors =
                    globalRemoteDonors[procI];

                forAll (procDonors, donorI)
                {
                    procDonorLookup.insert
                    (
                        procDonors[donorI].donorCell(),
                        donorI
                    );
                }
            }

            // globalDonorLookup will be used for fast lookup without search
            // HJ, 25/May/2015

            // Allocate the storage for addressing

            // Create global accept from processor addressing
            globalAcceptFromProcPtr_ = new labelListList(Pstream::nProcs());
            labelListList& globalAcceptProc = *globalAcceptFromProcPtr_;

            // Create global accept from cell addressing
            globalAcceptFromCellPtr_ = new labelListList(Pstream::nProcs());
            labelListList& globalAcceptCell = *globalAcceptFromCellPtr_;

            // For all processors, find which donor from which other
            // processor contains the data
            forAll (globalRemoteAcceptors, procI)
            {
                const donorAcceptorList& procAcceptors =
                    globalRemoteAcceptors[procI];

                // Resize addressing
                labelList& curAcceptFromProc = globalAcceptProc[procI];
                labelList& curAcceptFromCell = globalAcceptCell[procI];

                curAcceptFromProc.setSize(procAcceptors.size());
                curAcceptFromCell.setSize(procAcceptors.size());

                // Collect data from other processors
                forAll (procAcceptors, accI)
                {
                    const donorAcceptor& da = procAcceptors[accI];

                    // Get donor processor list of donors
                    const label donorProc = da.donorProcNo();

                    // Search to find which donor to pick up.

                    const Map<label>& procDonorLookup =
                        globalDonorLookup[donorProc];

                    // Find my donor cell index in the map
                    const Map<label>::const_iterator iter =
                        procDonorLookup.find(da.donorCell());

                    if (iter != procDonorLookup.end())
                    {
                        // Found donor.  Record addressing
                        curAcceptFromProc[accI] = da.donorProcNo();
                        curAcceptFromCell[accI] = iter();
                    }
                    else
                    {
                        FatalErrorIn
                        (
                            "void oversetMesh::calcInterpolationAddressing() const"
                        )   << "Cannot find pairing for acceptor " << da
                            << abort(FatalError);
                    }
                }
            }

            // Check parallel addressing
            forAll (globalAcceptProc, procI)
            {
                if
                (
                    min(globalAcceptProc[procI]) < 0
                 || min(globalAcceptCell[procI]) < 0
                )
                {
                    FatalErrorIn
                    (
                        "void oversetMesh::calcInterpolationAddressing() const"
                    )   << "Error in parallel addressing assembly "
                        << " for processor " << procI
                        << abort(FatalError);
                }
            }
        }
    }
    else
    {
        // Serial run: all donors and acceptors are local

        // Count local donors to local acceptors
        label nLocalAddr = 0;

        // Go through all acceptors and pick up the donor
        forAll (regions_, regionI)
        {
            // Acceptors
            const donorAcceptorList& curAcceptors =
                regions_[regionI].acceptors();

            // Create working hash set of extended neighbouring donors
            labelHashSet extDonors
            (
                polyMesh::facesPerCell_*
                Foam::pow(interpolationPtr_->extensionLevel() + 1, 3)
            );

            forAll (curAcceptors, aI)
            {
                // Get donor cell index
                const label& donorCellI = curAcceptors[aI].donorCell();

                // Grab local donor from the list
                locDonors[nLocalAddr] = donorCellI;
                locDonorAddr[nLocalAddr] = nLocalAddr;

                // Get neighbouring cells of this cells
                const labelList& curNbrs = cc[donorCellI];

                // Mark neighbouring donors based on extension level
                markNeighbouringDonors
                (
                    curNbrs,
                    cc,
                    donorMask,
                    extDonors,
                    interpolationPtr_->extensionLevel()
                );

                // Transfer the list from hash set into
                // localNeighbouringDonors list for this master donor
                // (acceptor)
                locNeiDonors[nLocalAddr] = extDonors.toc();

                // Clear the hash set, keeping the allocated memory
                extDonors.clear();

                // Increment the counter
                nLocalAddr++;
            }
        }
    }

    // Check serial addressing
    if (min(locDonorAddr) < 0)
    {
        FatalErrorIn("void oversetMesh::calcInterpolationAddressing() const")
            << "Error in local addressing assembly"
            << abort(FatalError);
    }

    // Force calculation of domain markup fields for post-processing
    // HJ, 9/Apr/2013
    oversetTypes();
}


void Foam::oversetMesh::clearOut() const
{
    deleteDemandDrivenData(acceptorCellsPtr_);
    deleteDemandDrivenData(donorCellsPtr_);
    deleteDemandDrivenData(holeCellsPtr_);

    deleteDemandDrivenData(oversetTypesPtr_);
    deleteDemandDrivenData(regionIDPtr_);

    deleteDemandDrivenData(gammaPtr_);
    deleteDemandDrivenData(gammaExtPtr_);
    deleteDemandDrivenData(sGammaPtr_);

    deleteDemandDrivenData(fringeFacesPtr_);
    deleteDemandDrivenData(fringeFaceCellsPtr_);
    deleteDemandDrivenData(fringeFaceFlipsPtr_);

    deleteDemandDrivenData(holeFacesPtr_);
    deleteDemandDrivenData(holeFaceCellsPtr_);
    deleteDemandDrivenData(holeFaceFlipsPtr_);
    deleteDemandDrivenData(holeInternalFacesPtr_);
    deleteDemandDrivenData(acceptorInternalFacesPtr_);

    deleteDemandDrivenData(mapPtr_);
    deleteDemandDrivenData(nDonorsToProcessorMapPtr_);
    deleteDemandDrivenData(donorsToAcceptorsMapPtr_);
}


void Foam::oversetMesh::markNeighbouringDonors
(
    const labelList& nbrCells,
    const labelListList& cellCells,
    const boolList& donorMask,
    labelHashSet& extendedNeighbouringDonors,
    label curLevel
) const
{
    if (curLevel > -1)
    {
        // For each neighbouring cell, loop through its neighbours
        forAll (nbrCells, i)
        {
            // Get neighbouring cell index
            const label& nbrCellI = nbrCells[i];

            // Check if the cell is eligible donor (not an acceptor or hole)
            if (donorMask[nbrCellI])
            {
                // Insert this neighbouring cell in the extended neighbouring
                // donor hash set
                extendedNeighbouringDonors.insert(nbrCellI);
            }

            // Recursive call with decremented level and neighbours of this
            // neighbouring cells
            markNeighbouringDonors
            (
                cellCells[nbrCellI],
                cellCells,
                donorMask,
                extendedNeighbouringDonors,
                curLevel - 1
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::oversetMesh::acceptorCells() const
{
    if (!acceptorCellsPtr_)
    {
        calcCellClassification();
    }

    return *acceptorCellsPtr_;
}


const Foam::labelList& Foam::oversetMesh::donorCells() const
{
    if (!donorCellsPtr_)
    {
        calcCellClassification();
    }

    return *donorCellsPtr_;
}


const Foam::labelList& Foam::oversetMesh::holeCells() const
{
    if (!holeCellsPtr_)
    {
        calcCellClassification();
    }

    return *holeCellsPtr_;
}


const Foam::volScalarField& Foam::oversetMesh::oversetTypes() const
{
    if (!oversetTypesPtr_)
    {
        calcDomainMarkup();
    }

    return *oversetTypesPtr_;
}


const Foam::labelList& Foam::oversetMesh::regionID() const
{
    if (!regionIDPtr_)
    {
        calcDomainMarkup();
    }

    return *regionIDPtr_;
}


const Foam::volScalarField& Foam::oversetMesh::gamma() const
{
    if (!gammaPtr_)
    {
        calcGamma();
    }

    return *gammaPtr_;
}


const Foam::volScalarField& Foam::oversetMesh::gammaExt() const
{
    if (!gammaExtPtr_)
    {
        calcGamma();
    }

    return *gammaExtPtr_;
}


const Foam::surfaceScalarField& Foam::oversetMesh::sGamma() const
{
    if (!sGammaPtr_)
    {
        calcGamma();
    }

    return *sGammaPtr_;
}


const Foam::labelList& Foam::oversetMesh::fringeFaces() const
{
    if (!fringeFacesPtr_)
    {
        calcFringeFaces();
    }

    return *fringeFacesPtr_;
}


const Foam::labelList& Foam::oversetMesh::fringeFaceCells() const
{
    if (!fringeFaceCellsPtr_)
    {
        calcFringeFaces();
    }

    return *fringeFaceCellsPtr_;
}


const Foam::boolList& Foam::oversetMesh::fringeFaceFlips() const
{
    if (!fringeFaceFlipsPtr_)
    {
        calcFringeFaces();
    }

    return *fringeFaceFlipsPtr_;
}


const Foam::labelList& Foam::oversetMesh::holeFaces() const
{
    if (!holeFacesPtr_)
    {
        calcHoleFaces();
    }

    return *holeFacesPtr_;
}


const Foam::labelList& Foam::oversetMesh::holeFaceCells() const
{
    if (!holeFaceCellsPtr_)
    {
        calcHoleFaces();
    }

    return *holeFaceCellsPtr_;
}


const Foam::boolList& Foam::oversetMesh::holeFaceFlips() const
{
    if (!holeFaceFlipsPtr_)
    {
        calcHoleFaces();
    }

    return *holeFaceFlipsPtr_;
}


const Foam::labelList& Foam::oversetMesh::holeInternalFaces() const
{
    if (!holeInternalFacesPtr_)
    {
        calcHoleFaces();
    }

    return *holeInternalFacesPtr_;
}


const Foam::labelList& Foam::oversetMesh::acceptorInternalFaces() const
{
    if (!acceptorInternalFacesPtr_)
    {
        calcFringeFaces();
    }

    return *acceptorInternalFacesPtr_;
}


const Foam::mapDistributeOversetDonors& Foam::oversetMesh::map() const
{
    if (!mapPtr_)
    {
        calcInterpolationAddressing();
    }

    return *mapPtr_;
}


const Foam::labelListList& Foam::oversetMesh::nDonorsToProcessorMap() const
{
    if (!nDonorsToProcessorMapPtr_)
    {
        calcInterpolationAddressing();
    }

    return *nDonorsToProcessorMapPtr_;
}


const Foam::labelListList& Foam::oversetMesh::donorsToAcceptorsMap() const
{
    if (!donorsToAcceptorsMapPtr_)
    {
        calcInterpolationAddressing();
    }

    return *donorsToAcceptorsMapPtr_;
}


// ************************************************************************* //
