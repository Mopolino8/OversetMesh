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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::oversetRegion::calcRegionCells() const
{
    if (regionCellsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcRegionCells() const")
            << "Region cells already calculated"
            << abort(FatalError);
    }

    // Create the octree search from all donor regions
    const regionSplit& rs = oversetMesh_.split();

    // Check region index
    if (index_ < 0 || index_ >= rs.nRegions())
    {
        FatalErrorIn("void oversetRegion::calcRegionCells() const")
            << "Invalid overset region index " << index_
            << ".  Available number of regions = " << rs.nRegions()
            << abort(FatalError);
    }

    // Count number of cells in region
    label nCells = 0;

    forAll (rs, cellI)
    {
        if (rs[cellI] == index_)
        {
            nCells++;
        }
    }

    if (nCells == 0)
    {
        FatalErrorIn("void oversetRegion::calcRegionCells() const")
            << "No cells found in region index " << index_
            << abort(FatalError);
    }

    regionCellsPtr_ = new labelList(nCells);
    labelList& rc = *regionCellsPtr_;

    // Reset counter and collect cells
    nCells = 0;
    
    forAll (rs, cellI)
    {
        if (rs[cellI] == index_)
        {
            rc[nCells] = cellI;
            nCells++;
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

    const label nRegions = oversetMesh_.split().nRegions();

    acceptorRegionsPtr_ = new labelList(nRegions);
    labelList& accRegions = *acceptorRegionsPtr_;

    label nAccRegions = 0;

    // Go through all regions apart from the current and check if
    // this region appears in the donor list
    for (label regionI = 0; regionI < nRegions; regionI++)
    {
        // Skip current region
        if (regionI == index())
        {
            continue;
        }

        const labelList& remoteDonors =
            oversetMesh_.region(regionI).donorRegions();

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


void Foam::oversetRegion::calcAcceptorCells() const
{
    if (acceptorCellsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcAcceptorCells() const")
            << "Acceptor cells already calculated"
            << abort(FatalError);
    }

    // Grab acceptor cells from fringe algorithm
    acceptorCellsPtr_ = new labelList(fringePtr_().acceptors());
}


void Foam::oversetRegion::calcDonorCells() const
{
    if (donorCellsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcDonorCells() const")
            << "Donor cells already calculated"
            << abort(FatalError);
    }

    // Algorithm:
    // - go through the hierarchy of donor regions and search for cell
    //   containing the acceptor

    const labelList& dr = donorRegions();

    // Get list of acceptor cell labels
    const labelList& a = acceptors();

    // Get cell centres
    const vectorField& cc = mesh_.cellCentres();

    donorCellsPtr_ = new labelList(a.size(), -1);
    labelList& d = *donorCellsPtr_;

    // Go through all donor regions
    forAll (dr, drI)
    {
        if (dr[drI] == index())
        {
            FatalErrorIn("void oversetRegion::calcDonorCells() const")
                << "Region " << index() << " specified as the donor "
                << "of itself.  List of donors: " << dr << nl
                << "This is not allowed: please check oversetMesh definition"
                << abort(FatalError);
        }

        // Get donor region tree
//         const indexedOctree<treeDataCell>& tree =
//             oversetMesh_.region(donorRegions_[drI]).cellTree();

//         scalar span = tree.bb().mag();
//         scalar greatSpan = sqr(GREAT);

        const labelList& curDonors =
            oversetMesh_.region(donorRegions_[drI]).eligibleDonors();

        // Search donor region with all unassigned cells
        forAll (a, aI)
        {
            if (d[aI] == -1)
            {
                const vector& curCentre = cc[a[aI]];

                // n-squared search - testing

                scalar minDistance = GREAT;
                label dc = -1;
                scalar magSqrDist;

                forAll (curDonors, donorI)
                {
                    magSqrDist = magSqr(cc[curDonors[donorI]] - curCentre);

                    if (magSqrDist < minDistance)
                    {
                        dc = curDonors[donorI];
                        minDistance = magSqrDist;
                    }
                }

                // Check for point in cell
                if (dc > -1)
                {
//                     if (mesh_.pointInCell(curCentre, dc))
                    if (mesh_.pointInCellBB(curCentre, dc))
                    {
                        d[aI] = dc;
                    }
                }
            }
        }
    }
    forAll (d, i)
    {
        if (d[i] == -1)
        {
            Info<< "Missing donor for cell " << a[i] << endl;
        }
    }

    // Check for unpaired donors
    if (!d.empty())
    {
        if (min(d) < 0)
        {
            FatalErrorIn("void oversetRegion::calcDonorCells() const")
                << "Region " << index() << " has unmatched acceptors"
                << abort(FatalError);
        }
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

    // Go through all regions apart from the current and mark all hole cells
    // using their hole boundary patch inside search

    const label nRegions = oversetMesh_.split().nRegions();

    for (label regionI = 0; regionI < nRegions; regionI++)
    {
        // Skip current region
        if (regionI == index())
        {
            continue;
        }

        const oversetRegion& otherRegion = oversetMesh_.region(regionI);

        if (!otherRegion.holePatches())
        {
            continue;
        }

     // Get reference to hole search
        const triSurfaceSearch& holeSearch = otherRegion.holeTriSurfSearch();

        boolList regionInside = holeSearch.calcInside(localC);

        // Note: hole mask has the size of all mesh cells and regionInside
        // only of the size of local region
        forAll (regionInside, i)
        {
            holeMask[rc[i]] |= regionInside[i];
        } 
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

    // Remove all acceptor cells from region
    const labelList& ac = acceptors();

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


void Foam::oversetRegion::calcCellTree() const
{
    if (cellTreePtr_)
    {
        FatalErrorIn("void oversetRegion::calcCellTree() const")
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

    cellTreePtr_ = new indexedOctree<treeDataCell>
    (
        treeDataCell
        (
            false,  //  Cache bb.  Reconsider for moving mesh cases
            mesh_,
            regionCells()
        ),
        overallBb,  // overall search domain
        8,          // maxLevel
        10,         // leafsize
        3.0         // duplicity
    );
}


void Foam::oversetRegion::clearOut()
{
    deleteDemandDrivenData(regionCellsPtr_);
    deleteDemandDrivenData(acceptorRegionsPtr_);

    deleteDemandDrivenData(acceptorCellsPtr_);
    deleteDemandDrivenData(donorCellsPtr_);
    deleteDemandDrivenData(holeCellsPtr_);
    deleteDemandDrivenData(eligibleDonorCellsPtr_);

    deleteDemandDrivenData(holeTriMeshPtr_);
    deleteDemandDrivenData(triSurfSearchPtr_);
    deleteDemandDrivenData(cellTreePtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetRegion::oversetRegion
(
    const fvMesh& mesh,
    const oversetMesh& oversetMesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    oversetMesh_(oversetMesh),
    index_(readLabel(dict.lookup("index"))),
    donorRegions_(dict.lookup("donorRegions")),
    holePatchNames_(dict.lookup("holePatches")),
    fringePtr_(),
    regionCellsPtr_(NULL),
    acceptorRegionsPtr_(NULL),

    acceptorCellsPtr_(NULL),
    donorCellsPtr_(NULL),
    holeCellsPtr_(NULL),
    eligibleDonorCellsPtr_(NULL),

    holeTriMeshPtr_(NULL),
    triSurfSearchPtr_(NULL),
    cellTreePtr_(NULL)
{
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
    if (!regionCellsPtr_)
    {
        calcRegionCells();
    }

    return *regionCellsPtr_;
}


const Foam::labelList& Foam::oversetRegion::acceptorRegions() const
{
    if (!acceptorRegionsPtr_)
    {
        calcAcceptorRegions();
    }

    return *acceptorRegionsPtr_;
}


const  Foam::labelList& Foam::oversetRegion::acceptors() const
{
    if (!acceptorCellsPtr_)
    {
        calcAcceptorCells();
    }

    return *acceptorCellsPtr_;
}


const  Foam::labelList& Foam::oversetRegion::donors() const
{
    if (!donorCellsPtr_)
    {
        calcDonorCells();
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


const Foam::triSurfaceMesh& Foam::oversetRegion::holeTriMesh() const
{
    if (!holeTriMeshPtr_)
    {
        // Collect all hole patches and make a triangular surface
        labelHashSet holePatches;

        forAll (holePatchNames_, nameI)
        {
            polyPatchID curHolePatch
            (
                holePatchNames_[nameI],
                mesh().boundaryMesh()
            );

            if (curHolePatch.active())
            {
                holePatches.insert(curHolePatch.index());
            }
            else
            {
                FatalErrorIn
                (
                    "const triSurfaceMesh& oversetRegion::holeTriMesh() const"
                )   << "Patch "  << holePatchNames_[nameI]
                    << " cannot be found.  Available patch names: "
                    << mesh().boundaryMesh().names()
                    << abort(FatalError);
            }
        }

        if (holePatches.empty())
        {
            FatalErrorIn
            (
                "const triSurfaceMesh& oversetRegion::holeTriMesh() const"
            )   << "No hole patches detected: cannot perform hole cutting"
                << abort(FatalError);
        }

        triSurface ts = triSurfaceTools::triangulate
        (
            mesh().boundaryMesh(),
            holePatches
        );

        // Invert faces
        triFaceList invertedFaces(ts.size());

        forAll (ts, tsI)
        {
            invertedFaces[tsI] = ts[tsI].reverseFace();
        }

        triSurface invertedTs
        (
            invertedFaces,
            ts.points()
        );

//         invertedTs.write(mesh().time().caseName() + ".stl");

        holeTriMeshPtr_ = new triSurfaceMesh
        (
            IOobject
            (
                "oversetMeshTri.ftr",
                mesh().time().constant(),  // instance
                "triSurface",              // local
                mesh(),                    // registry
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            invertedTs
        );

        
    }

    return *holeTriMeshPtr_;
}


const Foam::triSurfaceSearch& Foam::oversetRegion::holeTriSurfSearch() const
{
    if (!triSurfSearchPtr_)
    {
        triSurfSearchPtr_ = new triSurfaceSearch
        (
            holeTriMesh()
        );
    }

    return *triSurfSearchPtr_;
}


const Foam::indexedOctree<Foam::treeDataCell>&
Foam::oversetRegion::cellTree() const
{
    if (!cellTreePtr_)
    {
        calcCellTree();
    }

    return *cellTreePtr_;
}


// ************************************************************************* //
