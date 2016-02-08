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
    ANY WARRANTY; without even the implied wrranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "oversetMesh.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "demandDrivenData.H"
#include "oversetPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetMesh, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::oversetMesh::oversetMesh(const fvMesh& mesh)
:
    MeshObject<fvMesh, oversetMesh>(mesh),
    dict_
    (
        IOobject
        (
            "oversetMeshDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    regions_(),
    holePatchNames_(dict_.lookup("holePatches")),

    acceptorCellsPtr_(NULL),
    donorCellsPtr_(NULL),
    holeCellsPtr_(NULL),

    oversetTypesPtr_(NULL),
    regionIDPtr_(NULL),

    gammaPtr_(NULL),
    gammaExtPtr_(NULL),
    sGammaPtr_(NULL),

    fringeFacesPtr_(NULL),
    fringeFaceCellsPtr_(NULL),
    fringeFaceFlipsPtr_(NULL),
    holeFacesPtr_(NULL),
    holeFaceCellsPtr_(NULL),
    holeFaceFlipsPtr_(NULL),
    holeInternalFacesPtr_(NULL),
    acceptorInternalFacesPtr_(NULL),

    localDonorsPtr_(NULL),
    localNeighbouringDonorsPtr_(NULL),
    localDonorAddrPtr_(NULL),
    remoteDonorsPtr_(NULL),
    remoteNeighbouringDonorsPtr_(NULL),
    remoteAcceptorAddrPtr_(NULL),
    globalAcceptFromProcPtr_(NULL),
    globalAcceptFromCellPtr_(NULL),
    interpolationPtr_
    (
        oversetInterpolation::New
        (
            *this,
            dict_.subDict("interpolation")
        )
    )
{
    Info << "Creating oversetMesh" << endl;

    // Read regions
    PtrList<entry> regionEntries(dict_.lookup("regions"));

    regions_.setSize(regionEntries.size());
    forAll (regionEntries, regionI)
    {
        regions_.set
        (
            regionI,
            new oversetRegion
            (
                regionEntries[regionI].keyword(),
                regionI,
                mesh,
                *this,
                regionEntries[regionI].dict()
            )
        );
    }

    // Overset patch must come first for consistent handling of patch flux on
    // coupled boundaries (see oversetFvPatchField::patchFlux() member function)
    if (!isA<oversetPolyPatch>(mesh.boundaryMesh()[0]))
    {
        FatalErrorIn("oversetMesh::oversetMesh(const fvMesh& mesh)")
          << "Overset patch needs to come first for consistent reconstruction"
          << nl << " of fringe face fluxes on coupled boundaries."
          << nl << "First patch is: " << mesh.boundaryMesh()[0].name()
          << nl << "with type: " << mesh.boundaryMesh()[0].type()
          << abort(FatalError);
    }

    // Check for duplicate region names
    // TODO: HJ
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::oversetMesh::~oversetMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::oversetInterpolation& Foam::oversetMesh::interpolation() const
{
    return interpolationPtr_();
}


bool Foam::oversetMesh::movePoints() const
{
    // Perform appropriate updates on search and fringe
    // HJ, 3/Apr/2013
    Info<< "Overset mesh motion update" << endl;

    forAll (regions_, regionI)
    {
        regions_[regionI].update();
    }

    // Update interpolation (recalculate weights)
    interpolationPtr_->update();

    clearOut();

    return false;
}


bool Foam::oversetMesh::updateMesh(const mapPolyMesh&) const
{
    // Perform appropriate updates on search and fringe
    // HJ, 3/Apr/2013
    Info<< "Overset topo update" << endl;

    forAll (regions_, regionI)
    {
        regions_[regionI].update();
    }

    // Update interpolation (recalculate weights)
    interpolationPtr_->update();

    clearOut();

    return true;
}


// ************************************************************************* //
