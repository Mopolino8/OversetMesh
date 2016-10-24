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
#include "findRefCell.H"

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
    singleUpdate_
    (
        dict_.lookupOrDefault<Switch>("singleUpdatePerTimeStep", false)
    ),
    curTimeIndex_(-1),

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

    mapPtr_(NULL),
    remoteDonorToLocalAcceptorAddrPtr_(NULL),

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


void Foam::oversetMesh::setRefCell
(
    const volScalarField& field,
    const dictionary& dict,
    label& refCelli,
    scalar& refValue
) const
{
    // Set the reference cell using existing machinery
    Foam::setRefCell(field, dict, refCelli, refValue);

    // If the fringe is not coupled for this field (explicit overset treatment),
    // then we must avoid setting the reference cell because explicit overset
    // acts as a cell centred fixed value boundary condition. We achieve this by
    // invalidating refCelli, thus avoiding setting the reference for this field
    // when fvMatrix::setReference is called.
    // Note that we are relying on the fact that oversetPatch must come first in
    // the boundaryField as this is already necessary for correct flux
    // reconstruction.
    if (!field.boundaryField()[0].coupled())
    {
        refCelli = -1;
    }
}


bool Foam::oversetMesh::movePoints() const
{
    // Perform appropriate updates on search and fringe
    // HJ, 3/Apr/2013
    Info<< "Overset mesh motion update" << endl;

    // Get time index
    const label globalTimeIndex = mesh().time().timeIndex();

    // Update regions, interpolation and clear out the address if the time local
    // time index is smaller than global one
    if (curTimeIndex_ < globalTimeIndex)
    {
        forAll (regions_, regionI)
        {
            regions_[regionI].update();
        }

        interpolationPtr_->update();

        clearOut();

        // Update the local time step only if singleUpdate switch is turned on.
        // This update controls whether the overset addressing will be updated
        // once per time step or multiple times per time step.
        if (singleUpdate_)
        {
            curTimeIndex_ = globalTimeIndex;
        }
    }

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
