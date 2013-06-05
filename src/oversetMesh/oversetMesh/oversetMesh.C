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

    splitPtr_(NULL),
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

    fringeAddressingPtr_(NULL),
    interpolationPtr_(NULL)
{
    Info << "Creating oversetMesh" << endl;

    // Read regions
    List<dictionary> regionDicts(dict_.lookup("regions"));

    regions_.setSize(regionDicts.size());

    forAll (regionDicts, dictI)
    {
        regions_.set
        (
            dictI,
            new oversetRegion
            (
                mesh,
                *this,
                regionDicts[dictI]
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::oversetMesh::~oversetMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::regionSplit& Foam::oversetMesh::split() const
{
    if (!splitPtr_)
    {
        splitPtr_ = new regionSplit(mesh());

        // Check setup of regions
        if (splitPtr_->nRegions() != regions_.size())
        {
            FatalErrorIn
            (
                "const regionSplit& oversetMesh::split() const"
            )   << "Number of regions defined by oversetMesh does not "
                << "match the region split." << nl
                << "nRegions = " << regions_.size()
                << " split = " << splitPtr_->nRegions()
                << abort(FatalError);
        }
    }

    return *splitPtr_;
}


const Foam::oversetRegion& Foam::oversetMesh::region
(
    const label index
) const
{
    if (index < 0 || index >= split().nRegions())
    {
        FatalErrorIn("const oversetRegion& Foam::oversetMesh::region")
            << "Invalid index: should be between 0 and "
            << split().nRegions()
            << abort(FatalError);
    }

    forAll (regions_, rI)
    {
        if (regions_[rI].index() == index)
        {
            return regions_[rI];
        }
    }

    FatalErrorIn("const oversetRegion& Foam::oversetMesh::region")
        << "Cannot find region index "
        << split().nRegions()
        << abort(FatalError);

    // Dummy return to keep compiler happy
    return regions_[0];
}


const Foam::oversetInterpolation& Foam::oversetMesh::interpolation() const
{
    if (!interpolationPtr_)
    {
        interpolationPtr_ = oversetInterpolation::New
        (
            *this,
            dict_.subDict("interpolation")
        ).ptr();
    }

    return *interpolationPtr_;
}


bool Foam::oversetMesh::movePoints() const
{
    // Perform appropriate updates on search and fringe
    // HJ, 3/Apr/2013

    return false;
}


bool Foam::oversetMesh::updateMesh(const mapPolyMesh&) const
{
    // Perform appropriate updates on search and fringe
    // HJ, 3/Apr/2013

    return true;
}


// ************************************************************************* //
