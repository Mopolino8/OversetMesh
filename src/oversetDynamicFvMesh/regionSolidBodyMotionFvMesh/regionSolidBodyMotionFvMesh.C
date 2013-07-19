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

#include "regionSolidBodyMotionFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "transformField.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionSolidBodyMotionFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        regionSolidBodyMotionFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::regionSolidBodyMotionFvMesh::calcMotionMasks() const
{
    if (regionMotionMasksPtr_)
    {
        FatalErrorIn
        (
            "void regionSolidBodyMotionFvMesh::calcMotionMasks() const"
        )   << "Motion masks already calculated"
            << abort(FatalError);
    }

    regionMotionMasksPtr_ = new FieldField<Field, scalar>(rs_.nRegions());
    FieldField<Field, scalar>& rmm = *regionMotionMasksPtr_;

    const labelListList& cp = cellPoints();

    // Size all masks and set to zero
    for (label regionI = 0; regionI < rs_.nRegions(); regionI++)
    {
        rmm.set(regionI, new scalarField(nPoints(), scalar(0)));

//         scalarField& curMask = rmm[regionI];
//         curMask.setSize(nPoints(), scalar(0));
    }

    // Loop through all cells and mark appropriate mask
    const labelList& regionCells = rs_;

    forAll (regionCells, cellI)
    {
        // Grab mask belonging to the cell
        scalarField& curMask = rmm[regionCells[cellI]];

        // Mark all points of the cell
        const labelList& curPoints = cp[cellI];

        forAll (curPoints, pointI)
        {
            curMask[curPoints[pointI]] = 1;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolidBodyMotionFvMesh::regionSolidBodyMotionFvMesh
(
    const IOobject& io
)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    rs_(*this),
    regionSBMFs_(rs_.nRegions()),
    regionMotionMasksPtr_(NULL),
    undisplacedPoints_
    (
        IOobject
        (
            "points",
            time().constant(),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{
    // Read motion function for all regions
    List<dictionary> motionDicts(dynamicMeshCoeffs_.lookup("motionFunctions"));

    if (motionDicts.size() != rs_.nRegions())
    {
        FatalIOErrorIn
        (
            "regionSolidBodyMotionFvMesh::regionSolidBodyMotionFvMesh\n"
            "(\n"
            "    const IOobject& io\n"
            ")",
            dynamicMeshCoeffs_
        )   << "Number of motion SBMF dictionaries does not correspond to "
            << "the number of disconnected regions in the mesh." << nl
            << "nRegions: " << rs_.nRegions()
            << " dicts: " << motionDicts.size()
            << abort(FatalIOError);
    }

    forAll (regionSBMFs_, regionI)
    {
        regionSBMFs_.set
        (
            regionI,
            solidBodyMotionFunction::New(motionDicts[regionI], time())
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolidBodyMotionFvMesh::~regionSolidBodyMotionFvMesh()
{
    deleteDemandDrivenData(regionMotionMasksPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::scalarField&
Foam::regionSolidBodyMotionFvMesh::motionMask(const label regionI) const
{
    if (!regionMotionMasksPtr_)
    {
        calcMotionMasks();
    }

    return (*regionMotionMasksPtr_)[regionI];
}


bool Foam::regionSolidBodyMotionFvMesh::update()
{
    pointField newPoints(points().size(), vector::zero);

    for (label regionI = 0; regionI < rs_.nRegions(); regionI++)
    {
        newPoints += motionMask(regionI)*
            transform
            (
                regionSBMFs_[regionI].transformation(),
                undisplacedPoints_
            );
    }

    fvMesh::movePoints(newPoints);

    return false;
}


// ************************************************************************* //
