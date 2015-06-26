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

#include "faceCellsFringe.H"
#include "oversetRegion.H"
#include "polyPatchID.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceCellsFringe, 0);
    addToRunTimeSelectionTable(oversetFringe, faceCellsFringe, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceCellsFringe::calcAddressing() const
{
    if (acceptorsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::faceCellsFringe::calcAddressing() const"
        )   << "Addressing already calculated"
            << abort(FatalError);
    }

    // Get reference to region cell zone
    const cellZone& rcz = region().zone();

    // Make a hash set to collect acceptor points
    // Note: 2 patches connecting at the corner may create a duplicate,
    // which is filtered on insertion
    labelHashSet acceptorSet;

    // Find patches and mark cells
    forAll (patchNames_, nameI)
    {
        polyPatchID curFringePatch
        (
            patchNames_[nameI],
            mesh().boundaryMesh()
        );

        if (!curFringePatch.active())
        {
            FatalErrorIn
            (
                "void faceCellsFringe::calcAddressing() const"
            )   << "Fringe patch " << patchNames_[nameI]
                << " cannot be found"
                << abort(FatalError);
        }

        const unallocLabelList& curFaceCells =
            mesh().boundaryMesh()[curFringePatch.index()].faceCells();

        forAll (curFaceCells, fcI)
        {
            // Check if cell is in region zone
            if (rcz.whichCell(curFaceCells[fcI]) > -1)
            {
                // Found acceptor
                acceptorSet.insert(curFaceCells[fcI]);
            }
        }
    }

    // Collect acceptors
    acceptorsPtr_ = new labelList(acceptorSet.sortedToc());

    // Holes currently empty
    holesPtr_ = new labelList();
}


void Foam::faceCellsFringe::clearAddressing()
{
    deleteDemandDrivenData(holesPtr_);
    deleteDemandDrivenData(acceptorsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::faceCellsFringe::faceCellsFringe
(
    const fvMesh& mesh,
    const oversetRegion& region,
    const dictionary& dict
)
:
    oversetFringe(mesh, region, dict),
    patchNames_(dict.lookup("patches")),
    holesPtr_(NULL),
    acceptorsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceCellsFringe::~faceCellsFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::faceCellsFringe::holes() const
{
    if (!holesPtr_)
    {
        calcAddressing();
    }

    return *holesPtr_;
}


const Foam::labelList& Foam::faceCellsFringe::acceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }

    return *acceptorsPtr_;
}


void Foam::faceCellsFringe::update()
{
    Info<< "faceCellsFringe::update()" << endl;
}


// ************************************************************************* //
