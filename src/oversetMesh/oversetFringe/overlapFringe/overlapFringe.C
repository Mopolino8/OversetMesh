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
#include "oversetFringe.H"
#include "polyPatchID.H"
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

    HJ, HERE!!!

    // Find patches
    labelList cellType(mesh().nCells(), ACTIVE);
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
