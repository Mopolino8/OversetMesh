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

#include "manualFringe.H"
#include "oversetFringe.H"
#include "Time.H"
#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(manualFringe, 0);
    addToRunTimeSelectionTable(oversetFringe, manualFringe, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::manualFringe::calcAddressing() const
{
    if (holesPtr_ || acceptorsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::manualFringe::calcAddressing() const"
        )   << "Fringe addressing already calculated"
            << abort(FatalError);
    }

    holesPtr_ = new labelList
    (
        cellSet
        (
            mesh(),
            holesSetName_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ).toc()
    );

    acceptorsPtr_ = new labelList
    (
        cellSet
        (
            mesh(),
            acceptorsSetName_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ).toc()
    );
}


void Foam::manualFringe::clearAddressing()
{
    deleteDemandDrivenData(holesPtr_);
    deleteDemandDrivenData(acceptorsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::manualFringe::manualFringe
(
    const fvMesh& mesh,
    const oversetRegion& region,
    const dictionary& dict
)
:
    oversetFringe(mesh, region, dict),
    holesSetName_(dict.lookup("holes")),
    acceptorsSetName_(dict.lookup("acceptors")),
    holesPtr_(NULL),
    acceptorsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::manualFringe::~manualFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::manualFringe::holes() const
{
    if (!holesPtr_)
    {
        calcAddressing();
    }

    return *holesPtr_;
}


const Foam::labelList& Foam::manualFringe::acceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }

    return *acceptorsPtr_;
}


void Foam::manualFringe::update()
{
    Info<< "manualFringe::update()" << endl;
}


// ************************************************************************* //
