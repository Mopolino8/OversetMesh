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

#include "injectionInterpolation.H"
#include "oversetInterpolation.H"
#include "oversetMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(injectionInterpolation, 0);
    addToRunTimeSelectionTable
    (
        oversetInterpolation,
        injectionInterpolation,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::injectionInterpolation::calcAddressing() const
{
    if (addressingPtr_ || weightsPtr_)
    {
        FatalErrorIn
        (
            "void injectionInterpolation::calcAddressing() const"
        )   << "Addressing already calculated"
            << abort(FatalError);
    }


    addressingPtr_ = new labelListList(overset().acceptorCells().size());
    labelListList& addr = *addressingPtr_;

    weightsPtr_ = new scalarListList(overset().acceptorCells().size());
    scalarListList& w = *weightsPtr_;

    const labelList& donors = overset().fringeAddressing();

    forAll (addr, addrI)
    {
        // Single donor addressing
        addr[addrI].setSize(1);
        addr[addrI][0] = donors[addrI];

        w[addrI].setSize(1);
        w[addrI][0] = 1;
    }
}


void Foam::injectionInterpolation::clearAddressing() const
{
    deleteDemandDrivenData(addressingPtr_);
    deleteDemandDrivenData(weightsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::injectionInterpolation::injectionInterpolation
(
    const oversetMesh& overset,
    const dictionary& dict
)
:
    oversetInterpolation(overset, dict),
    addressingPtr_(NULL),
    weightsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::injectionInterpolation::~injectionInterpolation()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::injectionInterpolation::donors() const
{
    // Compact addressing: single donor for single cell injection
    return overset().donorCells();
}


const Foam::labelListList& Foam::injectionInterpolation::addressing() const
{
    if (!addressingPtr_)
    {
        calcAddressing();
    }

    return *addressingPtr_;
}


const Foam::scalarListList& Foam::injectionInterpolation::weights() const
{
    if (!weightsPtr_)
    {
        calcAddressing();
    }

    return *weightsPtr_;
}


void Foam::injectionInterpolation::update()
{
    Info<< "injectionInterpolation::update()" << endl;
}


// ************************************************************************* //
