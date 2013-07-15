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

#include "compositeFringe.H"
#include "oversetFringe.H"
#include "Time.H"
#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(compositeFringe, 0);
    addToRunTimeSelectionTable(oversetFringe, compositeFringe, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::compositeFringe::readBaseFringes(const dictionary& dict)
{
    if (!baseFringes_.empty())
    {
        baseFringes_.clear();
    }

    PtrList<entry> baseFringeEntries(dict.lookup("baseFringes"));
    baseFringes_.setSize(baseFringeEntries.size());

    forAll (baseFringes_, bfI)
    {
        baseFringes_.set
        (
            bfI,
            oversetFringe::New
            (
                this->mesh(),
                this->region(),
                baseFringeEntries[bfI].dict()
            )
        );
    }

}

void Foam::compositeFringe::calcAddressing() const
{
    if (holesPtr_ || acceptorsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::compositeFringe::calcAddressing() const"
        )   << "Fringe addressing already calculated"
            << abort(FatalError);
    }

    // Initialise hole and acceptor look-up
    boolList holes(mesh().nCells(), false);
    boolList acceptors(mesh().nCells(), false);

    // Go through all base fringes and record holes and acceptors
    forAll (baseFringes_, bfI)
    {
        const labelList& ch = baseFringes_[bfI].holes();

        forAll (ch, chI)
        {
            holes[ch[chI]] = true;
        }

        const labelList& ca = baseFringes_[bfI].acceptors();

        forAll (ca, caI)
        {
            acceptors[ca[caI]] = true;
        }
    }


    // Count and collect holes
    {
        label nHoles = 0;

        forAll (holes, hI)
        {
            if (holes[hI])
            {
                nHoles++;
            }
        }

        holesPtr_ = new labelList(nHoles);
        labelList& h = *holesPtr_;
        nHoles = 0;

        forAll (holes, hI)
        {
            if (holes[hI])
            {
                h[nHoles] = hI;
                nHoles++;
            }
        }
    }

    // Count and collect acceptors
    {
        label nAcceptors = 0;

        forAll (acceptors, aI)
        {
            if (acceptors[aI])
            {
                nAcceptors++;
            }
        }

        acceptorsPtr_ = new labelList(nAcceptors);
        labelList& a = *acceptorsPtr_;
        nAcceptors = 0;

        forAll (acceptors, aI)
        {
            if (acceptors[aI])
            {
                a[nAcceptors] = aI;
                nAcceptors++;
            }
        }
    }
}


void Foam::compositeFringe::clearAddressing()
{
    deleteDemandDrivenData(holesPtr_);
    deleteDemandDrivenData(acceptorsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::compositeFringe::compositeFringe
(
    const fvMesh& mesh,
    const oversetRegion& region,
    const dictionary& dict
)
:
    oversetFringe(mesh, region, dict),
    holesPtr_(NULL),
    acceptorsPtr_(NULL)
{
    // Read fringes
    readBaseFringes(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::compositeFringe::~compositeFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::compositeFringe::holes() const
{
    if (!holesPtr_)
    {
        calcAddressing();
    }

    return *holesPtr_;
}


const Foam::labelList& Foam::compositeFringe::acceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }

    return *acceptorsPtr_;
}


void Foam::compositeFringe::update()
{
    Info<< "compositeFringe::update()" << endl;
    forAll (baseFringes_, bfI)
    {
        baseFringes_[bfI].update();
    }
}


// ************************************************************************* //
