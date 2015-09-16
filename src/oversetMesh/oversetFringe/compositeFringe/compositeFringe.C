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
#include "oversetRegion.H"
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
    if (fringeHolesPtr_ || acceptorsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::compositeFringe::calcAddressing() const"
        )   << "Fringe addressing already calculated"
            << abort(FatalError);
    }

    // Get reference to region cell zone
    const cellZone& rcz = region().zone();

    // Make a hash set to collect acceptor points
    labelHashSet acceptorSet;

    // Make a hash set to collect hole points
    labelHashSet fringeHoleSet;

    // Go through all base fringes and record holes and acceptors
    forAll (baseFringes_, bfI)
    {
        const labelList& ch = baseFringes_[bfI].fringeHoles();

        forAll (ch, chI)
        {
            // Check if cell is in region zone
            if (rcz.whichCell(ch[chI]) > -1)
            {
                // Found acceptor
                acceptorSet.insert(ch[chI]);
            }
        }

        const labelList& ca = baseFringes_[bfI].acceptors();

        forAll (ca, caI)
        {
            // Check if cell is in region zone
            if (rcz.whichCell(ca[caI]) > -1)
            {
                // Found acceptor
                acceptorSet.insert(ca[caI]);
            }
        }
    }

    // Collect fringeHoles
    fringeHolesPtr_ = new labelList(fringeHoleSet.sortedToc());

    // Collect acceptors
    acceptorsPtr_ = new labelList(acceptorSet.sortedToc());
}


void Foam::compositeFringe::clearAddressing()
{
    deleteDemandDrivenData(fringeHolesPtr_);
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
    fringeHolesPtr_(NULL),
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

const Foam::labelList& Foam::compositeFringe::fringeHoles() const
{
    if (!fringeHolesPtr_)
    {
        calcAddressing();
    }

    return *fringeHolesPtr_;
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
    forAll (baseFringes_, bfI)
    {
        baseFringes_[bfI].update();
    }
}


// ************************************************************************* //
