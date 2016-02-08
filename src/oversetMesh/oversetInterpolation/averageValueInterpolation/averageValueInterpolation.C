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

#include "averageValueInterpolation.H"
#include "oversetInterpolation.H"
#include "oversetMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(averageValueInterpolation, 0);
    addToRunTimeSelectionTable
    (
        oversetInterpolation,
        averageValueInterpolation,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::averageValueInterpolation::calcWeights() const
{
    if (localWeightsPtr_ || remoteWeightsPtr_)
    {
        FatalErrorIn
        (
            "void averageValueInterpolation::calcWeights() const"
        )   << "Weights already calculated."
            << abort(FatalError);
    }

    // Get list of local donors
    const labelList& ld = overset().localDonors();

    // Get list of neighbouring donors (needed to set appropriate size of each
    // bottom most list of weights - for all donors)
    const labelListList& lnd = overset().localNeighbouringDonors();

    // Create a local weight field
    localWeightsPtr_ = new ScalarFieldField(ld.size());
    ScalarFieldField& localWeights = *localWeightsPtr_;

    // Loop through local donors, setting the weights
    forAll (ld, ldI)
    {
        // Get number of donors (master donor (1) + neighbouring donors)
        const label nDonors = 1 + lnd[ldI].size();

        // Calculate the weight field: each donor has the same influence as the
        // other (1/nDonors)
        localWeights.set
        (
            ldI,
            new scalarField
            (
                nDonors,
                1/scalar(nDonors)
            )
        );
    }

    // Handling remote donors for a parallel run
    if (Pstream::parRun())
    {
        // Get remote donor addressing
        const labelList& rd = overset().remoteDonors();
        const labelListList& rnd = overset().remoteNeighbouringDonors();

        // Create a global weights list for remote donor weights
        remoteWeightsPtr_ = new ListScalarFieldField(Pstream::nProcs());
        ListScalarFieldField& remoteWeights = *remoteWeightsPtr_;

        // Set the size of the weight field for this processor
        ScalarFieldField& myProcRemoteWeights =
            remoteWeights[Pstream::myProcNo()];
        myProcRemoteWeights.setSize(rd.size());

        // Loop through remote donors and calculate weights
        forAll (rd, rdI)
        {
            // Get number of donors (master donor (1) + neighbouring donors)
            const label nDonors = 1 + rnd[rdI].size();

            // Calculate the weight field: each donor has the same influence ast
            // the others (1/nDonors)
            myProcRemoteWeights.set
            (
                rdI,
                new scalarField
                (
                    nDonors,
                    1/scalar(nDonors)
                )
            );
        }

        // Gather remote weights (no need to scatter since the data is needed
        // only for the master processor).
        Pstream::gatherList(remoteWeights);
    }
}


void Foam::averageValueInterpolation::clearWeights() const
{
    deleteDemandDrivenData(localWeightsPtr_);
    deleteDemandDrivenData(remoteWeightsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::averageValueInterpolation::averageValueInterpolation
(
    const oversetMesh& overset,
    const dictionary& dict
)
:
    oversetInterpolation(overset, dict),
    localWeightsPtr_(NULL),
    remoteWeightsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::averageValueInterpolation::~averageValueInterpolation()
{
    clearWeights();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::oversetInterpolation::ScalarFieldField&
Foam::averageValueInterpolation::localWeights() const
{
    if (!localWeightsPtr_)
    {
        calcWeights();
    }

    return *localWeightsPtr_;
}


const Foam::oversetInterpolation::ListScalarFieldField&
Foam::averageValueInterpolation::remoteWeights() const
{
    // We cannot calculate the remoteWeights using usual lazy evaluation
    // mechanism since the data only exists on the master processor. Add
    // additional guards, VV, 8/Feb/2016.
    if (!Pstream::parRun())
    {
        FatalErrorIn
        (
            "const oversetInterpolation::ListScalarFieldField&\n"
            "averageValueInterpolation::remoteWeights() const"
        )   << "Attempted to calculate remoteWeights for a serial run."
            << "This is not allowed."
            << abort(FatalError);
    }
    else if (!Pstream::master())
    {
        FatalErrorIn
        (
            "const oversetInterpolation::ListScalarFieldField&\n"
            "averageValueInterpolation::remoteWeights() const"
        )   << "Attempted to calculate remoteWeights for a slave processor. "
            << "This is not allowed."
            << abort(FatalError);
    }
    else if (!remoteWeightsPtr_)
    {
        FatalErrorIn
        (
            "const oversetInterpolation::ListScalarFieldField&\n"
            "averageValueInterpolation::remoteWeights() const"
        )   << "Calculation of remoteWeights not possible because the data \n"
            << "exists only on the master processor. Please calculate \n"
            << "localWeights first (call .localWeights() member function)."
            << abort(FatalError);
    }

    return *remoteWeightsPtr_;
}


void Foam::averageValueInterpolation::update() const
{
    Info<< "averageValueInterpolation::update()" << endl;

    clearWeights();

    // Clear remote acceptor cell centres in the base class
    oversetInterpolation::update();
}


// ************************************************************************* //
