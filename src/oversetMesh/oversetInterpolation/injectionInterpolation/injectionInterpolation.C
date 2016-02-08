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

void Foam::injectionInterpolation::calcWeights() const
{
    if (localWeightsPtr_ || remoteWeightsPtr_)
    {
        FatalErrorIn
        (
            "void injectionInterpolation::calcWeights() const"
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
        // Set the size of this donor weight field (master donor (1) +
        // neighbouring donors). Set the size and all values to zero.
        localWeights.set
        (
            ldI,
            new scalarField(1 + lnd[ldI].size(), 0)
        );

        // Injection interpolation: set the first value to 1 (master donor
        // contribution), others (neighbouring donors) are 0.
        localWeights[ldI][0] = 1;
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
            // Allocate the storage for this donor weight field (master donor
            // (1) + neighbouring donors). Set the size and all values to zero.
            myProcRemoteWeights.set
            (
                rdI,
                new scalarField(1 + rnd[rdI].size(), 0)
            );

            // Injection interpolation: set the first value to 1 (master donor
            // contribution), others (neighbouring donors) are 0.
            myProcRemoteWeights[rdI][0] = 1;
        }

        // Gather remote weights (no need to scatter since the data is needed
        // only for the master processor).
        Pstream::gatherList(remoteWeights);
    }
}


void Foam::injectionInterpolation::clearWeights() const
{
    deleteDemandDrivenData(localWeightsPtr_);
    deleteDemandDrivenData(remoteWeightsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::injectionInterpolation::injectionInterpolation
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

Foam::injectionInterpolation::~injectionInterpolation()
{
    clearWeights();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::oversetInterpolation::ScalarFieldField&
Foam::injectionInterpolation::localWeights() const
{
    if (!localWeightsPtr_)
    {
        calcWeights();
    }

    return *localWeightsPtr_;
}


const Foam::oversetInterpolation::ListScalarFieldField&
Foam::injectionInterpolation::remoteWeights() const
{
    // We cannot calculate the remoteWeights using usual lazy evaluation
    // mechanism since the data only exists on the master processor. Add
    // additional guards, VV, 8/Feb/2016.
    if (!Pstream::parRun())
    {
        FatalErrorIn
        (
            "const oversetInterpolation::ListScalarFieldField&\n"
            "injectionInterpolation::remoteWeights() const"
        )   << "Attempted to calculate remoteWeights for a serial run."
            << "This is not allowed."
            << abort(FatalError);
    }
    else if (!Pstream::master())
    {
        FatalErrorIn
        (
            "const oversetInterpolation::ListScalarFieldField&\n"
            "injectionInterpolation::remoteWeights() const"
        )   << "Attempted to calculate remoteWeights for a slave processor. "
            << "This is not allowed."
            << abort(FatalError);
    }
    else if (!remoteWeightsPtr_)
    {
        FatalErrorIn
        (
            "const oversetInterpolation::ListScalarFieldField&\n"
            "injectionInterpolation::remoteWeights() const"
        )   << "Calculation of remoteWeights not possible because the data \n"
            << "exists only on the master processor. Please calculate \n"
            << "localWeights first (call .localWeights() member function)."
            << abort(FatalError);
    }

    return *remoteWeightsPtr_;
}


void Foam::injectionInterpolation::update() const
{
    Info<< "injectionInterpolation::update()" << endl;

    clearWeights();
}


// ************************************************************************* //
