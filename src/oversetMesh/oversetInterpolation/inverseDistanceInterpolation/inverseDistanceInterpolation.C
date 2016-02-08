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

#include "inverseDistanceInterpolation.H"
#include "oversetInterpolation.H"
#include "oversetMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(inverseDistanceInterpolation, 0);
    addToRunTimeSelectionTable
    (
        oversetInterpolation,
        inverseDistanceInterpolation,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::inverseDistanceInterpolation::calcWeights() const
{
    if (localWeightsPtr_ || remoteWeightsPtr_)
    {
        FatalErrorIn
        (
            "void inverseDistanceInterpolation::calcWeights() const"
        )   << "Weights already calculated."
            << abort(FatalError);
    }

    // Get list of local donors
    const labelList& ld = overset().localDonors();

    // Get list of neighbouring donors (needed to set appropriate size of each
    // bottom most list of weights - for all donors)
    const labelListList& lnd = overset().localNeighbouringDonors();

    // Get local donor addressing (local acceptor index for each local donor)
    const labelList& ldAddr = overset().localDonorAddr();

    // Get local acceptors and donors
    const labelList& acceptorCells = overset().acceptorCells();

    // Get necessary mesh data
    const vectorField& CC = overset().mesh().cellCentres();

    // Create a local weight field
    localWeightsPtr_ = new ScalarFieldField(ld.size());
    ScalarFieldField& localWeights = *localWeightsPtr_;

    // Loop through local donors, setting the weights
    forAll (ld, ldI)
    {
        // Set the size of this donor weight field (master donor (1) +
        // neighbouring donors). Set the size and all values to 1.
        localWeights.set
        (
            ldI,
            new scalarField(1 + lnd[ldI].size(), 1)
        );

        // Get acceptor cell centre
        const vector& accCC = CC[acceptorCells[ldAddr[ldI]]];

        // Inverse distance interpolation: weights are defined as inverted
        // distance from a given donor to the acceptor cell centre

        // Get reference to current weight field for this acceptor
        scalarField& curLocWeights = localWeights[ldI];

        // Calculate master donor first
        curLocWeights[0] /= mag(accCC - CC[ld[ldI]]) + SMALL;

        // Calculate neighbouring donors next
        const labelList& curNbrDonors = lnd[ldI];
        forAll (curNbrDonors, nbrI)
        {
            // Note nbrI + 1 index for subscripting because the first entry in
            // weights corresponds to the master donor
            curLocWeights[nbrI + 1] /=
                mag(accCC - CC[curNbrDonors[nbrI]]) + SMALL;
        }

        // Renormalize the weight field
        curLocWeights /= sum(curLocWeights);
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

        // Get the corresponding acceptor cell centres for remote donors on this
        // processor
        const vectorField& procRemAccCC = remoteAccCC()[Pstream::myProcNo()];

        // Set the size of the weight field for this processor
        ScalarFieldField& myProcRemoteWeights =
            remoteWeights[Pstream::myProcNo()];
        myProcRemoteWeights.setSize(rd.size());

        // Loop through remote donors and calculate weights
        forAll (rd, rdI)
        {
            // Allocate the storage for this donor weight field (master donor
            // (1) + neighbouring donors). Set the size and all values to 1.
            myProcRemoteWeights.set
            (
                rdI,
                new scalarField(1 + rnd[rdI].size(), 1)
            );

            // Get acceptor cell centre
            const vector& accCC = procRemAccCC[rdI];

            // Get reference to current weight field for this acceptor
            scalarField& curRemWeights = myProcRemoteWeights[rdI];

            // Calculate master donor first
            curRemWeights[0] /= mag(accCC - CC[rd[rdI]]);

            // Calculate neighbouring donors next
            const labelList& curNbrDonors = rnd[rdI];
            forAll (curNbrDonors, nbrI)
            {
                // Note nbrI + 1 index for subscripting because the first entry
                // in weights corresponds to the master donor
                curRemWeights[nbrI + 1] /= mag(accCC - CC[curNbrDonors[nbrI]]);
            }

            // Renormalize the weight field
            curRemWeights /= sum(curRemWeights);
        }

        // Gather remote weights (no need to scatter since the data is needed
        // only for the master processor).
        Pstream::gatherList(remoteWeights);
    }
}


void Foam::inverseDistanceInterpolation::clearWeights() const
{
    deleteDemandDrivenData(localWeightsPtr_);
    deleteDemandDrivenData(remoteWeightsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inverseDistanceInterpolation::inverseDistanceInterpolation
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

Foam::inverseDistanceInterpolation::~inverseDistanceInterpolation()
{
    clearWeights();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::oversetInterpolation::ScalarFieldField&
Foam::inverseDistanceInterpolation::localWeights() const
{
    if (!localWeightsPtr_)
    {
        calcWeights();
    }

    return *localWeightsPtr_;
}


const Foam::oversetInterpolation::ListScalarFieldField&
Foam::inverseDistanceInterpolation::remoteWeights() const
{
    // We cannot calculate the remoteWeights using usual lazy evaluation
    // mechanism since the data only exists on the master processor. Add
    // additional guards, VV, 8/Feb/2016.
    if (!Pstream::parRun())
    {
        FatalErrorIn
        (
            "const oversetInterpolation::ListScalarFieldField&\n"
            "inverseDistanceInterpolation::remoteWeights() const"
        )   << "Attempted to calculate remoteWeights for a serial run."
            << "This is not allowed."
            << abort(FatalError);
    }
    else if (!Pstream::master())
    {
        FatalErrorIn
        (
            "const oversetInterpolation::ListScalarFieldField&\n"
            "inverseDistanceInterpolation::remoteWeights() const"
        )   << "Attempted to calculate remoteWeights for a slave processor. "
            << "This is not allowed."
            << abort(FatalError);
    }
    else if (!remoteWeightsPtr_)
    {
        FatalErrorIn
        (
            "const oversetInterpolation::ListScalarFieldField&\n"
            "inverseDistanceInterpolation::remoteWeights() const"
        )   << "Calculation of remoteWeights not possible because the data \n"
            << "exists only on the master processor. Please calculate \n"
            << "localWeights first (call .localWeights() member function)."
            << abort(FatalError);
    }

    return *remoteWeightsPtr_;
}


void Foam::inverseDistanceInterpolation::update() const
{
    Info<< "inverseDistanceInterpolation::update()" << endl;

    clearWeights();

    // Clear remote acceptor cell centres in the base class
    oversetInterpolation::update();
}


// ************************************************************************* //
