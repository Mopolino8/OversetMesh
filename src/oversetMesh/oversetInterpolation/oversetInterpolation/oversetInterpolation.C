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

#include "oversetInterpolation.H"
#include "oversetMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetInterpolation, 0);
    defineRunTimeSelectionTable(oversetInterpolation, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::oversetInterpolation::calcRemoteAccCC() const
{
    // It is an error to attempt to calculate remote acceptors in a serial run
    if (!Pstream::parRun())
    {
        FatalErrorIn
        (
            "void oversetInterpolation::calcRemoteAccCC() const"
        )   << "Attempted to calculate remote acceptor cell centres for a"
               " serial run."
            << abort(FatalError);
    }
    else
    {
        // Algorithm:
        // 1. Loop through all overset regions and collect donorAcceptor data
        //    (where donor is on a local processor) where acceptors may be local
        //    or remote.
        //    Note: we need to loop in the exactly the same way as in
        //    oversetMesh::calcParallelAddressing() member function to get the
        //    correct addressing (subscripting by remote donor index).
        // 2. If the acceptor is remote, store acceptor cell centre in the
        //    remote acceptor cell centre list, indexed by remote donor index.
        // 3. Finally, gather/scatter the list across all processors.

        // Create remote (global) acceptor cell centres indexed by local
        // donor index on each processor
        remoteAccCCPtr_ = new ListVectorField(Pstream::nProcs());
        ListVectorField& remAccCC = *remoteAccCCPtr_;

        // Get the list for this processor
        vectorField& procRemAccCC = remAccCC[Pstream::myProcNo()];

        // Set the size of field for this processor to number of remote donors
        procRemAccCC.setSize(overset().remoteDonors().size());

        // Get overset regions
        const PtrList<oversetRegion>& regions = overset().regions();

        // Set the counter for remote donors
        label nRemoteDonors = 0;

        // Loop through all regions
        forAll (regions, regionI)
        {
            // Get donor list for this region
            const donorAcceptorList& curDonors = regions[regionI].donors();

            // Loop through donors on this processor
            forAll (curDonors, dI)
            {
                if (curDonors[dI].acceptorProcNo() != Pstream::myProcNo())
                {
                    // This is a remote donor (donor is on local processor and
                    // acceptor is on a remote processor), store the value of
                    // its acceptor cell centre
                    procRemAccCC[nRemoteDonors] = curDonors[dI].acceptorPoint();

                    // Increment the counter of remote donors
                    ++nRemoteDonors;
                }
            }
        }

        // Check if the number of remote donors counted here and the one counted
        // in oversetMesh::calcParallelAddressing() is the same
        if (nRemoteDonors != overset().remoteDonors().size())
        {
            FatalErrorIn
            (
                "void Foam::oversetInterpolation::calcRemoteAccCC() const"
            )   << "Number of remote donors while storing remote acceptor cell"
                << " centres: " << nRemoteDonors << nl
                << "Number of remote donors as calculated in"
                << " oversetMesh::calcParallelAddressing(): "
                << overset().remoteDonors().size() << nl
                << "Two numbers are different, something went wrong with"
                << " algorithm for storing acceptor cell centres in parallel."
                << abort(FatalError);
        }

        // Each processor has taken care of its own remoteDonors, storing the
        // values of acceptor cell centres indexed by remoteDonors. Perform
        // gather/scatter to distribute the data.
        Pstream::gatherList(remAccCC);
        Pstream::scatterList(remAccCC);
    }
}


void Foam::oversetInterpolation::clearRemoteAccCC() const
{
    deleteDemandDrivenData(remoteAccCCPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetInterpolation::oversetInterpolation
(
    const oversetMesh& overset,
    const dictionary& dict
)
:
    overset_(overset),
    extensionLevel_(dict.lookupOrDefault<label>("extensionLevel", 0)),
    remoteAccCCPtr_(NULL)
{
    // Sanity check
    if (extensionLevel_ < 0)
    {
        FatalErrorIn
        (
            "oversetInterpolation::oversetInterpolation"
        ) << "Negative extension level not allowed for overset interpolation."
          << nl << "Use 0 for compact stencil and above 0 for expanded stencil."
          << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::oversetInterpolation::~oversetInterpolation()
{
    clearRemoteAccCC();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::oversetInterpolation::ListVectorField&
Foam::oversetInterpolation::remoteAccCC() const
{
    if (!Pstream::parRun())
    {
        FatalErrorIn
        (
            "const oversetInterpolation::ListVectorField&\n"
            "oversetInterpolation::remoteAccCC() const"
        )   << "Attempted to calculate remote acceptor cell centres for a "
            << "serial run. This is not allowed."
            << abort(FatalError);
    }
    else if (!remoteAccCCPtr_)
    {
        calcRemoteAccCC();
    }

    return *remoteAccCCPtr_;
}


void Foam::oversetInterpolation::update() const
{
    Info<< "oversetInterpolation::update()" << endl;

    clearRemoteAccCC();
}


// ************************************************************************* //
