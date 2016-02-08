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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "oversetMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::oversetMesh::interpolate
(
    Field<Type>& accF,
    const Field<Type>& cellF
) const
{
    // Check sizes
    if
    (
        accF.size() != this->acceptorCells().size()
     || cellF.size() != this->mesh().nCells()
    )
    {
        FatalErrorIn
        (
            "void oversetMesh::donorToAcceptor\n"
            "(\n"
            "    Field<Type>& accF,\n"
            "    const Field<Type>& ,cellF\n"
            ") const"
        )   << "Size of fields does not correspond to interpolation" << nl
            << "Source field size = " << cellF.size()
            << " mesh size = " << this->mesh().nCells()
            << " target field size = " << accF.size()
            << " acceptor list size = " << this->acceptorCells().size()
            << abort(FatalError);
    }

    // Get local interpolation weights
    // For each acceptor cell, weights are stored in the following fashion:
    // - zeroth entry is the weight corresponding to the master donor cell
    // - other entries are weights corresponding to neighbouring donor cells
    // Hence, each acceptor cell has 1 + n values, where n is the number of
    // eligible neighbours of the master donor cell, and the 1 comes from the
    // master donor cell itself.
    const oversetInterpolation::ScalarFieldField& localWeights =
        interpolation().localWeights();

    // Get list of local donors
    const labelList& ld = this->localDonors();

    // Get local neighbouring donors
    const labelListList& lnd = this->localNeighbouringDonors();

    // Get local donor addressing
    const labelList& ldAddr = this->localDonorAddr();

    // Loop through local donors
    forAll (ld, ldI)
    {
        // Get weights for donor cells of this acceptor
        const scalarField& donorWeights = localWeights[ldI];

        // Get reference to the current acceptor value
        Type& curAccF = accF[ldAddr[ldI]];

        // First set master donor contribution
        curAccF = donorWeights[0]*cellF[ld[ldI]];

        // Add contributions from neighbouring donors
        const labelList& curNbrDonors = lnd[ldI];
        forAll (curNbrDonors, nbrI)
        {
            // Note nbrI + 1 index for subscripting because the first entry
            // in weights corresponds to the master donor
            curAccF += donorWeights[nbrI + 1]*cellF[curNbrDonors[nbrI]];
        }
    }

    if (Pstream::parRun())
    {
        // Get remote donor addressing
        const labelList& rd = this->remoteDonors();
        const labelListList& rnd = this->remoteNeighbouringDonors();

        // We will combine all the donor values (from master donor cell and
        // neighbouring donors) in a single list (similar to weights) in order
        // to have only one parallel communication instead of two
        List<FieldField<Field, Type> > globalRemoteDonors(Pstream::nProcs());

        // Fill in all remote donors (master + neighbouring donors)
        FieldField<Field, Type>& rdAll =
            globalRemoteDonors[Pstream::myProcNo()];

        // Set the size for each processor to correspond to the number of remote
        // donors (master donors)
        rdAll.setSize(rd.size());

        // For all master donors, allocate necessary storage for neighbouring
        // donors as well
        forAll (rdAll, rdI)
        {
            // Get a list of neighbouring donors
            const labelList& rndCur = rnd[rdI];

            // Allocate necessary storage, setting the size of this donor field
            // (master donor (1) + neighbouring donors) and initialising with
            // zero
            rdAll.set
            (
                rdI,
                new Field<Type>(1 + rndCur.size(), pTraits<Type>::zero)
            );

            // Get reference to all current remote donor values
            Field<Type>& rdAllCur = rdAll[rdI];

            // Populate the donor values to prepare for communication
            // First value is master donor value
            rdAllCur[0] = cellF[rd[rdI]];

            // Successive values are neighbouring donors for this master donor
            forAll (rndCur, nbrI)
            {
                rdAllCur[nbrI + 1] = cellF[rndCur[nbrI]];
            }
        }

        // Communicate to master
        Pstream::gatherList(globalRemoteDonors);

        // Prepare acceptor list
        List<List<Type> > globalRemoteAcceptors(Pstream::nProcs());

        // Master processor reorganises the donor data for each target
        if (Pstream::master())
        {
            // Get addressing
            const labelListList& globalAccProc = this->globalAcceptFromProc();
            const labelListList& globalAccCell = this->globalAcceptFromCell();

            // Resize and fill remote acceptor processor arrays for
            // all processors
            forAll (globalRemoteAcceptors, procI)
            {
                // Get processor addressing
                const labelList& procAccProc = globalAccProc[procI];
                const labelList& procAccCell = globalAccCell[procI];

                // Get processor acceptor data to fill in
                List<Type>& procRA = globalRemoteAcceptors[procI];

                // Resize results list
                procRA.setSize(procAccProc.size());

                // Get remote interpolation weights. For each processor, each
                // acceptor cell is associated with a master donor cell and its
                // eligible neighbours (live and other donor cells).
                const oversetInterpolation::ListScalarFieldField& remWeights =
                    interpolation().remoteWeights();

                forAll (procRA, accI)
                {
                    // Get remote processor index for remote donors
                    const label& donorProcIndex = procAccProc[accI];

                    // Get remote donor index
                    const label& donorIndex = procAccCell[accI];

                    // Get weights for this set of donors (associated with this
                    // acceptor)
                    const scalarField& donorWeights =
                        remWeights[donorProcIndex][donorIndex];

                    // Get all current remote donors
                    const Field<Type>& allDonorsForCurAcceptor =
                        globalRemoteDonors[donorProcIndex][donorIndex];

                    // Get reference of the current acceptor value and reset it
                    Type& curProcRA = procRA[accI];
                    curProcRA = pTraits<Type>::zero;

                    // Calculate the acceptor value from all weighted remote
                    // donor values
                    forAll (allDonorsForCurAcceptor, rdI)
                    {
                        curProcRA +=
                            donorWeights[rdI]*allDonorsForCurAcceptor[rdI];
                    }
                }
            }
        }

        // Communicate global acceptors to all processors
        Pstream::scatter(globalRemoteAcceptors);

        // Insert remote acceptors

        // Get data belonging to local processor
        const List<Type>& procRemoteAcceptors =
            globalRemoteAcceptors[Pstream::myProcNo()];

        // Get addressing
        const labelList& raAddr = this->remoteAcceptorAddr();

        // Insert remote acceptors
        forAll (raAddr, raI)
        {
            accF[raAddr[raI]] = procRemoteAcceptors[raI];
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::oversetMesh::interpolate
(
    const Field<Type>& cellF
) const
{
    tmp<Field<Type> > tresult(new Field<Type>(this->acceptorCells().size()));
    Field<Type>& result = tresult();

    interpolate(result, cellF);

    return tresult;
}


// ************************************************************************* //
