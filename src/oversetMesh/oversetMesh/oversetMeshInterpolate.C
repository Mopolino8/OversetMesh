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

    // Insert local donors into local acceptors

    // Get list of local donors
    const labelList& ld = this->localDonors();

    // Get local donor addressing
    const labelList& ldAddr = this->localDonorAddr();

    forAll (ld, ldI)
    {
        accF[ldAddr[ldI]] = cellF[ld[ldI]];
    }

    if (Pstream::parRun())
    {
        // Get remote donor addressing
        const labelList& rdAddr = this->remoteDonors();

        // Collect remote donors for communication
        List<List<Type> > globalRemoteDonors(Pstream::nProcs());

        // Fill in remote donors
        List<Type>& rd = globalRemoteDonors[Pstream::myProcNo()];
        rd.setSize(rdAddr.size());

        forAll (rdAddr, rdI)
        {
            rd[rdI] = cellF[rdAddr[rdI]];
        }

        // Communicate to master
        Pstream::gatherList(globalRemoteDonors);

        // Prepare acceptor list
        List<List<Type> > globalRemoteAcceptors(Pstream::nProcs());
    
        // Master processor reorganises the donor data for each targer
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

                forAll (procRA, accI)
                {
                    procRA[accI] =
                        globalRemoteDonors[procAccProc[accI]]
                        [procAccCell[accI]];
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
