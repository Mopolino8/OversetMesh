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
void Foam::oversetMesh::filterDonors
(
    Field<Type>& donorF,
    const Field<Type>& cellF
) const
{
    // Get donor list from interpolation
    const labelList& d = interpolation().donors();

    // Check sizes
    if
    (
        cellF.size() != mesh().nCells()
     || donorF.size() != d.size()
    )
    {
        FatalErrorIn
        (
            "void oversetMesh::filterDonors\n"
            "(\n"
            "    Field<Type>& donorF,\n"
            "    const Field<Type>& cellF,\n"
            ") const"
        )   << "Size of fields does not correspond to filtering" << nl
            << "Source field size = " << cellF.size()
            << " mesh size = " << mesh().nCells() << nl
            << " target field size = " << donorF.size()
            << " donor list size = " << d.size()
            << abort(FatalError);
    }

    // Collect donor cells
    forAll (donorF, i)
    {
        donorF[i] = cellF[d[i]];
    }
}


template<class Type>
void Foam::oversetMesh::donorToAcceptor
(
    Field<Type>& accF,
    const Field<Type>& donorF
) const
{
    // Interpolation with weights and addressing provided

    // Get donor list from interpolation.  Checking
    const labelList& d = interpolation().donors();

    // Get addressing and weights from interpolation
    const labelListList& addr = interpolation().addressing();
    const FieldField<Field, scalar>& weights = interpolation().weights();

    // Check sizes
    if
    (
        donorF.size() != d.size()
     || accF.size() != addr.size()
    )
    {
        FatalErrorIn
        (
            "void oversetMesh::donorToAcceptor\n"
            "(\n"
            "    Field<Type>& accF,\n"
            "    const Field<Type>& donorF,\n"
            ") const"
        )   << "Size of fields does not correspond to interpolation" << nl
            << "Source field size = " << donorF.size()
            << " donor list size = " << d.size()
            << " target field size = " << accF.size()
            << " acceptor list size = " << this->acceptorCells().size()
            << abort(FatalError);
    }

    forAll (accF, accI)
    {
        const labelList& nbr = addr[accI];
        const scalarField& w = weights[accI];
            
        accF[accI] = pTraits<Type>::zero;

        forAll (nbr, nI)
        {
            accF[accI] += donorF[nbr[nI]]*w[nI];
        }
    }
}


template<class Type>
void Foam::oversetMesh::interpolate
(
    Field<Type>& accF,
    const Field<Type>& cellF
) const
{
    Field<Type> donorF(interpolation().donors().size());

    filterDonors(donorF, cellF);

    donorToAcceptor(accF, donorF);
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
