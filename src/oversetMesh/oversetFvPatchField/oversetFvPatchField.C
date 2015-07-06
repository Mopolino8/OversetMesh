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

#include "oversetFvPatchField.H"
#include "oversetFvPatch.H"
#include "fvMatrix.H"
#include "demandDrivenData.H"

#include "IndirectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class Type2>
void oversetFvPatchField<Type>::setHoleValues(Field<Type2>& f) const
{
    const labelList& dc = oversetPatch_.overset().holeCells();

    forAll (dc, dcI)
    {
        f[dc[dcI]] = holeCellValue_;
    }
}


template<class Type>
template<class Type2>
void oversetFvPatchField<Type>::setAcceptorValues(Field<Type2>& f) const
{
    // Get acceptor values by interpolation
    Field<Type2> accValues = oversetPatch_.overset().interpolate(f);

    // Ger acceptor addressing
    const labelList& acceptors =
        oversetPatch_.overset().acceptorCells();

    // Check sizes
    if (accValues.size() != acceptors.size())
    {
        FatalErrorIn
        (
            "oversetFvPatchField<Type>::"
            "setAcceptorValues(Field<Type2>& f) const"
        )   << "Bad sizes: " << accValues.size()
            << " and " << acceptors.size()
            << abort(FatalError);
    }

    forAll (acceptors, accI)
    {
        f[acceptors[accI]] = accValues[accI];
    }
}


template<class Type>
void oversetFvPatchField<Type>::correctDiag
(
    fvMatrix<Type>& eqn
) const
{
    // Get access to diagonal
    scalarField& diag = eqn.diag();
    Field<Type>& source = eqn.source();

    const labelList& holeCells = oversetPatch_.overset().holeCells();
    const labelList& acceptorCells = oversetPatch_.overset().acceptorCells();

    label nLiveCells = diag.size() - holeCells.size() - acceptorCells.size();

    reduce(nLiveCells, sumOp<label>());

    // Estimate diagonal in live cells
    scalar liveDiag = 1;

    if (nLiveCells > 0)
    {
        liveDiag = gSumMag(diag)/nLiveCells;

        // Correct for sign
        liveDiag *= sign(gMax(diag));
    }
    else
    {
        FatalErrorIn
        (
            "void oversetFvPatchField<Type>::correctDiag\n"
            "(\n"
            "    fvMatrix<Type>& eqn\n"
            ") const"
        )   << "No live cells in matrix"
            << abort(FatalError);
    }

    // Fix diagonal if missing in hole cells
    forAll (holeCells, hcI)
    {
        if (mag(diag[holeCells[hcI]]) < SMALL)
        {
            diag[holeCells[hcI]] = liveDiag;
        }
    }

    // Fix diagonal if missing in acceptor cells
    forAll (acceptorCells, acI)
    {
        if (mag(diag[acceptorCells[acI]]) < SMALL)
        {
            diag[acceptorCells[acI]] = liveDiag;
        }
    }

    // Fix source in acceptor cells
    forAll (acceptorCells, acI)
    {
        source[acceptorCells[acI]] = pTraits<Type>::zero;
    }
}


template<class Type>
void oversetFvPatchField<Type>::correctOffDiag
(
    fvMatrix<Type>& eqn
) const
{
    // Kill off-diagonal coefficient in all hole and acceptor faces
    // Collect off-diagonal coefficients for all fringe faces

    // 1 Hole internal faces
    const labelList& holeInternalFaces =
        oversetPatch_.overset().holeInternalFaces();

    // 2 Acceptor internal faces
    const labelList& acceptorInternalFaces =
        oversetPatch_.overset().acceptorInternalFaces();

    // 3 Hole faces
    const labelList& holeFaces = oversetPatch_.overset().holeFaces();

    // 4 Fringe faces
    const labelList& fringeFaces = oversetPatch_.overset().fringeFaces();
//     const boolList& fringeFaceFlips =
//         oversetPatch_.overset().fringeFaceFlips();
    
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();

    if (eqn.symmetric())
    {
        scalarField& upper = eqn.upper();

//         Info<< "Symmetric correctOffDiag for field "
//             << this->dimensionedInternalField().name() << endl;

        // 1 Hole internal faces
        forAll (holeInternalFaces, hifI)
        {
            const label curFace = holeInternalFaces[hifI];

            if (curFace < upper.size())
            {
                upper[curFace] = 0;
            }
            else
            {
                // Coupled boundary
                label patchi = mesh.boundaryMesh().whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchi].empty())
                {
                    label patchFacei =
                        mesh.boundaryMesh()[patchi].whichFace(curFace);

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;
                }
            }
        }

        // 2 Acceptor internal faces
        forAll (acceptorInternalFaces, hifI)
        {
            const label curFace = acceptorInternalFaces[hifI];

            if (curFace < upper.size())
            {
                upper[curFace] = 0;
            }
            else
            {
                // Coupled boundary
                label patchi = mesh.boundaryMesh().whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchi].empty())
                {
                    label patchFacei =
                        mesh.boundaryMesh()[patchi].whichFace(curFace);

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;
                }
            }
        }

        // 3 Hole faces
        forAll (holeFaces, faceI)
        {
            const label curFace = holeFaces[faceI];

            if (curFace < upper.size())
            {
                upper[curFace] = 0;
            }
            else
            {
                // Coupled boundary
                label patchi = mesh.boundaryMesh().whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchi].empty())
                {
                    label patchFacei =
                        mesh.boundaryMesh()[patchi].whichFace(curFace);

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;
                }
            }
        }

        // 4 Fringe faces
        fringeUpperCoeffs_.setSize(fringeFaces.size(), 0);
        fringeLowerCoeffs_.setSize(fringeFaces.size(), 0);

        forAll (fringeFaces, fringeI)
        {
            const label curFace = fringeFaces[fringeI];

            if (curFace < upper.size())
            {
                // Internal face: kill the off-diagonal coefficient in
                // the live cell

                // Since the matrix is symmetric, there is no need to
                // distinguish between lower and upper coefficient
                fringeUpperCoeffs_[fringeI] = upper[curFace];
                fringeLowerCoeffs_[fringeI] = upper[curFace];

                upper[curFace] = 0;
            }
            else
            {
                // Coupled boundary
                label patchi = mesh.boundaryMesh().whichPatch(curFace);

                if (psi.boundaryField()[patchi].coupled())
                {
                    // Note: on coupled boundaries, all coefficients are
                    // identical.  We can take the first component
                    // HJ, 30/May/2013

                    label patchFacei =
                        mesh.boundaryMesh()[patchi].whichFace(curFace);

                    // For a coupled boundary
                    fringeUpperCoeffs_[fringeI] =
                        component(eqn.boundaryCoeffs()[patchi][patchFacei], 0);

                    fringeLowerCoeffs_[fringeI] =
                        component(eqn.internalCoeffs()[patchi][patchFacei], 0);

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;
                }
            }
        }
    }
    else if (eqn.asymmetric())
    {
        scalarField& upper = eqn.upper();
        scalarField& lower = eqn.lower();

//         Info<< "Asymmetric correctOffDiag for field "
//             << this->dimensionedInternalField().name() << endl;

        // 1 Hole internal faces
        forAll (holeInternalFaces, hifI)
        {
            const label curFace = holeInternalFaces[hifI];

            if (curFace < upper.size())
            {
                upper[curFace] = 0;
                lower[curFace] = 0;
            }
            else
            {
                // Coupled boundary
                label patchi = mesh.boundaryMesh().whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchi].empty())
                {
                    label patchFacei =
                        mesh.boundaryMesh()[patchi].whichFace(curFace);

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;
                }
            }
        }

        // 2 Acceptor internal faces
        forAll (acceptorInternalFaces, hifI)
        {
            const label curFace = acceptorInternalFaces[hifI];

            if (curFace < upper.size())
            {
                upper[curFace] = 0;
                lower[curFace] = 0;
            }
            else
            {
                // Coupled boundary
                label patchi = mesh.boundaryMesh().whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchi].empty())
                {
                    label patchFacei =
                        mesh.boundaryMesh()[patchi].whichFace(curFace);

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;
                }
            }
        }

        // 3 Hole faces
        forAll (holeFaces, faceI)
        {
            const label curFace = holeFaces[faceI];

            if (curFace < upper.size())
            {
                upper[curFace] = 0;
                lower[curFace] = 0;
            }
            else
            {
                // Coupled boundary
                label patchi = mesh.boundaryMesh().whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchi].empty())
                {
                    label patchFacei =
                        mesh.boundaryMesh()[patchi].whichFace(curFace);

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;
                }
            }
        }

        // 4 Fringe faces
        fringeUpperCoeffs_.setSize(fringeFaces.size(), 0);
        fringeLowerCoeffs_.setSize(fringeFaces.size(), 0);

        forAll (fringeFaces, fringeI)
        {
            const label curFace = fringeFaces[fringeI];

            if (curFace < upper.size())
            {
                // Internal face: kill the off-diagonal coefficient in
                // the live cell

                // Since the matrix is symmetric, there is no need to
                // distinguish between lower and upper coefficient
                fringeUpperCoeffs_[fringeI] = upper[curFace];
                fringeLowerCoeffs_[fringeI] = lower[curFace];

                upper[curFace] = 0;
                lower[curFace] = 0;
            }
            else
            {
                // Coupled boundary
                label patchi = mesh.boundaryMesh().whichPatch(curFace);

                if (psi.boundaryField()[patchi].coupled())
                {
                    // Note: on coupled boundaries, all coefficients are
                    // identical.  We can take the first component
                    // HJ.30/May/2013

                    label patchFacei =
                        mesh.boundaryMesh()[patchi].whichFace(curFace);

                    // For a coupled boundary
                    fringeUpperCoeffs_[fringeI] =
                        component(eqn.boundaryCoeffs()[patchi][patchFacei], 0);

                    fringeLowerCoeffs_[fringeI] =
                        component(eqn.internalCoeffs()[patchi][patchFacei], 0);

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    oversetPatch_(refCast<const oversetFvPatch>(p)),
    coupledFringe_(false),
    setHoleCellValue_(false),
    holeCellValue_(pTraits<Type>::zero),
    fringeUpperCoeffs_(),
    fringeLowerCoeffs_()
{}


template<class Type>
oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvPatchField<Type>(p, iF, f),
    oversetPatch_(refCast<const oversetFvPatch>(p)),
    coupledFringe_(false),
    setHoleCellValue_(false),
    holeCellValue_(pTraits<Type>::zero),
    fringeUpperCoeffs_(),
    fringeLowerCoeffs_()
{}


template<class Type>
oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    oversetPatch_(refCast<const oversetFvPatch>(p)),
    coupledFringe_(dict.lookup("coupledFringe")),
    setHoleCellValue_(dict.lookup("setHoleCellValue")),
    holeCellValue_(pTraits<Type>(dict.lookup("holeCellValue"))),
    fringeUpperCoeffs_(),
    fringeLowerCoeffs_()
{
    if (!isType<oversetFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "oversetFvPatchField<Type>::oversetFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    oversetPatch_(refCast<const oversetFvPatch>(p)),
    coupledFringe_(ptf.coupledFringe_),
    setHoleCellValue_(ptf.setHoleCellValue_),
    holeCellValue_(ptf.holeCellValue_),
    fringeUpperCoeffs_(),
    fringeLowerCoeffs_()
{
    if (!isType<oversetFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "oversetFvPatchField<Type>::oversetFvPatchField\n"
            "(\n"
            "    const oversetFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf
)
:
    oversetLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    oversetPatch_(refCast<const oversetFvPatch>(ptf.patch())),
    coupledFringe_(ptf.coupledFringe_),
    setHoleCellValue_(ptf.setHoleCellValue_),
    holeCellValue_(ptf.holeCellValue_),
    fringeUpperCoeffs_(),
    fringeLowerCoeffs_()
{}


template<class Type>
oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    oversetPatch_(refCast<const oversetFvPatch>(ptf.patch())),
    coupledFringe_(ptf.coupledFringe_),
    setHoleCellValue_(ptf.setHoleCellValue_),
    holeCellValue_(ptf.holeCellValue_),
    fringeUpperCoeffs_(),
    fringeLowerCoeffs_()
{}


// * * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * //

template<class Type>
oversetFvPatchField<Type>::~oversetFvPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > oversetFvPatchField<Type>::patchNeighbourField() const
{
    // Warning: returning own patch field, as "neighbour" does not exist
    return *this;
}


template<class Type>
void oversetFvPatchField<Type>::updateCoeffs()
{
    // Fix the value in hole cells.  Probably unnecessary.  HJ, 25/Jun/2013
    if (setHoleCellValue_)
    {
        // Get non-const access to internal field
        Field<Type>& psiI = const_cast<Field<Type>&>(this->internalField());

        this->setHoleValues(psiI);
    }

    // Clear fringe coefficients
    fringeUpperCoeffs_.clear();
    fringeLowerCoeffs_.clear();

    fvPatchField<Type>::updateCoeffs();
}


template<class Type>
void oversetFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if (coupledFringe_)
    {
        // Get non-constant access to internal field
        Field<Type>& psi = const_cast<Field<Type>&>(this->internalField());

        // Set acceptor values
        this->setAcceptorValues(psi);
    }
}


template<class Type>
void oversetFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (setHoleCellValue_)
    {
        // Get non-constant access to internal field
        Field<Type>& psi = const_cast<Field<Type>&>(this->internalField());

        this->setHoleValues(psi);
    }
}


template<class Type>
void oversetFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& eqn
)
{
    // Eliminate unnecesary off-diagonal coefficients
    this->correctOffDiag(eqn);

    // Build matrix diagonal for cells where it is missing
    this->correctDiag(eqn);

    // Set values in hole cells
    const labelList& holeCells = oversetPatch_.overset().holeCells();

    // Correct equation for hole cells
    Field<Type> holeCellsPsi
    (
        holeCells.size(),
        holeCellValue_
    );

    eqn.setValues(holeCells, holeCellsPsi);
}


template<class Type>
void oversetFvPatchField<Type>::transformCoupleField
(
    scalarField& f,
    const direction cmpt
) const
{}


template<class Type>
void oversetFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix& m,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011
    if (coupledFringe_)
    {
        // Get non-const access to psi
        scalarField& psi = const_cast<scalarField&>(psiInternal);

        // Set acceptor values in field
        this->setAcceptorValues(psi);

        // Do multiplication on fringe faces

        const unallocLabelList& own = m.lduAddr().lowerAddr();
        const unallocLabelList& nei = m.lduAddr().upperAddr();

        const labelList& fringeFaces = oversetPatch_.overset().fringeFaces();
        const boolList& fringeFaceFlips =
            oversetPatch_.overset().fringeFaceFlips();

        if (switchToLhs)
        {
            forAll (fringeFaces, fringeI)
            {
                const label& curFace = fringeFaces[fringeI];

                // Multiplication is only done for internal fringe faces
                // HJ, 1/May/2015
                if (curFace < own.size())
                {
                    // Get addressing
                    const label& o = own[curFace];
                    const label& n = nei[curFace];

                    // Note change of sign in multiplication, because
                    // fringe coefficients belong to A
                    // HJ, 22/May/2013
                    if (fringeFaceFlips[fringeI])
                    {
                        // Face pointing into live cell
                        // Add overset off-diagonal contribution to live cell
                        result[n] -= fringeLowerCoeffs_[fringeI]*psi[o];
                    }
                    else
                    {
                        // Face pointing out of live cell
                        // Add overset off-diagonal contribution to live cell
                        result[o] -= fringeUpperCoeffs_[fringeI]*psi[n];
                    }
                }
            }
        }
        else
        {
            forAll (fringeFaces, fringeI)
            {
                const label& curFace = fringeFaces[fringeI];

                // Multiplication is only done for internal fringe faces
                // HJ, 1/May/2015
                if (curFace < own.size())
                {
                    // Get addressing
                    const label& o = own[curFace];
                    const label& n = nei[curFace];

                    // Note change of sign in multiplication, because
                    // fringe coefficients belong to A
                    // HJ, 22/May/2013
                    if (fringeFaceFlips[fringeI])
                    {
                        // Face pointing into of live cell
                        // Add overset off-diagonal contribution to live cell
                        result[n] += fringeLowerCoeffs_[fringeI]*psi[o];
                    }
                    else
                    {
                        // Face pointing out live cell
                        // Add overset off-diagonal contribution to live cell
                        result[o] += fringeUpperCoeffs_[fringeI]*psi[n];
                    }
                }
            }
        }

        // Do acceptor cells
        // Get diagonal
        const scalarField& diag = m.diag();

        const labelList& acceptorCells =
            oversetPatch_.overset().acceptorCells();

        forAll (acceptorCells, acI)
        {
            const label& curCell = acceptorCells[acI];

            result[curCell] = -diag[curCell]*psi[curCell];
        }
    }
    else
    {
        FatalErrorIn
        (
            "void oversetFvPatchField<Type>::updateInterfaceMatrix\n"
            "(\n"
            "    const scalarField& psiInternal,\n"
            "    scalarField& result,\n"
            "    const lduMatrix&,\n"
            "    const scalarField& coeffs,\n"
            "    const direction cmpt,\n"
            "    const Pstream::commsTypes commsType\n"
            ") const"
        )   << "Attempting implicit update for field "
            << this->dimensionedInternalField().name()
            << " on patch " << this->patch().name()
            << " with coupledFringe switched off"
            << abort(FatalError);
    }
}


template<class Type>
void oversetFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix& m,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{}


template<class Type>
void oversetFvPatchField<Type>::patchFlux
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& pFlux,
    const fvMatrix<Type>& matrix
) const
{
    // Correct fluxes on fringe faces
    const unallocLabelList& own = matrix.lduAddr().lowerAddr();
    const unallocLabelList& nei = matrix.lduAddr().upperAddr();

    const labelList& fringeFaces = oversetPatch_.overset().fringeFaces();

    const Field<Type>& psi = matrix.psi().internalField();
    Field<Type>& fluxIn = pFlux.internalField();

    // Note that fringe corects coefficients on internal faces
    forAll (fringeFaces, fringeI)
    {
        // Get addressing
        const label& faceI = fringeFaces[fringeI];

        // Multiplication is only done for internal fringe faces
        // HJ, 1/May/2015
        if (faceI < own.size())
        {
            // HJ, check signs!!!
            fluxIn[faceI] = fringeUpperCoeffs_[fringeI]*psi[nei[faceI]]
                - fringeLowerCoeffs_[fringeI]*psi[own[faceI]];
        }
    }
}


template<class Type>
void oversetFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.writeKeyword("coupledFringe") << coupledFringe_
        << token::END_STATEMENT << nl;

    os.writeKeyword("setHoleCellValue")
        << setHoleCellValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("holeCellValue")
        << holeCellValue_ << token::END_STATEMENT << nl;

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
