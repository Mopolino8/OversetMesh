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

Class
    oversetInterpolation

Description
    Virtual base class for overset interpolation selection.
    Handles cell centres of global (remote) acceptors indexed by remote donors
    for parallel runs using higher order interpolation schemes (inverse
    distance, least squares...). Stored this way for easier use in derived
    classes.

    Derived classes handle cached weights which are stored in two parts:
    1. localWeights - where donors are on the same processor as acceptor
    2. remoteWeights - where donors are on remote processor
       Note: only master processor holds remoteWeights.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

SourceFiles
    oversetInterpolation.C
    newOversetInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef oversetInterpolation_H
#define oversetInterpolation_H

#include "dictionary.H"
#include "typeInfo.H"
#include "labelList.H"
#include "scalarField.H"
#include "vectorField.H"
#include "scalarFieldField.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class oversetMesh;

/*---------------------------------------------------------------------------*\
                    Class oversetInterpolation Declaration
\*---------------------------------------------------------------------------*/

class oversetInterpolation
{
public:

    // Public data types

        typedef FieldField<Field, scalar> ScalarFieldField;
        typedef List<ScalarFieldField> ListScalarFieldField;

        typedef List<vectorField> ListVectorField;


private:

    // Private data

        //- Overset mesh reference
        const oversetMesh& overset_;

        //- Run time selectable donor extension level:
        // 0 = one level of neighbours are used (immediate neighbours)
        // 1 = two levels of neighbours are used (neighbours of neighbours)
        // 2 = ...
        // Extension level is used by oversetMesh class when marking local and
        // remote donors for given acceptor. Default value is 0.
        const label extensionLevel_;


        // Demand-driven data

            //- Cell centres of global acceptors indexed by local donor indices
            //  for easier calculation of weights in parallel runs for higher
            //  order interpolation schemes
            mutable ListVectorField* remoteAccCCPtr_;


    // Private Member Functions

        //- Disallow copy construct
        oversetInterpolation(const oversetInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const oversetInterpolation&);


        //- Calculate remote acceptor cell centres
        void calcRemoteAccCC() const;

        //- Clear remote acceptor cell centres
        void clearRemoteAccCC() const;


public:

    //- Runtime type information
    TypeName("oversetInterpolation");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            oversetInterpolation,
            dictionary,
            (
                const oversetMesh& overset,
                const dictionary& dict
            ),
            (overset, dict)
        );


    // Selectors

        //- Return a reference to the selected thermal model
        static autoPtr<oversetInterpolation> New
        (
            const oversetMesh& overset,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        oversetInterpolation
        (
            const oversetMesh& overset,
            const dictionary& dict
        );


    // Destructor

        virtual ~oversetInterpolation();


    // Member Functions

        // Access

            //- Return reference to overset mesh
            const oversetMesh& overset() const
            {
                return overset_;
            }

            //- Return extension level
            label extensionLevel() const
            {
                return extensionLevel_;
            }

            //- Return local donor weights
            virtual const ScalarFieldField& localWeights() const = 0;

            //- Return remote donor weights
            virtual const ListScalarFieldField& remoteWeights() const = 0;

            //- Return remote acceptor cell centres
            virtual const ListVectorField& remoteAccCC() const;


        //- Update the interpolation
        virtual void update() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
