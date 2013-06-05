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
    injectionInterpolation

Description
    Injection overset interpolation algorithm

    Donor cells are collected directly

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    injectionInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef injectionInterpolation_H
#define injectionInterpolation_H

#include "oversetInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class injectionInterpolation Declaration
\*---------------------------------------------------------------------------*/

class injectionInterpolation
:
    public oversetInterpolation
{
    // Private data

        // Demand-driven data

            //- Addressing
            mutable labelListList* addressingPtr_;

            //- Weights
            mutable scalarListList* weightsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        injectionInterpolation(const injectionInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const injectionInterpolation&);


        //- Calculate addressing
        void calcAddressing() const;

        //- Clear addressing
        void clearAddressing() const;


public:

    //- Runtime type information
    TypeName("injection");


    // Constructors

        //- Construct from dictionary
        injectionInterpolation
        (
            const oversetMesh& overset,
            const dictionary& dict
        );


    // Destructor

        virtual ~injectionInterpolation();


    // Member Functions

            //- Return interpolation donor list
            virtual const labelList& donors() const;

            //- Return acceptor addressing
            virtual const labelListList& addressing() const;

            //- Return acceptor weights
            virtual const scalarListList& weights() const;


        //- Update the interpolation
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
