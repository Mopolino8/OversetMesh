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
    For each acceptor, the donor data is extracted from the cell where
    the donor centroid falls
    Acceptor value is equal to the donor value, irrespective of the local
    position

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

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

            // Note: for each master donor (acceptor), we have a list of scalars
            // for all donors: master + neighbouring. Master weight is the first
            // [index 0], followed by neighbouring weights.

            //- Weights for local donor/acceptor pairs (all donors and acceptor
            //  are on the same processor)
            mutable ScalarFieldField* localWeightsPtr_;

            // Master processor data (only accessible on master processor)

                //- Weights for remote donor/acceptor pairs (all donors are on a
                //  remote processor). List for all processors is stored on
                //  master processor (trading off memory for CPU performance)
                mutable ListScalarFieldField* remoteWeightsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        injectionInterpolation(const injectionInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const injectionInterpolation&);


        //- Calculate weights
        void calcWeights() const;

        //- Clear weights
        void clearWeights() const;


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

        //- Return local donor weights
        virtual const ScalarFieldField& localWeights() const;

        //- Return remote donor weights
        virtual const ListScalarFieldField& remoteWeights() const;

        //- Update the interpolation
        virtual void update() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
