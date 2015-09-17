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
    overlapFringe

Description
    Overlap fringe algorithm

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    overlapFringe.C

\*---------------------------------------------------------------------------*/

#ifndef overlapFringe_H
#define overlapFringe_H

#include "oversetFringe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class overlapFringe Declaration
\*---------------------------------------------------------------------------*/

class overlapFringe
:
    public oversetFringe
{
    // Private data

        //- Fringe hole cells
        mutable labelList* fringeHolesPtr_;

        //- Acceptor cells
        mutable labelList* acceptorsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        overlapFringe(const overlapFringe&);

        //- Disallow default bitwise assignment
        void operator=(const overlapFringe&);


        // Calculate hole-acceptor addressing
        void calcAddressing() const;

        // Clear hole-acceptor addressing
        void clearAddressing() const;


public:

    //- Runtime type information
    TypeName("overlap");


    // Constructors

        //- Construct from dictionary
        overlapFringe
        (
            const fvMesh& mesh,
            const oversetRegion& region,
            const dictionary& dict
        );


    // Destructor

        virtual ~overlapFringe();


    // Member Functions

        //- Return list of hole cells
        //  Fringe hole cells are collected in addition to geometric hole
        //  cells, which fall outside of all donor regions
        virtual const labelList& fringeHoles() const;

        //- Return list of acceptor cells
        virtual const labelList& acceptors() const;

        //- Update the fringe
        virtual void update() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
