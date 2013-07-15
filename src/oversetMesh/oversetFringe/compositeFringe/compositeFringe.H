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
    compositeFringe

Description
    Composite fringe selection algorithm, consisting of multiple
    fringes on a single region

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    compositeFringe.C

\*---------------------------------------------------------------------------*/

#ifndef compositeFringe_H
#define compositeFringe_H

#include "oversetFringe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class compositeFringe Declaration
\*---------------------------------------------------------------------------*/

class compositeFringe
:
    public oversetFringe
{
    // Private data

        //- List of primitive selection algorithms
        PtrList<oversetFringe> baseFringes_;

        //- Hole cells
        mutable labelList* holesPtr_;

        //- Acceptor cells
        mutable labelList* acceptorsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        compositeFringe(const compositeFringe&);

        //- Disallow default bitwise assignment
        void operator=(const compositeFringe&);


        //- Read base fringes
        void readBaseFringes(const dictionary& dict);

        // Calculate hole and acceptor addressing
        void calcAddressing() const;

        // Clear addressing
        void clearAddressing();


public:

    //- Runtime type information
    TypeName("composite");


    // Constructors

        //- Construct from dictionary
        compositeFringe
        (
            const fvMesh& mesh,
            const oversetRegion& region,
            const dictionary& dict
        );


    // Destructor

        virtual ~compositeFringe();


    // Member Functions

        //- Return list of deactivated (hole) cells
        virtual const labelList& holes() const;

        //- Return list of acceptor cells
        virtual const labelList& acceptors() const;

        //- Update the fringe
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
