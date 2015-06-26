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
    oversetFringe

Description
    Virtual base class for overset fringe selection.

    oversetFringe
    - is owned by a particular oversetRegion
    - defines the hole cells for its region
    - defines acceptor cells for its region

    Donor cells are calculated by oversetRegion, using a list of
    acceptor cells and 
    

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    oversetFringe.C
    newOversetFringe.C

\*---------------------------------------------------------------------------*/

#ifndef oversetFringe_H
#define oversetFringe_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class oversetRegion;

/*---------------------------------------------------------------------------*\
                         Class oversetFringe Declaration
\*---------------------------------------------------------------------------*/

class oversetFringe
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to oversetRegion
        const oversetRegion& region_;


    // Private Member Functions

        //- Disallow copy construct
        oversetFringe(const oversetFringe&);

        //- Disallow default bitwise assignment
        void operator=(const oversetFringe&);


public:

    //- Runtime type information
    TypeName("oversetFringe");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            oversetFringe,
            dictionary,
            (
                const fvMesh& mesh,
                const oversetRegion& region,
                const dictionary& dict
            ),
            (mesh, region, dict)
        );


    // Selectors

        //- Return a reference to the selected thermal model
        static autoPtr<oversetFringe> New
        (
            const fvMesh& mesh,
            const oversetRegion& region,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        oversetFringe
        (
            const fvMesh& mesh,
            const oversetRegion& region,
            const dictionary& dict
        );


    // Destructor

        virtual ~oversetFringe()
        {}


    // Member Functions

        // Access

            //- Return reference to mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return access to region
            const oversetRegion& region() const;


        //- Return list of deactivated (hole) cells
        virtual const labelList& holes() const = 0;

        //- Return list of acceptor cells
        virtual const labelList& acceptors() const = 0;

        //- Update the fringe
        virtual void update() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
