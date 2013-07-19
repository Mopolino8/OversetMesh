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

Class
    Foam::regionSolidBodyMotionFvMesh

Description
    Solid-body motion of the mesh specified by a run-time selectable
    motion function.

SourceFiles
    regionSolidBodyMotionFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef regionSolidBodyMotionFvMesh_H
#define regionSolidBodyMotionFvMesh_H

#include "dynamicFvMesh.H"
#include "dictionary.H"
#include "pointIOField.H"
#include "regionSplit.H"
#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class regionSolidBodyMotionFvMesh Declaration
\*---------------------------------------------------------------------------*/

class regionSolidBodyMotionFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Dictionary of motion control parameters
        dictionary dynamicMeshCoeffs_;

        //- Region split
        regionSplit rs_;

        //- Motion control function
        PtrList<solidBodyMotionFunction> regionSBMFs_;

        //- Region motion masks
        mutable FieldField<Field, scalar>* regionMotionMasksPtr_;

        //- Reference points which are transformed
        pointIOField undisplacedPoints_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionSolidBodyMotionFvMesh(const regionSolidBodyMotionFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const regionSolidBodyMotionFvMesh&);


        //- Calculate motion masks
        void calcMotionMasks() const;


public:

    //- Runtime type information
    TypeName("regionSolidBodyMotionFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit regionSolidBodyMotionFvMesh(const IOobject& io);


    // Destructor

        virtual ~regionSolidBodyMotionFvMesh();


    // Member Functions

        //- Return motion masks for regions
        const scalarField& motionMask(const label regionI) const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
