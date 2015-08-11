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
    Foam::oversetMesh

Description
   Overset mesh class

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    oversetMesh.C
    oversetMeshAddressing.C
    oversetMeshInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef oversetMesh_H
#define oversetMesh_H

#include "MeshObject.H"
#include "foamTime.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "oversetRegion.H"
#include "oversetInterpolation.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class oversetMesh Declaration
\*---------------------------------------------------------------------------*/

class oversetMesh
:
    public MeshObject<fvMesh, oversetMesh>
{
public:

    // Public enumerations

    enum cellType
    {
        HOLE = -1,
        ACTIVE = 0,
        DONOR = 1,
        ACCEPTOR = 2
    };


private:

    // Private data

        //- Coefficient dictionary
        IOdictionary dict_;

        //- List of overset regions
        PtrList<oversetRegion> regions_;

        //- Hole patch names
        wordList holePatchNames_;


    // Demand-driven data

        // Cell clasification

            //- Acceptor cell labels
            mutable labelList* acceptorCellsPtr_;

            //- Donor cell labels
            mutable labelList* donorCellsPtr_;

            //- Hole cell labels
            mutable labelList* holeCellsPtr_;


        // Domain markup fields

            //- Return overset type indicator field
            mutable volScalarField* oversetTypesPtr_;

            //- Region ID: region index for each cell
            mutable labelList* regionIDPtr_;


        // Overset discretisation support

            //- Live cells indicator
            mutable volScalarField* gammaPtr_;

            //- Non-hole cells indicator, marking live and acceptor cells
            mutable volScalarField* gammaExtPtr_;

            //- Fluid faces indicator, marking faces between live cells
            mutable surfaceScalarField* sGammaPtr_;

            //- List of faces for which one neighbour is a live cell
            //  and another neighbour is an acceptor cell
            mutable labelList* fringeFacesPtr_;

            //- List of fringe cell index for each fringe face
            mutable labelList* fringeFaceCellsPtr_;

            //- List of fringe face flip:
            //  false if fringe face points into fringe cell
            //  (out of the live cell)
            //  true if fringe face points into a live cell
            mutable boolList* fringeFaceFlipsPtr_;

            //- List of faces between hole and acceptor cells
            mutable labelList* holeFacesPtr_;

            //- List of fringe cell index for each hole face
            mutable labelList* holeFaceCellsPtr_;

            //- List of hole face flip:
            //  false if hole face points into fringe cell
            //  (out of the hole cell)
            //  true if hole face points into a live cell
            mutable boolList* holeFaceFlipsPtr_;

            //- List of faces in the region bounded by hole faces
            mutable labelList* holeInternalFacesPtr_;

            //- List of faces between two acceptor cells
            mutable labelList* acceptorInternalFacesPtr_;


        // Interpolation addressing

            //- Local donors cells: cells where the donor is on the
            //  same processor as the acceptor
            mutable labelList* localDonorsPtr_;

            //- Local donors addressing: acceptor index for each
           //  local donor
            mutable labelList* localDonorAddrPtr_;

            //- Remote donors cells: a local donor cell for which the
            //  acceptor is on a remote processor
            mutable labelList* remoteDonorsPtr_;

            //- Remote acceptor addressing: local acceptor cell where the
            //  donor is on a remote processor
            mutable labelList* remoteAcceptorAddrPtr_;


            // Master processor data re-mapping addressing
            // This data exists only on master processor

                //- Parallel data information processor addressing
                //  For each global acceptor, record donor processor index
                mutable labelListList* globalAcceptFromProcPtr_;

                //- Parallel data information cell addressing
                //  For each global acceptor, record donor cell index
                mutable labelListList* globalAcceptFromCellPtr_;


            //- Overset interpolation
            mutable oversetInterpolation* interpolationPtr_;


    // Private member functions

        //- Calculate cell classification
        void calcCellClassification() const;

        //- Calculate domain type fields
        void calcDomainMarkup() const;

        //- Calculate gamma fields
        void calcGamma() const;

        //- Calculate fringe faces
        void calcFringeFaces() const;

        //- Calculate hole faces
        void calcHoleFaces() const;

        //- Calculate parallel fringe interpolation addressing
        void calcParallelAddressing() const;


        // Storage management

            //- Clear all demand-driven data
            void clearOut() const;


public:

    // Declare name of the class and its debug switch
    TypeName("oversetMesh");


    // Constructors

        //- Construct given an fvMesh
        explicit oversetMesh(const fvMesh&);

        //- Selection
        static const oversetMesh& New(const fvMesh& mesh)
        {
            return MeshObject<fvMesh, oversetMesh>::New(mesh);
        }


    //- Destructor
    virtual ~oversetMesh();


    // Member functions

        //- Return list of overset regions
        const PtrList<oversetRegion>& regions() const
        {
            return regions_;
        }

        //- Return hole patch names
        const wordList& holePatchNames() const
        {
            return holePatchNames_;
        }


        // Cell clasification

            //- Return acceptor cells
            const labelList& acceptorCells() const;

            //- Return donor cells
            const labelList& donorCells() const;

            //- Return hole cells
            const labelList& holeCells() const;

            //- Return cell type indicator field
            // -1 = hole
            //  0 = active cell
            //  1 = donor
            //  2 = acceptor
            const volScalarField& oversetTypes() const;

            //- Return region indicator
            const labelList& regionID() const;


        // Overset discretisation support

            //- Return fluid cells indicator, marking only live cells
            const volScalarField& gamma() const;

            //- Return fluid cells indicator, marking live and acceptor cells
            const volScalarField& gammaExt() const;

            //- Return fluid faces indicator, marking faces between live cells
            const surfaceScalarField& sGamma() const;

            //- Return list of faces for which one neighbour is an live cell
            //  and another neighbour is an acceptor cell
            const labelList& fringeFaces() const;

            //- Return fringe cell index for each fringe face
            const labelList& fringeFaceCells() const;

            //- Return list of fringe face flip:
            //  false if fringe face points into fringe cell
            //  (out of the live cell)
            //  true if fringe face points into a live cell
            const boolList& fringeFaceFlips() const;

            //- Return list of faces between hole and acceptor cells
            const labelList& holeFaces() const;

            //- Return fringe cell index for each hole face
            const labelList& holeFaceCells() const;

            //- Return list of fringe face flip:
            //  false if fringe face points into fringe cell
            //  (out of the hole cell)
            //  true if fringe face points into a hole cell
            const boolList& holeFaceFlips() const;

            //- Return list of internal faces in the region bounded by
            //  hole faces
            const labelList& holeInternalFaces() const;

            //- Return list of faces between two acceptor cells
            const labelList& acceptorInternalFaces() const;


        // Fringe and interpolation

           //- Return local donor cells: cells where donor and
           //  acceptor are on the same processor
           const labelList& localDonors() const;

           //- Return local donor addressing: acceptor index for each
           //  local donor
           const labelList& localDonorAddr() const;

           //- Return remote donor cells: cells where the donor is local
           //  and acceptor is on a remote processor
           const labelList& remoteDonors() const;

           //- Return remote acceptor addressing: local acceptor index
           //  for which donor is on a remote processor
           const labelList& remoteAcceptorAddr() const;


            // Master processor data re-mapping addressing
            // This data exists only on master processor

                //- Parallel data information processor addressing
                //  For each global acceptor, record donor processor index
               const labelListList& globalAcceptFromProc() const;

                //- Parallel data information cell addressing
                //  For each global acceptor, record donor cell index
                const labelListList& globalAcceptFromCell() const;
    

            //- Return interpolation
            const oversetInterpolation& interpolation() const;

            //- Interpolate to acceptors
            template<class Type>
            void interpolate
            (
                Field<Type>& accF,
                const Field<Type>& cellField
            ) const;


            //- Interpolate to acceptors
            template<class Type>
            tmp<Field<Type> > interpolate
            (
                const Field<Type>& cellField
            ) const;



        // Mesh object update functions

            //- Update after mesh motion
            virtual bool movePoints() const;

            //- Update after topological change
            virtual bool updateMesh(const mapPolyMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "oversetMeshInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
