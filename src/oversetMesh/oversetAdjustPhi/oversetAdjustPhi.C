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

#include "oversetAdjustPhi.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "oversetMesh.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::oversetAdjustPhi
(
    surfaceScalarField& phi,
    volVectorField& U
)
{
    const fvMesh& mesh = phi.mesh();

    // If the mesh is moving, adjustment needs to be calculated on
    // relative fluxes.  HJ, 13/Feb/2009
    if (mesh.moving())
    {
        fvc::makeRelative(phi, U);
    }

    Info<< "Overset adjust phi" << endl;

    // Get overset mesh
    const oversetMesh& om = oversetMesh::New(mesh);

    // Get addressing to fringe faces
    const labelList& fringeFaces = om.fringeFaces();
    const boolList& fringeFaceFlips = om.fringeFaceFlips();

    // Get internal owner-neighbour addressing
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Get region split to identify separate mesh components
    const regionSplit& rs = om.split();

    // Sum up incoming and outgoing flux
    scalarField regionFringeIn(rs.nRegions(), 0);
    scalarField regionFringeOut(rs.nRegions(), 0);

    scalarField regionFringeMagBalance(rs.nRegions(), 0);

    // Adjust fluxes per region
    scalarField& phiIn = phi.internalField();

    forAll (fringeFaces, ffI)
    {
        const label curFace = fringeFaces[ffI];
        const bool curFlip = fringeFaceFlips[ffI];

        // Get region index
        const label curRegion = rs[owner[curFace]];

        // Debug check
        if (rs[owner[curFace]] != rs[neighbour[curFace]])
        {
            FatalErrorIn
            (
                "void Foam::oversetAdjustPhi\n"
                "(\n"
                "    surfaceScalarField& phi,\n"
                "    volVectorField& U\n"
                ")"
            )   << "Region index different for owner and neighbour of face "
                << curFace
                << abort(FatalError);
        }

        if (!curFlip)
        {
            // Correct orientation
            if (phiIn[curFace] < 0)
            {
                // Incoming flux
                regionFringeIn[curRegion] -= phiIn[curFace];
            }
            else
            {
                // Outgoing flux
                regionFringeOut[curRegion] += phiIn[curFace];
            }
        }
        else
        {
            // Flipped face; opposite sign of flux
            if (phiIn[curFace] > 0)
            {
                // Incoming flux
                regionFringeIn[curRegion] += phiIn[curFace];
            }
            else
            {
                // Outgoing flux
                regionFringeOut[curRegion] -= phiIn[curFace];
            }
        }
    }

    Info<< "Region fringe balance: in = " << regionFringeIn
        << " out = " << regionFringeOut
        << " balance = " << regionFringeOut - regionFringeIn
        << endl;

    // Calculate region flux correction
    scalarField regionFluxScale = regionFringeIn/(regionFringeOut + SMALL);

    Info<< "regionFluxScale: " << regionFluxScale << endl;

    // Go through all fringe faces on each region and balance the fluxes
    forAll (fringeFaces, ffI)
    {
        const label curFace = fringeFaces[ffI];
        const bool curFlip = fringeFaceFlips[ffI];

        // Get region index
        const label curRegion = rs[owner[curFace]];

        // Scale outgoing flux to match the incoming one
        if (!curFlip)
        {
            // Correct orientation
            if (phiIn[curFace] > 0)
            {
                // Scale outgoing flux
                phiIn[curFace] *= regionFluxScale[curRegion];
            }
        }
        else
        {
            // Flipped face; opposite sign of flux
            if (phiIn[curFace] < 0)
            {
                // Scale outgoing flux
                phiIn[curFace] *= regionFluxScale[curRegion];
            }
        }
    }

    // If the mesh is moving, adjustment needs to be calculated on
    // relative fluxes.  Now reverting to absolute fluxes.  HJ, 13/Feb/2009
    if (mesh.moving())
    {
        fvc::makeAbsolute(phi, U);
    }
}


// ************************************************************************* //
