/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    potentialDyMOversetFoam

Description
    Transient solver for potential flow with dynamic overset mesh.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "oversetMesh.H"
#include "oversetFvPatchFields.H"
#include "oversetAdjustPhi.H"
#include "globalOversetAdjustPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createFields.H"
#   include "initTotalVolume.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
#       include "readPISOControls.H"
#       include "checkTotalVolume.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        bool meshChanged = mesh.update();
        reduce(meshChanged, orOp<bool>());

#       include "createOversetMasks.H"

        // Update moving wall velocity boundary condition and calculate the flux
        U.correctBoundaryConditions();
        phi == (linearInterpolate(U) & mesh.Sf());

        // Resetting pressure field
        p.internalField() = 0;

#       include "volContinuity.H"
#       include "meshCourantNo.H"

        // Solve potential flow equations

        oversetAdjustPhi(phi, U); // Fringe flux adjustment
        globalOversetAdjustPhi(phi, U, p); // Global flux adjustment

        for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
        {
            p.storePrevIter();

            Info<< "Initial flux contour continuity error = "
                << mag(sum(phi.boundaryField()))
                << endl;

            fvScalarMatrix pEqn
            (
                fvm::laplacian
                (
                    dimensionedScalar
                    (
                        "1",
                        dimTime/p.dimensions()*dimensionSet(0, 2, -2, 0, 0),
                        1
                    ),
                    p
                )
             ==
                fvc::div(phi)
            );

            // Adjust non-orthogonal fringe fluxes if necessary
            om.correctNonOrthoFluxes(pEqn, U);

            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();

            if (nonOrth != nNonOrthCorr)
            {
                p.relax();
            }
            else
            {
                phi -= pEqn.flux();
#               include "oversetContinuityErrs.H"
            }
        }

        // Update div phi field for visualisation purposes
        oversetDivPhi = cellOversetMask*fvc::div(phi);

        // Note: do not reconstruct the velocity field for next time step
//        U = fvc::reconstruct(phi);
//        U.correctBoundaryConditions();

        Info<< "Interpolated U error = "
            << (
                   sqrt
                   (
                       sum
                       (
                           sqr
                           (
                               faceOversetMask*
                               (
                                   (fvc::interpolate(U) & mesh.Sf())
                                 - phi
                               )
                           )
                       )
                   )/sum(mesh.magSf())
               ).value()
            << endl;

        // Calculate velocity magnitude
        {
            volScalarField magU = mag(U);

            Info<< "mag(U): max: " << gMax(magU.internalField())
                << " min: " << gMin(magU.internalField()) << endl;
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
