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

Application
    scalarTransportOversetFoam

Description
    Solves a transport equation for a passive scalar with support for
    overset meshes.

    Experimental: no overset mesh changes required

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    simpleControl simple(mesh);

#   include "createFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

#   include "CourantNo.H"

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        fvScalarMatrix TEqn
        (
            fvm::ddt(T)
          + fvm::div(phi, T)
          - fvm::laplacian(DT, T)
        );

        TEqn.solve();

//         volScalarField Tresidual
//         (
//             "Tresidual",
//             T
//         );
//         Tresidual.writeOpt() = IOobject::AUTO_WRITE;

//         Tresidual.internalField() = TEqn.residual();
//         Tresidual.boundaryField() == 0;

//         Info<< "residual " << gSumMag(Tresidual.internalField()) << endl;

        fvVectorMatrix VEqn
        (
            fvm::ddt(V)
          + fvm::div(phi, V)
          - fvm::laplacian(DV, V)
        );

        VEqn.solve();

//         return 0;

//         VEqn.solve();

//         volVectorField Vresidual
//         (
//             "Vresidual",
//             V
//         );
//         Vresidual.writeOpt() = IOobject::AUTO_WRITE;

//         Vresidual.internalField() = VEqn.residual();
//         Vresidual.boundaryField() == vector::zero;

//         Info<< "residual " << gSumMag(Vresidual.internalField()) << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
