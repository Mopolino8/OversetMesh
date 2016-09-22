{
    p.boundaryField().updateCoeffs();

    rAU = 1.0/UEqn.A();
    rAU.correctBoundaryConditions(); // Overset update
    surfaceScalarField rAUf = fvc::interpolate(rAU);

    U = rAU*UEqn.H();
    U.correctBoundaryConditions(); // Overset update

    surfaceScalarField phiU
    (
        "phiU",
        fvc::interpolate(U) & mesh.Sf()
    );

    phi = phiU;

    // Adjust overset fluxes
    oversetAdjustPhi(phi, U); // Fringe flux adjustment
    globalOversetAdjustPhi(phi, U, p); // Global flux adjustment

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUf, p) == fvc::div(phi)
        );

        // Adjust non-orthogonal fringe fluxes if necessary
        om.correctNonOrthoFluxes(pEqn, U);

        pEqn.setReference(pRefCell, pRefValue);

        if
        (
            (
                pimple.corrPISO() == 1
             || pimple.corrPISO() == pimple.nCorrPISO()
            )
         && pimple.finalNonOrthogonalIter()
        )
        {
            pEqn.solve(mesh.solutionDict().solver("pFinal"));
        }
        else
        {
            pEqn.solve();
        }

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }

#   include "oversetContinuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    if (!pimple.finalIter())
    {
        p.relax();
    }

#   include "movingMeshContinuityErrs.H"

    U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
}
