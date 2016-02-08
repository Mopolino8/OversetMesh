{
    if (nOuterCorr != 1)
    {
        pd.storePrevIter();
    }

    rAU = 1.0/UEqn.A();
    rAU.correctBoundaryConditions(); // Overset update
    surfaceScalarField rAUf = fvc::interpolate(rAU);

    U = rAU*UEqn.H();
    U.correctBoundaryConditions(); // Overset update

    surfaceScalarField phiU
    (
        "phiU",
        faceOversetMask*(fvc::interpolate(U) & mesh.Sf())
    );

    // Adjust overset fluxes
    oversetAdjustPhi(phi, U); // Fringe flux adjustment
    globalOversetAdjustPhi(phi, U, p); // Global flux adjustment

    phi = phiU
      + faceOversetMask*
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf();


    for(int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rAUf, pd) == fvc::div(phi)
        );

        // Adjust non-orthogonal fringe fluxes if necessary
        om.correctNonOrthoFluxes(pdEqn, U);

        pdEqn.setReference(pRefCell, pRefValue);

        if (corr == nCorr - 1 && nonOrth == nNonOrthCorr)
        {
            pdEqn.solve(mesh.solutionDict().solver(pd.name() + "Final"));
        }
        else
        {
            pdEqn.solve(mesh.solutionDict().solver(pd.name()));
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pdEqn.flux();
        }
    }

    // Explicitly relax pressure except for last corrector
    if (oCorr != nOuterCorr - 1)
    {
        pd.relax();
    }

    U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();
}
