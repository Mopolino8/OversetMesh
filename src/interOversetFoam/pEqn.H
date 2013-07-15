{
    if (nOuterCorr != 1)
    {
        pd.storePrevIter();
    }

    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    U = rUA*UEqn.H();
    // Overset update
    U.correctBoundaryConditions();

    surfaceScalarField phiU
    (
        "phiU",
        faceOversetMask*(fvc::interpolate(U) & mesh.Sf())
    );

    // Adjust overset fluxes
    oversetAdjustPhi(phiU, U);
    adjustPhi(phiU, U, pd);

    phi = phiU
      + faceOversetMask*
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rUAf*mesh.magSf();


    for(int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAf, pd) == fvc::div(phi)
        );

        pdEqn.setReference(pdRefCell, pdRefValue);

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

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();
}
