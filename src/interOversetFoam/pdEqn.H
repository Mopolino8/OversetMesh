{
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

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rAUf, pd) == fvc::div(phi)
        );

        // Adjust non-orthogonal fringe fluxes if necessary
        om.correctNonOrthoFluxes(pdEqn, U);

        pdEqn.setReference(pRefCell, pRefValue);

        pdEqn.solve
        (
            mesh.solutionDict().solver(pd.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pdEqn.flux();
        }
    }

    U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();
}
