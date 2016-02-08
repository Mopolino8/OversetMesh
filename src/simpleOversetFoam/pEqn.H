    rAU = 1.0/UEqn().A();
    rAU.correctBoundaryConditions(); // Overset update

    U = rAU*UEqn().H();
    U.correctBoundaryConditions(); // Overset update

    UEqn.clear();
    phi = faceOversetMask*(fvc::interpolate(U) & mesh.Sf());

    // Adjust overset fluxes
    oversetAdjustPhi(phi, U); // Fringe flux adjustment
    globalOversetAdjustPhi(phi, U, p); // Global flux adjustment

    // Non-orthogonal pressure corrector loop
    for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) == fvc::div(phi)
        );

        // Adjust non-orthogonal fringe fluxes if necessary
        om.correctNonOrthoFluxes(pEqn, U);

        pEqn.setReference(pRefCell, pRefValue);

        // Retain the residual from the first iteration
        if (nonOrth == 0)
        {
            eqnResidual = pEqn.solve().initialResidual();
            maxResidual = max(eqnResidual, maxResidual);
        }
        else
        {
            pEqn.solve();
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pEqn.flux();
        }
    }

#   include "oversetContinuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    U -= rAU*fvc::grad(p);
    U.correctBoundaryConditions();

