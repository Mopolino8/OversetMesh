    volScalarField AU = UEqn().A();

    U = UEqn().H()/AU;
    // Overset update
    U.correctBoundaryConditions();

    UEqn.clear();
    phi = faceOversetMask*(fvc::interpolate(U) & mesh.Sf());

    // Adjust overset fluxes
    oversetAdjustPhi(phi, U);
    adjustPhi(phi, U, p);

    // Non-orthogonal pressure corrector loop
    for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(1.0/AU, p) == fvc::div(phi)
        );

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
    U -= fvc::grad(p)/AU;
    U.correctBoundaryConditions();

