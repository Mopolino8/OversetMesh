//
// createOversetMask.H
// ~~~~~~~~~~~~
    const oversetMesh& om = oversetMesh::New(mesh);

    const volScalarField& cellOversetMask = om.gamma();

    const volScalarField& cellOversetMaskExt = om.gammaExt();

    const surfaceScalarField& faceOversetMask = om.sGamma();

    if (runTime.outputTime())
    {
        om.gamma().write();
        om.gammaExt().write();
        om.sGamma().write();
        om.oversetTypes().write();
    }
