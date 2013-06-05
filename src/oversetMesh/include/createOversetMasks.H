//
// createOversetMask.H
// ~~~~~~~~~~~~

    const oversetMesh& om = oversetMesh::New(mesh);

    const volScalarField& cellOversetMask = om.gamma();
    cellOversetMask.write();

    const volScalarField& cellOversetMaskExt = om.gammaExt();
    cellOversetMaskExt.write();

    const surfaceScalarField& faceOversetMask = om.sGamma();
    faceOversetMask.write();

    // Write regions and overset types
    om.regionID().write();
    om.oversetTypes().write();
