//
// dynamicOversetMask.H
// ~~~~~~~~~~~~

    const oversetMesh& om = oversetMesh::New(mesh);

    const surfaceScalarField& faceOversetMask = om.sGamma();
