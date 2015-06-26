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
    om.oversetTypes().write();

    // Create region ID field
    volScalarField regionIndex
    (
        IOobject
        (
            "regionIndex",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    );

    {
        scalarField& regionIndexIn = regionIndex.internalField();

        const labelList& regionID = om.regionID();

        forAll (regionIndexIn, cellI)
        {
            regionIndexIn[cellI] = regionID[cellI];
        }


        // Update boundary values
        forAll (regionIndex.boundaryField(), patchI)
        {
            if (!regionIndex.boundaryField()[patchI].coupled())
            {
                regionIndex.boundaryField()[patchI] =
                    regionIndex.boundaryField()[patchI].patchInternalField();
            }
        }
    }

    regionIndex.write();
