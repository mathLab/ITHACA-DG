
    Foam::Info
        << "Create mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

   _mesh = autoPtr<dgMesh>
    (
        new dgMesh
        (
            Foam::IOobject
            (
                Foam::dgMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        )
    );
    Foam::dgMesh& mesh = _mesh();
