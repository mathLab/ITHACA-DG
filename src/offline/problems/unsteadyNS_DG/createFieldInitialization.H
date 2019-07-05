_U = autoPtr<dgVectorField>
(
	new dgVectorField
	(
		IOobject
		(
			"U",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
	mesh
	)
);
dgVectorField& U = _U();

Info << "Reading field p\n" << endl;
_p = autoPtr<dgScalarField>
(
	new dgScalarField
	(
		IOobject
		(
			"p",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
       )
     );
dgScalarField& p = _p();


//     Info<< "Reading field p\n" << endl;
//     dgScalarField p
//     (
//         IOobject
//         (
//             "p",
//             runTime.timeName(),
//             mesh,
//             IOobject::MUST_READ,
//             IOobject::AUTO_WRITE
//         ),
//         mesh
//     );

/*
    Info<< "Reading field U\n" << endl;
    dgVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );*/


