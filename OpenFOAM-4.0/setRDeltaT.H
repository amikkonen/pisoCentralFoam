{
    volScalarField& rDeltaT = trDeltaT.ref();

    scalar rDeltaTSmoothingCoeff
    (
        runTime.controlDict().lookupOrDefault<scalar>
        (
            "rDeltaTSmoothingCoeff",
            0.02
        )
    );

    scalar rDeltaTDampingCoeff
    (
        runTime.controlDict().lookupOrDefault<scalar>
        (
            "rDeltaTDampingCoeff",
            1.0
        )
    );

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);



    // Seems like there is a bug here somewhere.
    // if CourantType == "cellCourant" gives the following error: 
    //    --> FOAM FATAL ERROR: 
    //Arguments of max have different dimensions
    //     dimensions : [0 0 -1 0 0 0 0] and [1 -3 -1 0 0 0 0]


    //    From function Foam::dimensionSet Foam::max(const Foam::dimensionSet&, const Foam::dimensionSet&)
    //    in file dimensionSet/dimensionSet.C at line 231.

    //FOAM aborting


    if ( CourantType == "cellCourant" )
    {
        // Set the reciprocal time-step from the local Courant number

        rDeltaT.ref() = max
        (
	  1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
	  
	      // I'm guessing the error is here.
          (
	      0.5 * fvc::surfaceSum(mag(phi))()()
	      /mesh.V()
	      )/ maxCo
        );
    }
    else if ( CourantType == "faceCourant" )
    {
        if (mesh.nInternalFaces())
        {
            
            surfaceScalarField phiOverRho =
            mag
            (
                phi / (p_pos * psi_pos * a_pos + p_neg * psi_neg * a_neg)
            );
            
            surfaceScalarField SfUfbyDelta =
                mesh.surfaceInterpolation::deltaCoeffs()*phiOverRho*2.0;
            
            forAll(mesh.V(), iCell)
            {
                const labelList& cellFaces = mesh.cells()[iCell];
                
                scalar maxCellRDeltaT = -1.0;
                scalar faceRDeltaT = maxCellRDeltaT;
                label faceId = -1;
                forAll(cellFaces, jFace)
                {
                    faceId = cellFaces[jFace];
                    if (mesh.isInternalFace(faceId))
                    {
                        faceRDeltaT = SfUfbyDelta.internalField()[faceId] / mesh.magSf().internalField()[faceId] / maxCo;
                    }
                    if (faceRDeltaT > maxCellRDeltaT)
                    {
                        maxCellRDeltaT = faceRDeltaT;
                    }
                }
                rDeltaT.ref()[iCell] = maxCellRDeltaT;
            }
            forAll(SfUfbyDelta.boundaryField(), kPatch)
            {
                if(polyPatch::constraintType(mesh.boundaryMesh().types()[kPatch]))
                {
                }
                else
                {
                    const fvsPatchField<scalar>& pSfUfbyDelta = SfUfbyDelta.boundaryField()[kPatch];
                    const labelList& patchCells = mesh.boundary()[kPatch].faceCells();
                    scalar faceRDeltaT = -1;
                    forAll(pSfUfbyDelta, qFace)
                    {
                        faceRDeltaT = pSfUfbyDelta[qFace] / mesh.magSf().boundaryField()[kPatch][qFace] / maxCo;
                        if (faceRDeltaT > rDeltaT.internalField()[patchCells[qFace]])
                        {
                            rDeltaT.ref()[patchCells[qFace]] = faceRDeltaT;
                        }
                    }
                }
            }
        }
        
    }
    else
    {
        FatalErrorIn
        (
            "setRDeltaT.H:"
        )
        << "Wrong type of Courant criterion: " << CourantType
        << endl << " must be one of:" 
        << endl << "1) cellCourant"
        << endl << "2) faceCourant"
        << endl << abort(FatalError);
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/(rDeltaT.primitiveField() + 1.0e-12))
        << ", " << gMax(1/(rDeltaT.primitiveField() + 1.0e-12)) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT =
            rDeltaT0
           *max(rDeltaT/rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

        Info<< "Damped flow time scale min/max = "
            << gMin(1/rDeltaT.primitiveField())
            << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
    }
}
