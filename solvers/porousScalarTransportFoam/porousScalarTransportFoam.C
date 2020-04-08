/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    porouScalarTransportFoam

Description
    Solves the transport equation for a passive scalar
    in porous media with dispersion coefficient model

Author
    Pierre Horgue

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "multiscalarMixture.H"
#include "reactionModel.H"
#include "sourceEventFile.H"
#include "patchEventFile.H"
#include "outputEventFile.H"
#include "eventFlux.H"
#include "EulerD3dt3Scheme.H"
#include "EulerD2dt2Scheme.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readReactions.H"
    #include "readTimeControls.H"
    #include "readEvent.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    bool redoTimeStep = false;
    bool determineCellSizes = true;
    int breakLoop = 0;

    const auto& ref = composition.Y(0);     // (see if this can be put in somewhere else)
    volScalarField cellSizes(       // h in TauD equation (see if this can be put somewhere else)
        IOobject(
            word("cellSizes"),
            ref.mesh().time().timeName(),
            ref.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        ref.mesh(),
        dimensionedScalar("",dimless,2.0)
    );
  
    while (runTime.run())
    {

        scalar maxDCVariation = runTime.controlDict().lookupOrDefault<scalar>("variationMax",0.0125);
        scalar kValueTauD = runTime.controlDict().lookupOrDefault<scalar>("kValueTauD",0.005);

        #include "CourantNo.H"
        forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateIndex(runTime.timeOutputValue());
        if (outputEventIsPresent) outputEvent.updateIndex(runTime.timeOutputValue());
        forAll(sourceEventList,sourceEventi) sourceEventList[sourceEventi]->updateIndex(runTime.timeOutputValue());

        redoTimeStep = false;

        //#include "setDeltaT.H"      // Not used anymore, new deltaT adjust is done in solveReactiveTransport.H
        runTime++;

        do{
            if(redoTimeStep){
                Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
                
            }

            Info << "Time = " << runTime.timeName() << nl << endl;
            #include "solveReactiveTransport.H"

        } while (redoTimeStep);

        #include "CmassBalance.H"

        #include "eventWrite.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
