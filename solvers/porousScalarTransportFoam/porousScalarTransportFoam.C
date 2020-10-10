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

#include "loggerLevels.H"
const int LOG_LEVEL = 1;
LoggerLevel<1> InfoL1;
LoggerLevel<2> InfoL2;

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
    scalar maxDCVariation = runTime.controlDict().lookupOrDefault<scalar>("variationMax",0.0125);
    scalar maxDeltaTIncrease = runTime.controlDict().lookupOrDefault<scalar>("maxDeltaTIncrease",0.2);
    scalar deltaTDecrease = runTime.controlDict().lookupOrDefault<scalar>("deltaTDecrease",0.2);
    scalar kValueTauD = runTime.controlDict().lookupOrDefault<scalar>("kValueTauD",0.005);
    
    bool redoTimeStep = false;
    int breakLoop = 0;

    const auto& ref = composition.Y(0);     // (see if this can be put in somewhere else)
    const auto& mesh_temp = ref.mesh();
    volScalarField cellSizes(       // h in TauD equation (see if this can be put somewhere else)
        IOobject(
            word("cellSizes"), mesh_temp.time().timeName(), mesh_temp, IOobject::NO_READ, IOobject::NO_WRITE
            ), ref.mesh(), dimensionedScalar("",dimless,2.0)
        );

    scalar domThickness=0.0;
    const pointField& pp = mesh_temp.points();
    edgeList allEdges = mesh_temp.edges();
    forAll(allEdges, edgeI)
    {
        if (pp[allEdges[edgeI].start()].z() != pp[allEdges[edgeI].end()].z())
        {
            domThickness = fabs(pp[allEdges[edgeI].start()].z() - pp[allEdges[edgeI].end()].z());
            break;
        }
    }
    InfoL1 << "domThickness used: " << domThickness << endl;
    
    surfaceScalarField faceAreas = mesh_temp.magSf();
    surfaceScalarField faceLenghtsSqr = sqr(faceAreas/domThickness);
    surfaceVectorField unitaryNormalToFace = (mesh_temp.Sf()/faceAreas);
    surfaceScalarField tauDField_base(
        IOobject(
            word("vectorTauD"),
            mesh_temp.time().timeName(),
            mesh_temp,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_temp, 
        (kValueTauD)
    );
    tauDField_base *= faceLenghtsSqr;
                    
    scalar refValA = runTime.controlDict().lookupOrDefault<scalar>("refValA",500.0);
    scalar refValB = runTime.controlDict().lookupOrDefault<scalar>("refValB",10.0);
    scalar refValC = runTime.controlDict().lookupOrDefault<scalar>("refValC",10.0);
    double referenceVal[] = {refValA, refValB, refValC};
    
    while (runTime.run())
    {
        #include "CourantNo.H"
        forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateIndex(runTime.timeOutputValue());
        if (outputEventIsPresent) outputEvent.updateIndex(runTime.timeOutputValue());
        forAll(sourceEventList,sourceEventi) sourceEventList[sourceEventi]->updateIndex(runTime.timeOutputValue());

        redoTimeStep = false;

        //#include "setDeltaT.H"      // Not used anymore, new deltaT adjust is done in solveReactiveTransport.H
        runTime++;

        do{
            if(redoTimeStep){
                InfoL2<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
            }

            InfoL2 << "Time = " << runTime.timeName() << nl << endl;
            #include "solveReactiveTransport.H"

        } while (redoTimeStep);

        #include "CmassBalance.H"

        #include "eventWrite.H"

        InfoL2<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        
    }
	
	InfoL1<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	    << nl << endl;
    InfoL1<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
