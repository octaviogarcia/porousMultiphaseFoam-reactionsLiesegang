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

Description
    Timestep for groundwaterTransport2DFoam solver 

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{

    //- C-variation control
    scalar dtForC = VGREAT;
    scalar dtForhwater = VGREAT;
    if (timeScheme == "backward")
    {
        if (hwatermax > SMALL) dtForhwater = Foam::pow(3*truncationError*(hwatermax+VSMALL)/(dhwater3dT3max+VSMALL),1./3.);
        forAll(composition.Y(), speciesi)
        {
            if (Cmax[speciesi] > SMALL)
            {
                dtForC = min
                    (
                        dtForC,
                        Foam::pow(3*truncationError*(Cmax[speciesi]+VSMALL)/(dC3dT3max[speciesi]+VSMALL),1./3.)
                    );
            }
        }
        d3dt3Operator.storeDeltaT00(runTime.deltaT0Value());
    }
    else if (timeScheme == "CrankNicolson")
    {
        if (hwatermax > SMALL) dtForhwater = Foam::pow(12*truncationError*(hwatermax+VSMALL)/(dhwater3dT3max+VSMALL),1./3.);
        forAll(composition.Y(), speciesi)
        {
            if (Cmax[speciesi] > SMALL)
            {
                dtForC = min
                    (
                        dtForC,
                        Foam::pow(12*truncationError*(Cmax[speciesi]+VSMALL)/(dC3dT3max[speciesi]+VSMALL),1./3.)
                    );
            }
        }
        d3dt3Operator.storeDeltaT00(runTime.deltaT0Value());
    }
    else if (timeScheme == "Euler")
    {
        if (hwatermax > SMALL) dtForhwater = Foam::pow(2*truncationError*(hwatermax+VSMALL)/(dhwater2dT2max+VSMALL),1./2.);
        forAll(composition.Y(), speciesi)
        {
            if (Cmax[speciesi] > SMALL)
            {
                dtForC = min
                    (
                        dtForC,
                        Foam::pow(2*truncationError*(Cmax[speciesi]+VSMALL)/(dC2dT2max[speciesi]+VSMALL),1./2.)
                    );
            }
        }
    }

    scalar newDeltaT = min(min(dtForhwater,dtForC), 1.2*runTime.deltaTValue());
    runTime.setDeltaT
        (
            min
            (
                newDeltaT,
                maxDeltaT
            )
        );

    scalar timeOfNextEvent = GREAT;
    if (eventTimeTracking)
    {
        if (infiltrationEventIsPresent) timeOfNextEvent = min(timeOfNextEvent,infiltrationEvent.currentEventEndTime());
        forAll(sourceEventList,sourceEventi) timeOfNextEvent = min(timeOfNextEvent,sourceEventList[sourceEventi]->currentEventEndTime());
    }
    if (outputEventIsPresent) timeOfNextEvent = min(timeOfNextEvent,outputEvent.currentEventEndTime());

    scalar timeToNextEvent = timeOfNextEvent-runTime.timeOutputValue();
    scalar nSteps =  timeToNextEvent/runTime.deltaTValue();
    if ((nSteps < labelMax) && (nSteps != 0))
    {
        const label nStepsToNextEvent = label(max(nSteps, 1) + 0.99);
        runTime.setDeltaTNoAdjust(timeToNextEvent/nStepsToNextEvent);
    }

    //- To handle close event times (inferior to current timestep)
    if (nSteps == 0)
    {
        scalar timeToCloseEvent = GREAT;
        if (eventTimeTracking)
        {
            if (infiltrationEventIsPresent)
            {
                if (infiltrationEvent.currentEventEndTime() != runTime.timeOutputValue())
                {
                    timeToCloseEvent = min(timeToCloseEvent,infiltrationEvent.currentEventEndTime()-runTime.timeOutputValue());
                }
            }
            forAll(sourceEventList,sourceEventi) 
            {
                if (sourceEventList[sourceEventi]->currentEventEndTime() != runTime.timeOutputValue())
                {
                    timeToCloseEvent = min(timeToCloseEvent,sourceEventList[sourceEventi]->currentEventEndTime()-runTime.timeOutputValue());
                }
            }
            forAll(patchEventList,patchEventi)
            {
                if (patchEventList[patchEventi]->currentEventEndTime() != runTime.timeOutputValue())
                {
                    timeToCloseEvent = min(timeToCloseEvent,patchEventList[patchEventi]->currentEventEndTime()-runTime.timeOutputValue());
                }
            }
        }
        if (outputEventIsPresent)
        {
            if (outputEvent.currentEventEndTime() != runTime.timeOutputValue())
            {
                timeToCloseEvent = min(timeToCloseEvent,outputEvent.currentEventEndTime()-runTime.timeOutputValue());
            }
        }
        runTime.setDeltaTNoAdjust(min(runTime.deltaTValue(),timeToCloseEvent));
    }

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
           
}    

// ************************************************************************* //
