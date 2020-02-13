{
    //- update dispersion tensor coefficients and source terms
    forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateValue(runTime);
    forAll(sourceEventList,sourceEventi) sourceEventList[sourceEventi]->updateValue(runTime);
    composition.correct(U, Saturation, eps);

    scalar Cresidual;
    if(reactionModel->needsSubcycling())
    {
        reactionModel->setPrevMaxC(
                            max(
                                (composition.Y( reactionModel->getInfluencedSpecieHS() )).internalField())
                                .value()
                            );

        scalar iterTransport = 0;
        do
        {
            ++iterTransport;

            reactionModel->correct(false);
            #include "CEqn.H"
            
            Info<< "Reactive transport iteration : " << iterTransport << ", Residual = " << Cresidual << endl; 

        } while(iterTransport < 7); 
		//while(Cresidual > reactionTol);
		

        if(reactionMassConservative && !reactionModel->alwaysMassConservative())
        {
            reactionModel->correct(true);
            #include "CEqn.H"

            Info<< "Reactive transport mass conservation correction, Final residual = " << Cresidual << endl;

            if(Cresidual > reactionTol)
            {
                Warning()
                    << "Final residual after mass conservation correction exceeds tolerance"
                    << endl;
            }
        }
    }
    else
    {
        reactionModel->correct(reactionMassConservative);
        #include "CEqn.H"
    }

    reactionModel->postTransport();

    forAll(composition.Y(), speciesi)
    {
        const auto& speciesName = composition.species()[speciesi];
        const auto& C = composition.Y(speciesi);
        const auto& R = composition.R(speciesi);

        dCdTmax[speciesi] = max(mag(fvc::ddt(C))).value();
        if (timeScheme == "Euler")
        {
            volScalarField dC2dT2 = d2dt2Operator.fvcD2dt2(C);
            dC2dT2max[speciesi] = 0;
            forAll(dC2dT2, celli)
            {
                if(mag(dC2dT2[celli]) > dC2dT2max[speciesi])
                {
                    Cmax[speciesi] = C[celli];
                    dC2dT2max[speciesi] = mag(dC2dT2[celli]);
                }
            }
        }
        else
        {
            volScalarField dC3dT3 = d3dt3Operator.d3dt3(C);
            dC3dT3max[speciesi] = 0;
            forAll(dC3dT3, celli)
            {
                if(mag(dC3dT3[celli]) > dC3dT3max[speciesi])
                {
                    Cmax[speciesi] = C[celli];
                    dC3dT3max[speciesi] = mag(dC3dT3[celli]);
                }
            }
        }
        Info<< "Concentration Min(" << speciesName << ") = " << min(C).value() 
            << " Max(" << speciesName << ") = " << max(C.internalField()).value()
            << " mass(" << speciesName << ") = " << fvc::domainIntegrate(R*C*Saturation*eps).value()
            << " dCmax = " << dCdTmax[speciesi]*runTime.deltaTValue()
            << endl;
        
        if(speciesi == reactionModel->getInfluencedSpecieHS()){
            updateDeltaT=true;
            
            //-v
            scalar redCoef = reactionModel->getRedCoef();
            scalar cs = reactionModel->getCsValue();
            scalar cMax = max(C.internalField()).value();
            scalar newDeltaT = runTime.deltaTValue();

            if (cMax >9.0)
            {
                initialIncrease=false;
            }

            if(!initialIncrease && cMax < cs){
                Info<< "-> Computing new deltaT according max concentration..." << endl;
                scalar cMaxPrev = reactionModel->getPrevMaxC();

                scalar difference = cs-cMax;
                Info<< "-> cmaxprev: " << cMaxPrev << endl;
                Info<< "-> cs-cmax: " << difference << endl;
                if (difference<0.0){
                    difference = difference * -1.0;
                }

                dimensionedScalar redPercentage = Foam::exp(dimensionedScalar("",dimless,(-(difference*redCoef))));
                
                Info<< "-> redPercentage: " << redPercentage.value() << endl;
                Info<< "-> Previous deltaT = " <<  runTime.deltaTValue() << endl;

                if(difference > 0.05){  //while cmax < cs-0.05
                    if(cMax-cMaxPrev>0 && cMax > cs-1.0){   // c is increasing and close to saturation
                        newDeltaT = newDeltaT*(1.0-redPercentage.value());
                    }
                    else if(cMax-cMaxPrev>0 && cMax < cs-1.0){  //c is increasing and far from saturation
                        newDeltaT = newDeltaT+newDeltaT*(redPercentage.value());
                    }
                    else if(cMax-cMaxPrev<0 && cMax > cs-1.0){  //c is decreasing and close to saturation
                        newDeltaT = newDeltaT;
                    }
                    else if(cMax-cMaxPrev<0 && cMax < cs-1.0){  //c is decreasing and far from saturation
                        newDeltaT = newDeltaT+newDeltaT*(redPercentage.value());
                    }

                    runTime.setDeltaT( min(newDeltaT, maxDeltaT) );
                }
                else{       //when difference between cs and cMax is already too small

                    newDeltaT = minDeltaT;
                    runTime.setDeltaT(minDeltaT);
                    
                }
            }
            else{   //increments before first bands (before completing first increase) 

                    newDeltaT = newDeltaT*1.2;
                    Info<< "-> Previous deltaT = " <<  runTime.deltaTValue() << endl;
                    runTime.setDeltaT( min(newDeltaT, maxDeltaT) );
            }

            Info<< "-> New deltaT = " <<  runTime.deltaTValue() << endl;
            //-^
        }
    }

}
