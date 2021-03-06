{    
    //- update dispersion tensor coefficients and source terms
    forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateValue(runTime);
    forAll(sourceEventList,sourceEventi) sourceEventList[sourceEventi]->updateValue(runTime);
    composition.correct(U, Saturation, eps);

    scalar Cresidual;
    bool actualOkStop = true;

    if(reactionModel->needsSubcycling())
    {
        scalar iterTransport = 0;
        do
        {
            ++iterTransport;

            actualOkStop = true;
            reactionModel->correct(false);
            #include "CEqn.H"
            
            InfoL2 << "Reactive transport iteration : " << iterTransport << ", Residual = " << Cresidual << endl; 
            
            if(actualOkStop){
                InfoL2<< "No more iterations needed..." << endl;
            }

        } while(Cresidual > reactionTol); 

        if(reactionMassConservative && !reactionModel->alwaysMassConservative())
        {
            InfoL2 << " Mass conservation:" << endl;
            reactionModel->correct(true);
            #include "CEqn.H"

            InfoL2<< "Reactive transport mass conservation correction, Final residual = " << Cresidual << endl;

            if(Cresidual > reactionTol)
            {
                Warning() << "Final residual after mass conservation correction exceeds tolerance" << endl;
            }
        }
    }
    else
    {
        reactionModel->correct(reactionMassConservative);
        #include "CEqn.H"
    }

    reactionModel->postTransport();

    scalar tauD = VGREAT;

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
        InfoL2<< "Concentration Min(" << speciesName << ") = " << min(C).value() 
            << " Max(" << speciesName << ") = " << max(C.internalField()).value()
            << " mass(" << speciesName << ") = " << fvc::domainIntegrate(R*C*Saturation*eps).value()
            << " dCmax = " << dCdTmax[speciesi]*runTime.deltaTValue()
            << endl;
        
        if(adjustTimeStep){
            auto Deff = composition.Deff(speciesi)[5].xx();

            if(Deff != 0.0){
                const auto& maxCref = referenceVal[speciesi];
                const auto& Cref = composition.Y(speciesi);

                InfoL2 << setprecision(15);


                // formula: (k*maxCref*h^2)/(Deff*|GradC|)
                scalar tempTauD = 0.0;

                surfaceScalarField tauDField(tauDField_base);
                surfaceScalarField magGradCParallel = mag(fvc::interpolate(fvc::grad(Cref)) ^ unitaryNormalToFace).ref() + dimensionedScalar("", dimensionSet(0,-4,0,0,1,0,0), 1.0e-50);
                tauDField *= maxCref/(magGradCParallel*Deff);

                tempTauD = min(tauDField).value();


                InfoL2<< "tauD vectorized specie: " << speciesi << " value: " << tempTauD << endl;
                tauD = min(tempTauD, tauD);       // between species looks for the smallest
            }
            else //if(speciesi == composition.species().size()-1)
            {
                scalar newDeltaT = runTime.deltaTValue();
                const auto& Dact = composition.Y(speciesi);
                const auto& Dprev = composition.Y(speciesi).oldTime();

                auto variation = max(Dact-Dprev).value();
                
                InfoL2 << endl << "tauD selected: " << tauD << endl;
                InfoL2 << "D Variation: " << variation << endl;

                if(variation > maxDCVariation){
                    InfoL2<< "-> Max D variation exceded ..."<< endl;
                    breakLoop=100;      //amount of timesteps taken with constant deltaT to avoid a loop where deltaT is increased and decreased continuously, can be generalized

                    forAll(composition.Y(), speciesj){      // Reset fields to previous value.
                        composition.Y(speciesj) = composition.Y(speciesj).oldTime();
                    }

                    runTime.setTime(runTime.timeOutputValue() - runTime.deltaTValue(), runTime.timeIndex());
                    
                    newDeltaT *= 1.0-deltaTDecrease;       // calculate new deltaT to achieve less variation....if deltaTDecrease = 0.2 => newDeltaT *= 0.8;
                    InfoL2<< "deltaT reduced 20% ..."<< endl << endl;

                    runTime.setTime(runTime.timeOutputValue() + min(newDeltaT, maxDeltaT), runTime.timeIndex());
                    
                    redoTimeStep=true;      // redo last timestep
                }
                else if (!breakLoop || (variation < maxDCVariation*0.75)){  // variation is acceptable to keep increasing deltaT (can be generalized)
                                                                                    //if maxDeltaTIncrease=0.2
                    if (tauD > ((1.0+maxDeltaTIncrease) * newDeltaT)){              //   if tauD > (1.0+0.2) * newDeltaT
                        newDeltaT *= (1.0+maxDeltaTIncrease);                       //      newDeltaT *= (1.0+0.2);
                        InfoL2<< "deltaT increased 20% ..."<< endl << endl;
                    }else if (tauD < (1.0+maxDeltaTIncrease) * newDeltaT){          //   if tauD < (1.0+0.2) * newDeltaT
                        newDeltaT = tauD;
                        InfoL2<< "using tauD ..."<< endl << endl;
                    }
                    redoTimeStep=false;
                }
                else{
                    breakLoop--;
                    InfoL2<< "Decreasing loop counter. Actual value: " << breakLoop << endl;
                    if(!breakLoop){
                        InfoL2<< "Loop broken ... " << endl;
                    }
                    redoTimeStep=false;
                }
                runTime.setDeltaT( min(newDeltaT, maxDeltaT) );
                InfoL2<< "deltaT = " <<  runTime.deltaTValue() << endl << endl;
            }
        }
    }

}
