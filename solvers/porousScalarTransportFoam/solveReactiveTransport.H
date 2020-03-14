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
            
            Info << "Reactive transport iteration : " << iterTransport << ", Residual = " << Cresidual << endl; 
            
            if(actualOkStop){
                Info<< "No more iterations needed..." << endl;
            }

        } while(Cresidual > reactionTol); 

        if(reactionMassConservative && !reactionModel->alwaysMassConservative())
        {
            Info << " Mass conservation:" << endl;
            reactionModel->correct(true);
            #include "CEqn.H"

            Info<< "Reactive transport mass conservation correction, Final residual = " << Cresidual << endl;

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

    scalar tauD=VGREAT;

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
        
        if(adjustTimeStep){
            auto Deff = composition.Deff(speciesi)[5].xx();

            if(Deff != 0.0){
                const auto& Cref = composition.Y(speciesi);
                double referenceVal[] = {500.0, 10.0, 10.0};
                const auto& maxCref = referenceVal[speciesi];
                auto maxmagGradC = max(mag(fvc::grad(Cref))).value();

                auto h = Cref.mesh().delta().ref()[5].x();
                Info << "h value: " << h << endl;
                auto hsq = pow(h,2);
                Info << "h squared: " << hsq << endl;
                auto k = kValueTauD;
                Info << "k value: " << k << endl;
                scalar tempTauD = 0.0;

                if(maxmagGradC>0.0){
                    auto fraction = maxCref/maxmagGradC;
                    tempTauD = k*(fraction * (hsq/Deff));
                }else{
                    tempTauD = VGREAT;
                }
                Info << "tauD specie: " << speciesi << " value: " << tempTauD << endl;
                tauD = tempTauD < tauD ? tempTauD : tauD;
            }
            else{

                scalar newDeltaT = runTime.deltaTValue();
                const auto& Dact = composition.Y(speciesi);
                const auto& Dprev = composition.Y(speciesi).oldTime();

                auto variation = max(Dact-Dprev).value();

                Info << endl << "tauD selected: " << tauD << endl;
                Info << "D Variation: " << variation << endl;

                if(variation > maxDCVariation){
                    Info<< "-> Max D variation exceded ..."<< endl;
                    breakLoop=100;      //amount of timesteps taken with constant deltaT to avoid a loop where deltaT is increased and decreased continuously, can be generalized

                    forAll(composition.Y(), speciesj){      // Reset fields to previous value.
                        composition.Y(speciesj) = composition.Y(speciesj).oldTime();
                    }

                    runTime.setTime(runTime.timeOutputValue() - runTime.deltaTValue(), runTime.timeIndex());
                    
                    newDeltaT *= 0.8;       // calculate new delta T to achieve less variation
                    Info<< "deltaT reduced 20% ..."<< endl << endl;

                    runTime.setTime(runTime.timeOutputValue() + min(newDeltaT, maxDeltaT), runTime.timeIndex());
                    
                    redoTimeStep=true;      // redo last timestep
                }
                else if (!breakLoop || (variation < maxDCVariation*0.75)){  // variation is acceptable to keep increasing deltaT (can be generalized)
                    if (tauD > (1.2 * newDeltaT)){
                        newDeltaT *= 1.2;
                        Info<< "deltaT increased 20% ..."<< endl << endl;
                    }else if (tauD < 1.2 * newDeltaT){
                        newDeltaT = tauD;
                        Info<< "using tauD ..."<< endl << endl;
                    }
                    redoTimeStep=false;
                }
                else{
                    breakLoop--;
                    Info<< "Decreasing loop counter. Actual value: " << breakLoop << endl;
                    if(!breakLoop){
                        Info<< "Loop broken ... " << endl;
                    }
                    redoTimeStep=false;
                }
                runTime.setDeltaT( min(newDeltaT, maxDeltaT) );
                Info<< "deltaT = " <<  runTime.deltaTValue() << endl << endl;
            }
        }
    }

}
