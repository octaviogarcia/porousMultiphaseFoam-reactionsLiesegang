{
    //- update dispersion tensor coefficients and source terms
    forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateValue(runTime);
    forAll(tracerSourceEventList,tracerSourceEventi) tracerSourceEventList[tracerSourceEventi]->updateValue(runTime);
    composition.correct(Utheta, theta);

    scalar Cresidual;
    if(reactionModel->needsSubcycling())
    {
        scalar iterTransport = 0;
        do
        {
            ++iterTransport;

            reactionModel->correct(false);
            #include "CEqn.H"
            
            Info<< "Reactive transport iteration : " << iterTransport << ", Residual = " << Cresidual << endl; 

        } while(Cresidual > reactionTol);

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
            << " Max(" << speciesName << ") = " << max(C).value()
            << " d" << speciesName << "dTmax = " << dCdTmax[speciesi]
            << " mass(" << speciesName << ") = " << fvc::domainIntegrate(C*theta).value()
            << " dCmax = " << dCdTmax[speciesi]*runTime.deltaTValue()
            << endl;

    }

}
