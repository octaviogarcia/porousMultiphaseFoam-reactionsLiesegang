{
    Cresidual = 0;

    //- update water flux
    phiWater = phi * fvc::interpolate(hwater);

    forAll(composition.Y(), speciesi)
    {
        auto& C = composition.Y(speciesi);
        const auto& R = composition.R(speciesi);
        const auto& Deff = composition.Deff(speciesi);
        const auto& lambda = composition.lambda(speciesi);
        const auto& sourceTerm = composition.sourceTerm(speciesi);
        const auto& reactionTerm = reactionModel->reactionTerm(speciesi);

        fvScalarMatrix CEqn
        (
            eps * R * fvm::ddt(hwater,C)
            + fvm::div(phiWater, C, "div(phi,C)")
            - fvm::laplacian(eps*hwater*Deff, C, "laplacian(Deff,C)")
            ==
            - sourceTerm * zScale
            - eps * R * hwater * fvm::Sp(lambda,C)
            - fvm::Sp(seepageTerm,C)
        );

        C.storePrevIter();

        CEqn.solve("C");

        dimensionedScalar vSmall("vSmall", C.dimensions(), VSMALL);

        Cresidual = max(Cresidual, gMax(mag(((C-C.prevIter())/(C.prevIter()+vSmall)))->internalField()));
    }
}

