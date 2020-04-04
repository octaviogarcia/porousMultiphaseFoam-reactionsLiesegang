{
    Cresidual = 0;

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
            eps * R * Saturation * fvm::ddt(C)
            + fvm::div(phi, C, "div(phi,C)")    //performance
            - fvm::laplacian(eps*Saturation*Deff, C, "laplacian(Deff,C)")
            ==
            - sourceTerm    //performance
            - eps * R * Saturation * fvm::Sp(lambda,C)    //performance
            + eps * R * Saturation * reactionTerm
        );

        C.storePrevIter();

        const solverPerformance& sp = CEqn.solve("C");

        actualOkStop = actualOkStop && (sp.nIterations()<1 ? true : false);
        
        dimensionedScalar vSmall("vSmall", C.dimensions(), VSMALL);
		
        Cresidual = max(Cresidual, gMax(mag(((C-C.prevIter())/(C.prevIter()+vSmall)))->internalField()));
    }
}
