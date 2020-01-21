{
    Cresidual = 0;

    forAll(composition.Y(), speciesi)
    {
        auto& C = composition.Y(speciesi);
        const auto& R = composition.R(speciesi);
        const auto& Deff = composition.Deff(speciesi);
        const auto& lambda = composition.lambda(speciesi);
        const auto& sourceTerm_tracer = composition.sourceTerm(speciesi);
        const auto& reactionTerm = reactionModel->reactionTerm(speciesi);

        fvScalarMatrix CEqn
        (
            R * fvm::ddt(theta,C)
            + fvm::div(phi, C, "div(phi,C)")
            - fvm::laplacian(theta*Deff, C, "laplacian(Deff,C)")
            ==
            - sourceTerm_tracer
            - R * theta * fvm::Sp(lambda,C)
        );

        C.storePrevIter();

        CEqn.solve("C");

        dimensionedScalar vSmall("vSmall", C.dimensions(), VSMALL);

        Cresidual = max(Cresidual, gMax(mag(((C-C.prevIter())/(C.prevIter()+vSmall)))->internalField()));
    }
}
