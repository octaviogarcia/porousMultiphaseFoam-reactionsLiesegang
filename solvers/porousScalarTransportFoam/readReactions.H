const auto& reactionModelType = transportProperties.lookupOrDefault<word>("reactionModel", "noReactions");

const IOdictionary reactions
(
    IOobject
    (
        "reactions",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

autoPtr<reactionModel> reactionModel = reactionModel::New(composition, reactions);

scalar reactionTol = 0;

if(reactionModel->needsSubcycling())
{
	reactionTol = mesh.solutionDict().subOrEmptyDict("reactions").lookupType<scalar>("tolerance");
}

bool reactionMassConservative = mesh.solutionDict().subOrEmptyDict("reactions").lookupOrDefault<bool>("massConservative",false);
