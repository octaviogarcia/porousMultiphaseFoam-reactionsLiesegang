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

\*---------------------------------------------------------------------------*/

#include "noReactions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionModels
{
defineTypeNameAndDebug(noReactions, 0);

addToRunTimeSelectionTable
(
    reactionModel,
    noReactions,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionModels::noReactions::noReactions
(
    basicMultiComponentMixture& composition,
    const dictionary& reactions
)
	:
	Y_(composition.Y())
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void Foam::reactionModels::noReactions::correct(bool massConsevative)
{}


Foam::tmp<Foam::fvScalarMatrix> Foam::reactionModels::noReactions::reactionTerm(const label speciesi)
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y_[speciesi], Y_[speciesi].dimensions()*dimVol/dimTime));

	return tSu;
}

bool Foam::reactionModels::noReactions::needsSubcycling() const
{
    return false;
}

bool Foam::reactionModels::noReactions::alwaysMassConservative() const
{
    return true;
}

void Foam::reactionModels::noReactions::postTransport()
{}


// ************************************************************************* //
