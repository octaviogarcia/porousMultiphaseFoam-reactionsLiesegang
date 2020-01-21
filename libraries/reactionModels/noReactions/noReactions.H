/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::noReactions

Description
    A reaction model that does nothing.

SourceFiles
    noReactions.C

Developer
    Gabriel S. Gerlero

\*---------------------------------------------------------------------------*/

#ifndef noReactions_H
#define noReactions_H

#include "reactionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionModels
{

/*---------------------------------------------------------------------------*\
                   Class noReactions Declaration
\*---------------------------------------------------------------------------*/

class noReactions
: public reactionModel
{
public:

    //- Runtime type information
    TypeName("noReactions");

    // Constructors

    noReactions
    (
        basicMultiComponentMixture& composition,
        const dictionary& reactions
    );

    // Member Functions

    void correct(bool massConservative) override;

    tmp<fvScalarMatrix> reactionTerm(const label speciesi) override;

    bool needsSubcycling() const override;

    bool alwaysMassConservative() const override;

    void postTransport() override;

private:
    const PtrList<volScalarField>& Y_;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif