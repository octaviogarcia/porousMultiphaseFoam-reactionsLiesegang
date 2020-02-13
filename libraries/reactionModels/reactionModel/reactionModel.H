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
    Foam::reactionModel

Description
    An abstract base class for reaction models.

SourceFiles
    reactionModel.C
    reactionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef reactionModel_H
#define reactionModel_H

#include "basicMultiComponentMixture.H"
#include "fvScalarMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class reactionModel Declaration
\*---------------------------------------------------------------------------*/

class reactionModel
{
public:

    //- Runtime type information
    TypeName("reactionModel");


    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        reactionModel,
        dictionary,
        (
            basicMultiComponentMixture& composition,
            const dictionary& reactions
        ),
        (composition, reactions)
    );


    // Selectors

    //- Return a reference to the selected dispersion model
    static autoPtr<reactionModel> New
    (
        basicMultiComponentMixture& composition,
        const dictionary& reactions,
        const word& defaultModel = "noReactions"
    );


    // Member Functions

    //- Whether there might be benefit in solving the entire set of transport 
    // equations more than once per time step.
    virtual bool needsSubcycling() const = 0;

    //- Whether this reaction model is always mass conservative regardless of 
    // the value of the massConservative parameter passed to correct(...)
    virtual bool alwaysMassConservative() const = 0;

    //- To be called before each time the set of transport equations is solved.
    virtual void correct(bool massConservative=false) = 0;

    //- To be called when constructing the transport equation for each species,
    // after the trasnport equation for the previous species (if any) has
    // already been solved. Returns the reaction term that will be plugged into
    // the equation.
    virtual tmp<fvScalarMatrix> reactionTerm(const label speciesi) = 0;

    //- To be called when there are no more transport equations to be solved
    // for the current time step.
    virtual void postTransport() = 0;

    virtual double getRedCoef()=0;

    virtual double getCsValue()=0;

    virtual int getInfluencedSpecieHS()=0;

    virtual void setPrevMaxC(scalar pcm)=0;

    virtual double getPrevMaxC()=0;

    //- Destructor
    virtual ~reactionModel() = default;

    //- Disallow copy construction & copy assignment
    reactionModel(const reactionModel&) = delete;
    void operator=(const reactionModel&) = delete;


protected:

    // Constructors

    reactionModel() = default;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
