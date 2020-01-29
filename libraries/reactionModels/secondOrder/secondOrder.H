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
    Foam::firstOrder

Description
    General second-order non-equilibrium reaction model.

SourceFiles
    secondOrder.C

Developers
    Gabriel S. Gerlero, Pierre Horgue, Pablo A. Kler

\*---------------------------------------------------------------------------*/

#ifndef secondOrder_H
#define secondOrder_H

#include "reactionModel.H"
#include "specieCoeffs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionModels
{

/*---------------------------------------------------------------------------*\
                   Class fristOrder Declaration
\*---------------------------------------------------------------------------*/

class secondOrder
: public reactionModel
{
public:

    //- Runtime type information
    TypeName("secondOrder");

    // Constructors

    secondOrder
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
    //- Add reaction with its order obtained from lhs exponents. Returns the order
    label addReaction 
    (
        const List<specieCoeffs>& lhs,
        const List<specieCoeffs>& rhs,
        const dimensionedScalar& k
    );

    //- Add a second order reaction (first order on both a and b)
    void addSecondOrderReaction 
    (
        const List<specieCoeffs>& lhs,
        const List<specieCoeffs>& rhs,
        const dimensionedScalar& k,
        label a,
        label b
    );
    //- Add a second order reaction (second order on a)
    void addSecondOrderReaction 
    (
        const List<specieCoeffs>& lhs,
        const List<specieCoeffs>& rhs,
        const dimensionedScalar& k,
        label a
    );

    //- Add a first order reaction (first order on a)
    void addFirstOrderReaction 
    (
        const List<specieCoeffs>& lhs,
        const List<specieCoeffs>& rhs,
        const dimensionedScalar& k,
        label a
    );

    //- Add a zeroth order reaction
    void addZerothOrderReaction
    (
        const List<specieCoeffs>& lhs,
        const List<specieCoeffs>& rhs,
        const dimensionedScalar& k
    ); 

    tmp<fvScalarMatrix> computeReactionTerm
    (
        const label speciesi,
        bool implicit
    ) const;

    const PtrList<volScalarField>& Y_;
    List<dimensionedScalar> k0_; // Zeroth-order rate constants
    List<List<dimensionedScalar>> k1_; // First-order rate constants
    List<List<List<dimensionedScalar>>> k2_; // Second-order rate constants
    List<List<dimensionedScalar>> sk1_; // First-order rate constants
    List<List<List<dimensionedScalar>>> sk2_; // Second-order rate constants

    bool massConservative_ = false;
    PtrList<fvScalarMatrix> massConservativeReactionTerms_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif