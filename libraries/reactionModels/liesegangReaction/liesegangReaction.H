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
    liesegangReaction.C

Developers
    Gabriel S. Gerlero, Pierre Horgue, Pablo A. Kler

\*---------------------------------------------------------------------------*/

#ifndef liesegangReaction_H
#define liesegangReaction_H

#include "reactionModel.H"
#include "specieCoeffs.H"

#include "loggerLevels.H"
const int LOG_LEVEL = 1;
LoggerLevel<1> InfoL1;
LoggerLevel<2> InfoL2;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionModels
{

/*---------------------------------------------------------------------------*\
                   Class liesegangReaction Declaration
\*---------------------------------------------------------------------------*/

class liesegangReaction
: public reactionModel
{
public:

    //- Runtime type information
    TypeName("liesegangReaction");

    // Constructors

    liesegangReaction
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

    void calcCellsInRadius();
    //- Add reaction with its order obtained from lhs exponents. Returns the order
    label addReaction 
    (
        const List<specieCoeffs>& lhs,
        const List<specieCoeffs>& rhs,
        const dimensionedScalar& k,
        const label liesegang
    );

    //- Add a second order reaction (first order on both a and b)
    void addSecondOrderReaction 
    (
        const List<specieCoeffs>& lhs,
        const List<specieCoeffs>& rhs,
        const dimensionedScalar& k,
        const label liesegang,
        label a,
        label b
    );

    //- Add a first order reaction (first order on a)
    void addFirstOrderReaction 
    (
        const List<specieCoeffs>& lhs,
        const List<specieCoeffs>& rhs,
        const dimensionedScalar& k,
        const label liesegang,
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

    void heaviside2InternalField(
		const DimensionedField<scalar, Foam::volMesh>& cField,
		const DimensionedField<scalar, Foam::volMesh>& csField
	);
		
    const PtrList<volScalarField>& Y_;
    List<dimensionedScalar> k0_; // Zeroth-order rate constants
    List<List<dimensionedScalar>> k1_; // First-order rate constants
    List<List<List<dimensionedScalar>>> k2_; // Second-order rate constants

    List<List<dimensionedScalar>> sk1; // First-order rate special constants 
    List<List<List<dimensionedScalar>>> sk2; // Second-order rate special constants

    bool massConservative_ = false;
    PtrList<fvScalarMatrix> massConservativeReactionTerms_;
    
    const fvMesh& mesh;

    label influencedSpecieK1=-1;
    label influencedSpecieK2=-1;

    bool flag1st=true;

    volScalarField binary;          //Field with check results on concentration levels (below rho (1) or above or equal (0))
    volScalarField heaviField;      //Field with results of applying Heaviside Step function through Foam::neg() to each cell
    volVectorField cellsSizes;
	
    dimensionedScalar rho;
    dimensionedScalar cradius;      //Value that indicates the cristalization radius = how many METERS to look around to check rho.
    dimensionedScalar cs_scalar = dimensionedScalar("",dimensionSet(0,-3,0,0,1,0,0),0.);
    label influencedSpecieHS;
    List<List<label>> inRadius;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
