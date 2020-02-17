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

#include "secondOrder.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionModels
{
defineTypeNameAndDebug(secondOrder, 0);

addToRunTimeSelectionTable
(
    reactionModel,
    secondOrder,
    dictionary
);
} // End namespace reactionModels


namespace
{
//- Check if a reacting species has its concentration in dimMoles/dimVol
void checkDims 
(
    const Foam::basicMultiComponentMixture& composition,
    const Foam::label speciesi
)
{
    static const auto expected = dimMoles/dimVol;

    if(composition.Y(speciesi).dimensions() != expected)
    {
        FatalErrorIn("secondOrder.C")
            << "Species "
            << composition.species()[speciesi]
            << " appears in a reaction but its "
            << "concentration field has dimensions "
            << composition.Y(speciesi).dimensions()
            << ". Expected dimensions: "
            << expected
            << " (= moles/vol)"
            << abort(FatalError);
    }
}

} // End unnamed namespace
} // End namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionModels::secondOrder::secondOrder
(
    basicMultiComponentMixture& composition,
    const dictionary& reactions
)
	:
    Y_(composition.Y()),
    k0_(Y_.size()),
    k1_(Y_.size(), List<dimensionedScalar>(Y_.size())),
    k2_(Y_.size(), List<List<dimensionedScalar>>(Y_.size(), List<dimensionedScalar>(Y_.size()))),
    sk1_(Y_.size(), List<dimensionedScalar>(Y_.size())),
    sk2_(Y_.size(), List<List<dimensionedScalar>>(Y_.size(), List<dimensionedScalar>(Y_.size())))
{
    //- Set the dimensions of all reaction coefficients
    forAll(Y_, speciesi)
    {
        k0_[speciesi].dimensions().reset(dimMoles/dimVol/dimTime);

        forAll(Y_, speciesj)
        {
            k1_[speciesi][speciesj].dimensions().reset(dimless/dimTime);

            forAll(Y_, speciesk)
            {
                k2_[speciesi][speciesj][speciesk].dimensions().reset(dimVol/dimMoles/dimTime);
            }
        }
    }

    //- Read list of reactions
    auto reactionList = reactions.subOrEmptyDict("reactions");

    Info<< "Reading reactions..." << endl << endl;
    
    forAllConstIter(dictionary, reactionList, iter)
    {
        const auto& reactionName = iter().keyword();
        const auto& reaction = reactionList.subDict(reactionName);

        Info<< "Reaction " << reactionName << endl
            << "{" << endl
            << "    " << reaction.lookupType<string>("reaction") << endl;

        List<specieCoeffs> lhs;
        List<specieCoeffs> rhs;

        specieCoeffs::setLRhs
        (
            IStringStream(reaction.lookup("reaction"))(),
            composition.species(),
            lhs,
            rhs
        );

        forAll(lhs, i)
        {
            checkDims(composition, lhs[i].index);
        }

        forAll(rhs, i)
        {
            checkDims(composition, rhs[i].index);
        } 

        const dimensionedScalar& kf = reaction.lookup("kf");

        auto order = addReaction(lhs, rhs, kf);

        Info<< "    [->]: kf " << kf.value() << ", order " << order << endl;

        if(reaction.found("kb")) 
        {
            const dimensionedScalar& kb = reaction.lookup("kb");
            order = addReaction(rhs, lhs, kb); //- Add the inverse reaction

            Info<< "    [<-]: kb " << kb.value() << ", order " << order << endl;
        }

        Info<< "}" << endl << endl;
    }

    //- If no reactions list found, read any reaction rate constants given with each species
    if(!reactions.found("reactions"))
    {
        Info<< "===> No reactions found. Reading rate constants directly..." << endl;

        forAll(Y_, speciesi)
        {
            const auto& speciesDict = reactions.optionalSubDict(composition.species()[speciesi]);

            scalar K0 = speciesDict.lookupOrDefault("K0", 0);

            scalarList K1 = speciesDict.lookupOrDefault("K1", scalarList());
            if (K1.empty())
            {
                K1.resize(Y_.size(), 0);
            }
            else if (K1.size() != Y_.size())
            {
                FatalErrorIn("secondOrder.C")
                    << "K1 list found for species "
                    << composition.species()[speciesi]
                    << " with wrong size "
                    << K1.size()
                    << ". Expected size: "
                    << Y_.size()
                    << abort(FatalError);
            }

            scalarList K2 = speciesDict.lookupOrDefault("K2", scalarList());
            if (K2.empty() || Y_.size()*(Y_.size()+1)/2)
            {
                K2.resize(Y_.size()*Y_.size(), 0);
            }
            else if (K2.size() != Y_.size()*Y_.size())
            {
                FatalErrorIn("secondOrder.C")
                    << "K2 list found for species "
                    << composition.species()[speciesi]
                    << " found with wrong size " << K2.size()
                    <<  ". Expected size: either "
                    << Y_.size()*(Y_.size()+1)/2
                    <<  " or "
                    << Y_.size()*Y_.size()
                    << abort(FatalError);
            }

            if(K0)
            {
                checkDims(composition, speciesi);
                //Note: negative sign for compatibility with other code
                k0_[speciesi].value() = -K0;
            }

            forAll(Y_, speciesj)
            {
                auto value = K1[speciesj]; 

                if(value)
                {
                    checkDims(composition, speciesi);
                    checkDims(composition, speciesj);
                    //Note: negative sign for compatibility with other code
                    k1_[speciesi][speciesj].value() = -value; 
                }

                forAll(Y_, speciesk)
                {
                    auto value = K2[speciesj*Y_.size() + speciesk];

                    if(value)
                    {
                        checkDims(composition, speciesi);
                        checkDims(composition, speciesj);
                        checkDims(composition, speciesk);
                        //Note: negative sign for compatibility with other code
                        k2_[speciesi][speciesj][speciesk].value() = -value;
                    }
                }
            }
        }

        Info<< "Finished reading rate constants directly." << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void Foam::reactionModels::secondOrder::correct(bool massConservative)
{
    massConservative_ = massConservative;

    if(massConservative_)
    {
        //-To be mass conservative, reaction terms need to be precomputed and fully explicit
        massConservativeReactionTerms_.resize(Y_.size());

        forAll(Y_, speciesi) 
        {
            massConservativeReactionTerms_.set
            (
                speciesi,
                computeReactionTerm(speciesi, false)
            );
        }
    }
}


Foam::tmp<Foam::fvScalarMatrix> Foam::reactionModels::secondOrder::reactionTerm
(
    const label speciesi
)
{
    if(massConservative_)
    {
        return massConservativeReactionTerms_[speciesi];
    }
    else
    {
        return computeReactionTerm(speciesi, true);
    }
}

bool Foam::reactionModels::secondOrder::needsSubcycling() const
{
    return true;
}

bool Foam::reactionModels::secondOrder::alwaysMassConservative() const
{
    return false;
}


void Foam::reactionModels::secondOrder::postTransport()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::reactionModels::secondOrder::addReaction
(
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs,
    const dimensionedScalar& k
)
{
    labelList orderIndices;

    forAll(lhs, i)
    {
        auto exponent = floor(lhs[i].exponent);
        if(lhs[i].exponent != exponent)
        {
            FatalErrorIn("secondOrder.C")
                << "Non-integer exponents are not supported by this reaction model"
                << abort(FatalError);
        }

        //- Append this species' index as many times as the exponent dictates
        orderIndices.setSize(orderIndices.size() + exponent, lhs[i].index);
    }

    switch(orderIndices.size())
    {
        case 2:
            addSecondOrderReaction(lhs, rhs, k, orderIndices[0], orderIndices[1]);
            break;

        case 1:
            addFirstOrderReaction(lhs, rhs, k, orderIndices[0]);
            break;

        case 0:
            addZerothOrderReaction(lhs, rhs, k);
            break;

        default:
            FatalErrorIn("secondOrder.C")
                << "Reaction of order " << orderIndices.size()
                << " not supported by this reaction model. "
                << "Only reactions of orders 0, 1 and 2 are supported"
                << abort(FatalError);
            break;
    }

    return orderIndices.size();
}

void Foam::reactionModels::secondOrder::addSecondOrderReaction
(
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs,
    const dimensionedScalar& k,
    label a,
    label b
)
{
    forAll(lhs, i)
    {
        k2_[lhs[i].index][a][b] -= lhs[i].stoichCoeff * k;
    }

    forAll(rhs, i)
    {
        k2_[rhs[i].index][a][b] += rhs[i].stoichCoeff * k;
    }
}

void Foam::reactionModels::secondOrder::addSecondOrderReaction
(
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs,
    const dimensionedScalar& k,
    label a
)
{
    addSecondOrderReaction(lhs, rhs, k, a, a);
}

void Foam::reactionModels::secondOrder::addFirstOrderReaction
(
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs,
    const dimensionedScalar& k,
    label a
)
{
    forAll(lhs, i)
    {
        k1_[lhs[i].index][a] -= lhs[i].stoichCoeff * k;
    }

    forAll(rhs, i)
    {
        k1_[rhs[i].index][a] += rhs[i].stoichCoeff * k;
    }
}

void Foam::reactionModels::secondOrder::addZerothOrderReaction
(
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs,
    const dimensionedScalar& k
)
{
    forAll(lhs, i)
    {
        k0_[lhs[i].index] -= lhs[i].stoichCoeff * k;
    }

    forAll(rhs, i)
    {
        k0_[rhs[i].index] += rhs[i].stoichCoeff * k;
    }
}


Foam::tmp<Foam::fvScalarMatrix> Foam::reactionModels::secondOrder::computeReactionTerm
(
    const label speciesi,
    bool implicit
) const
{
    tmp<fvScalarMatrix> tTerm(new fvScalarMatrix(Y_[speciesi], Y_[speciesi].dimensions()*dimVol/dimTime));
    fvScalarMatrix& term = tTerm.ref();

    term += k0_[speciesi];

    forAll(Y_, speciesj)
    {
        if(implicit && speciesj == speciesi)
        {
            term += fvm::Sp(k1_[speciesi][speciesi], Y_[speciesi]);

            forAll(Y_, speciesk)
            {
                term += fvm::Sp(k2_[speciesi][speciesi][speciesk]*Y_[speciesk], Y_[speciesi]);
            }
        }
        else
        {
            term += k1_[speciesi][speciesj]*Y_[speciesj];

            forAll(Y_, speciesk)
            {
                if(implicit && speciesk == speciesi)
                {
                    term += fvm::Sp(k2_[speciesi][speciesj][speciesi]*Y_[speciesj], Y_[speciesi]);
                }
                else
                {
                    term += k2_[speciesi][speciesj][speciesk]*Y_[speciesj]*Y_[speciesk];
                }
            }
        }
    }

    return tTerm;
}

// ************************************************************************* //
