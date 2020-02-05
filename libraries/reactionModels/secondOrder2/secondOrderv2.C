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

#include "secondOrderv2.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmSup.H"
#include "OStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionModels
{
defineTypeNameAndDebug(secondOrderv2, 0);

addToRunTimeSelectionTable
(
    reactionModel,
    secondOrderv2,
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
        FatalErrorIn("secondOrderv2.C")
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

// * * * * * * * * * * * * * * * * AUX  * * * * * * * * * * * * * * //
template<class Type, template<class> class PatchField, class GeoMesh>
void printField(const Foam::GeometricField<Type,PatchField,GeoMesh>& dfield){
	const auto& field = dfield.internalField();
    forAll(field, index){        
		if(index % 10 == 0){
	        Foam::Info << Foam::endl;
		}
		Foam::Info<< " " << field[index];
	}
	Foam::Info << Foam::endl;
}

//No estoy seguro si openfoam maneja floats o doubles
template<typename F>
F sigmoidAbs(F x,F steepness){
	x*=steepness;
	const F negOne_to_one = x/(1+abs(x));
	const F zero_to_one = (negOne_to_one + 1)/2.0;
	return zero_to_one;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionModels::secondOrderv2::secondOrderv2
(
    basicMultiComponentMixture& composition,
    const dictionary& reactions
)
	:
    Y_(composition.Y()),
    k0_(Y_.size()),
    k1_(Y_.size(), List<dimensionedScalar>(Y_.size())),
    k2_(Y_.size(), List<List<dimensionedScalar>>(Y_.size(), List<dimensionedScalar>(Y_.size()))),
    sk1(Y_.size(), List<dimensionedScalar>(Y_.size())),
    sk2(Y_.size(), List<List<dimensionedScalar>>(Y_.size(), List<dimensionedScalar>(Y_.size()))),
    binary(IOobject(
             	    word("binary"),
             	    Y_[0].mesh().time().timeName(), Y_[0].mesh(),
             		IOobject::NO_READ, IOobject::NO_WRITE
            	),
            	Y_[0].mesh(), dimensionedScalar("",dimless,1)),
    heaviField(binary),
    rho("", dimensionedScalar("",dimensionSet(0,-3,0,0,1,0,0),0.)),
    cs_scalar("", rho),
	//@HACK ver si hay alguna forma de iniciarlo mejor
	cradius("",
		dimensionedScalar(
			"",
			dimensionSet(0,1,0,0,0,0,0),
			Y_[0].mesh().delta().ref()[5].x()
		)
	),
	steepness(4096)
{
	// Asi se puede scar la cantidad de celdas, si podemos sacar la longitud del mesh o pasarlo como argumento
	// Seria mas facil
    // Info<< "Cells " << Y_[0].mesh().cells().size() << endl;
    //- Set the dimensions of all reaction coefficients
    forAll(Y_, speciesi)
    {
        k0_[speciesi].dimensions().reset(dimMoles/dimVol/dimTime);

        forAll(Y_, speciesj)
        {
            k1_[speciesi][speciesj].dimensions().reset(dimless/dimTime);
            sk1[speciesi][speciesj].dimensions().reset(dimless/dimTime);

            forAll(Y_, speciesk)
            {
                k2_[speciesi][speciesj][speciesk].dimensions().reset(dimVol/dimMoles/dimTime);
                sk2[speciesi][speciesj][speciesk].dimensions().reset(dimVol/dimMoles/dimTime);
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
        const dimensionedScalar& liesegang = reaction.lookupOrDefault("liesegang",dimensionedScalar("",dimless,0.));

        dimensionedScalar rho_value = reaction.lookupOrDefault("rho",dimensionedScalar("",dimensionSet(0,-3,0,0,1,0,0),0.));
		
        if(rho_value.value()!=0){
			rho = rho_value;
        }
		
		dimensionedScalar cradius_value = reaction.lookupOrDefault("cradius",dimensionedScalar("",dimensionSet(0,1,0,0,0,0,0),0.));
		
		if(cradius_value.value()!=0){
			cradius = cradius_value;
		}
		
		dimensionedScalar steepness_value = reaction.lookupOrDefault("steepness",dimensionedScalar("",dimensionSet(0,0,0,0,0,0,0),0.));
		if(steepness_value.value() != 0){
			steepness = steepness_value.value();
		}

        Info<<"   rho Loaded for reaction: " << rho.value() << endl;
		Info<<"   cradius Loaded for reaction: " << cradius.value() << endl;

        auto order = addReaction(lhs, rhs, kf, liesegang.value());

        Info<< "    [->]: kf " << kf.value() << ", order " << order << endl;
        Info<< "    Liesegang: " << liesegang << endl;

        dimensionedScalar cs_value = reaction.lookupOrDefault("cs",dimensionedScalar("",dimensionSet(0,-3,0,0,1,0,0),0.));
        
        if(cs_value.value()!=0){
			cs_scalar = cs_value;
        }

        Info<< "    cs Loaded for reaction: " << cs_scalar.value() << endl; 
                
        if(reaction.found("kb")) 
        {
            const dimensionedScalar& kb = reaction.lookup("kb");
            order = addReaction(rhs, lhs, kb, liesegang.value()); //- Add the inverse reaction
            Info<< "    [<-]: kb " << kb.value() << ", order " << order << endl;
        }

        Info<< "}" << endl << endl;

        Info<< "Loaded Reaction: " << endl;     //TODO: delete this
        OStringStream reactionOSStream;         //<-|
        Info<< specieCoeffs::reactionStr(       //  |
                reactionOSStream                //  |
                ,composition.species()          //  |
                ,lhs                            //  |
                ,rhs) << endl << endl;          //<-| 
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
                FatalErrorIn("secondOrderv2.C")
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
                FatalErrorIn("secondOrderv2.C")
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


void Foam::reactionModels::secondOrderv2::correct(bool massConservative)
{
    if(!flag1st){
        return;
    }
	
    Info << "rho value:" << rho.value() << endl << endl;

	// Para nuestro caso, influenced species K1 y K2 son el lado derecho
	// de las reacciones C = D osea D.
	if(influencedSpecieK1 != influencedSpecieK2){
		FatalErrorIn("secondOrderv2.C")
		<< "InfluencedSpecieK1,K2 valores distintos."
		<< abort(FatalError);
	}
	volScalarField cs( // Field with c* amount for each cell
		IOobject(
			word("cs"),
			Y_[0].mesh().time().timeName(), 
			Y_[0].mesh(),
			IOobject::NO_READ, 
			IOobject::NO_WRITE
		),
		Y_[0].mesh(), 
		cs_scalar
	);
	
	binary.ref().field() = 1;
	auto binaryField = binary.ref();
    auto csField = cs.ref();
	const auto& fieldTargetMass = Y_[influencedSpecieK1].internalField();
	
	//@hack;
	// [cells] = RADIUS [m] / DELTAX [m/cell]
	const double r = (cradius / Y_[0].mesh().delta().ref()[5].x()).value();
	Info << "Cell check radius " << r << endl;
	//if r = 1.25
	const int int_r = round(r); // 2
	
	forAll(fieldTargetMass, cell){
		//if r = 1.25 => int_r = 2 => we check [cell-2,cell-1,cell,cell+1,cell+2]
		const int begin = max(cell - int_r, 0);
		const int end = min(cell + int_r, fieldTargetMass.size() - 1);
		for(int i = begin;i<=end && csField[cell] > 0.0;i++){
			//converts [cell-2,cell-1,cell,cell+1,cell+2] to [2,1,0,1,2]
			const int inside_pos = abs(i-cell);
			const double overTargetMass = fieldTargetMass[i] - rho.value();
			const int is_not_center = inside_pos != 0;
			const double s = is_not_center*sigmoidAbs(overTargetMass,steepness);
			
			csField[cell] = (1-s)*csField[cell];
		}
		binary[cell] = sigmoidAbs(rho.value() - fieldTargetMass[cell],steepness);
	}
	
	auto cField = Y_[influencedSpecieHS].internalField();
	heaviside2InternalField(cField, csField);
	
	//printField(heaviField);
	flag1st=false;

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


Foam::tmp<Foam::fvScalarMatrix> Foam::reactionModels::secondOrderv2::reactionTerm
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

bool Foam::reactionModels::secondOrderv2::needsSubcycling() const
{
    return true;
}

bool Foam::reactionModels::secondOrderv2::alwaysMassConservative() const
{
    return false;
}


void Foam::reactionModels::secondOrderv2::postTransport()
{
    flag1st=true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::reactionModels::secondOrderv2::addReaction
(
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs,
    const dimensionedScalar& k,
    const label liesegang
)
{
    labelList orderIndices;

    forAll(lhs, i)
    {
        auto exponent = floor(lhs[i].exponent);
        if(lhs[i].exponent != exponent)
        {
            FatalErrorIn("secondOrderv2.C")
                << "Non-integer exponents are not supported by this reaction model"
                << abort(FatalError);
        }

        //- Append this species' index as many times as the exponent dictates
        orderIndices.setSize(orderIndices.size() + exponent, lhs[i].index);
    }

    forAll(rhs, i)
    {
        auto exponent = floor(rhs[i].exponent);
        if(rhs[i].exponent != exponent)
        {
            FatalErrorIn("secondOrderv2.C")
                << "Non-integer exponents are not supported by this reaction model"
                << abort(FatalError);
        }

        //- Append this species' index as many times as the exponent dictates
        orderIndices.setSize(orderIndices.size() + exponent, rhs[i].index);
    }

    switch(orderIndices.size())
    {
        case 2:
            addSecondOrderReaction(lhs, rhs, k, liesegang, orderIndices[0], orderIndices[1]);
            break;

        case 1:
            addFirstOrderReaction(lhs, rhs, k, liesegang, orderIndices[0]);
            break;

        case 0:
            addZerothOrderReaction(lhs, rhs, k);
            break;

        default:
            FatalErrorIn("secondOrderv2.C")
                << "Reaction of order " << orderIndices.size()
                << " not supported by this reaction model. "
                << "Only reactions of orders 0, 1 and 2 are supported"
                << abort(FatalError);
            break;
    }

    return orderIndices.size();
}

void Foam::reactionModels::secondOrderv2::addSecondOrderReaction
(
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs,
    const dimensionedScalar& k,
    const label liesegang,
    label a,
    label b
)
{
    if(liesegang == 2){
        forAll(lhs, i)
        {
            sk2[lhs[i].index][a][b] -= lhs[i].stoichCoeff * k;
        }

        forAll(rhs, i)
        {
            sk2[rhs[i].index][a][b] += rhs[i].stoichCoeff * k;
            influencedSpecieK2 = rhs[i].index;
        }
    }
    else{
        forAll(lhs, i)
        {
            k2_[lhs[i].index][a][b] -= lhs[i].stoichCoeff * k;
        }

        forAll(rhs, i)
        {
            k2_[rhs[i].index][a][b] += rhs[i].stoichCoeff * k;
        }
    }
}

void Foam::reactionModels::secondOrderv2::addFirstOrderReaction
(
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs,
    const dimensionedScalar& k,
    const label liesegang,
    label a
)
{
    if(liesegang == 1){
        forAll(lhs, i){
            sk1[lhs[i].index][a] -= lhs[i].stoichCoeff * k;
            influencedSpecieHS = lhs[i].index;
        }
        forAll(rhs, i){
            sk1[rhs[i].index][a] += rhs[i].stoichCoeff * k;
            influencedSpecieK1 = rhs[i].index;
        }
    }
    else{
        forAll(lhs, i){
            k1_[lhs[i].index][a] -= lhs[i].stoichCoeff * k;
        }
        forAll(rhs, i){
            k1_[rhs[i].index][a] += rhs[i].stoichCoeff * k;
        }
    }
}

void Foam::reactionModels::secondOrderv2::addZerothOrderReaction
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


Foam::tmp<Foam::fvScalarMatrix> Foam::reactionModels::secondOrderv2::computeReactionTerm
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
            term += fvm::Sp(sk1[speciesi][speciesi] * heaviField * binary, Y_[speciesi]);

            forAll(Y_, speciesk)
            {
                term += fvm::Sp(k2_[speciesi][speciesi][speciesk]*Y_[speciesk], Y_[speciesi]);
                term += fvm::Sp(sk2[speciesi][speciesi][speciesk]*Y_[speciesk]*binary, Y_[speciesi]);
            }
        }
        else
        {
            term += k1_[speciesi][speciesj]*Y_[speciesj];
            term += sk1[speciesi][speciesj]*Y_[speciesj] * heaviField * binary;

            forAll(Y_, speciesk)
            {
                if(implicit && speciesk == speciesi)
                {
                    term += fvm::Sp(k2_[speciesi][speciesj][speciesi]*Y_[speciesj], Y_[speciesi]);
                    term += fvm::Sp(sk2[speciesi][speciesj][speciesi]*Y_[speciesj]*binary, Y_[speciesi]);
                }
                else
                {
                    term += k2_[speciesi][speciesj][speciesk]*Y_[speciesj]*Y_[speciesk];
                    term += sk2[speciesi][speciesj][speciesk]*Y_[speciesj]*Y_[speciesk]*binary;
                }
            }
        }
    }

    return tTerm;
}

void Foam::reactionModels::secondOrderv2::heaviside2InternalField
(
    const DimensionedField<scalar, Foam::volMesh>& cField,
    const DimensionedField<scalar, Foam::volMesh>& csField
)
{
	auto& hIntfield = heaviField.ref();
    forAll(cField,cell){
		// Esto es literalmente x > 0? 1 : 0
        hIntfield[cell]=Foam::pos(cField[cell]-csField[cell]);
		// Con un sigmoide no funciona, supongo que tiene que ser muy steep o no funciona
		// Entonces creo que es mejor dejarlo con un escalon nomas.
		//hIntfield[cell]=sigmoidAbs(cField[cell]-csField[cell],steepness);
    }
}

// ************************************************************************* //
