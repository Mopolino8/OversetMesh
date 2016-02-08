/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::donorSuitability::donorSuitability

Description
    An abstract base class for donor suitability algorithms. This class is
    used by overlapFringe algorithm to optimise overlap given donor suitability
    function.

    Donor suitability criteria is given in isDonorSuitable member function. The
    threshold value represents the admissible relative difference between donor
    and acceptor donor suitability value. For example, if using cellVolumes
    donor suitability, specifying threshold = 0.1 means that the donor's cell
    volume needs to be within 10% of the acceptor's cell volume.

SourceFiles
    donorSuitability.C
    newAcceptorSet.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef donorSuitability_H
#define donorSuitability_H

#include "dictionary.H"
#include "labelList.H"
#include "donorAcceptorList.H"
#include "scalarField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class overlapFringe;

namespace donorSuitability
{

/*---------------------------------------------------------------------------*\
                         Class donorSuitability Declaration
\*---------------------------------------------------------------------------*/

class donorSuitability
{
    // Private data

        //- Reference to overlapFringe algorithm
        const overlapFringe& overlapFringe_;

        //- Subdictionary for acceptor set specification
        const dictionary coeffDict_;

        //- Donor Suitability Function (DSF), a list of scalar fields for each
        //  processor. In case of serial run, it is a list with a single
        //  entry. Each field is defined in cell centres and is used for overlap
        //  minimisation.
        List<scalarField> dsf_;

        //- Threshold criterion for donor suitability
        scalar threshold_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        donorSuitability(const donorSuitability&);

        //- Disallow default bitwise assignment
        void operator=(const donorSuitability&);


protected:

    // Protected Member Functions

        // Access functions

            //- Return access to the donor suitability function
            List<scalarField>& donorSuitabilityFunction()
            {
                return dsf_;
            }

            //- Return access to the threshold criterion
            scalar& threshold()
            {
                return threshold_;
            }


        // Constructor helper for derived classes

            //- Perform gather scatter of the donor suitability function for
            //  parallel runs
            void combineDonorSuitabilityFunction(const scalarField& localDsf);


public:

    //- Runtime type information
    TypeName("donorSuitability");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            donorSuitability,
            dictionary,
            (
                const overlapFringe& overlapFringeAlgorithm,
                const dictionary& dict
            ),
            (overlapFringeAlgorithm, dict)
        );


    // Constructors

        //- Construct from components
        donorSuitability
        (
            const overlapFringe& overlapFringeAlgorithm,
            const dictionary& dict
        );


    // Selectors

        //- Return an autoPtr to the selected donorSuitability
        static autoPtr<donorSuitability> New
        (
            const overlapFringe& overlapFringeAlgorithm,
            const dictionary& dict
        );


    //- Destructor
    // Note: the destructor is pure virtual to make this class abstract
    virtual ~donorSuitability() = 0;


    // Member Functions

        // Access

            //- Return a const reference to overlapFringe algorithm
            const overlapFringe& overlapFringeAlgorithm() const;

            //- Return coefficient dictionary
            const dictionary& coeffDict() const;


        // Overlap minimisation via Donor Suitability Function

            //- Is this donor/acceptor pair suitable
            bool isDonorSuitable
            (
                const donorAcceptor& daPair
            ) const
            {
                // If donor is not found for this acceptor, return false. If it
                // is found, return suitability criteria based on donor
                // suitability function
                if (!daPair.donorFound())
                {
                    return false;
                }
                else
                {
                    // Get donor suitability function in acceptor/donor pair
                    const scalar& dsfAcceptor =
                        dsf_[daPair.acceptorProcNo()][daPair.acceptorCell()];
                    const scalar& dsfDonor =
                        dsf_[daPair.donorProcNo()][daPair.donorCell()];

                    return
                        mag
                        (
                            (dsfDonor - dsfAcceptor)/
                            (dsfAcceptor + SMALL)
                        ) < threshold_;
                }
            }

            //- Perform overlap minimisation given donor acceptor list
            void rejectUnsuitableDonors
            (
                donorAcceptorList& acceptorCandidates
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace donorSuitability
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
