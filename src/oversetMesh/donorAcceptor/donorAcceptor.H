/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::donorAcceptor

Description
    Class holds donor and acceptor data.  Used for searching and communications

SourceFiles

Author
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef donorAcceptor_H
#define donorAcceptor_H

#include "label.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class donorAcceptor;

inline bool operator==(const donorAcceptor&, const donorAcceptor&);
inline bool operator!=(const donorAcceptor&, const donorAcceptor&);

Istream& operator>>(Istream&, donorAcceptor&);
Ostream& operator<<(Ostream&, const donorAcceptor&);


/*---------------------------------------------------------------------------*\
                        Class donorAcceptor Declaration
\*---------------------------------------------------------------------------*/

class donorAcceptor
{
    // Private data

        // Acceptor side

            //- Acceptor cell number
            label acceptorCell_;

            //- Acceptor processor number
            label acceptorProcNo_;

            //- Acceptor location
            point acceptorPoint_;


        // Donor side

            //- Acceptor cell number
            label donorCell_;

            //- Donor processor number
            label donorProcNo_;


public:

    // Constructors

        //- Construct null
        inline donorAcceptor()
        {}

        //- Construct from acceptor data
        inline donorAcceptor
        (
            const label acceptorCell,
            const label acceptorProcNo,
            const point& acceptorPoint
        )
        :
            acceptorCell_(acceptorCell),
            acceptorProcNo_(acceptorProcNo),
            acceptorPoint_(acceptorPoint),
            donorCell_(-1),
            donorProcNo_(-1)
        {}

        //- Construct from Istream
        inline donorAcceptor(Istream& is)
        :
            acceptorCell_(readLabel(is)),
            acceptorProcNo_(readLabel(is)),
            acceptorPoint_(is),
            donorCell_(readLabel(is)),
            donorProcNo_(readLabel(is))
        {}

        //- Copy constructor - default


    //- Destructor - default


    // Member Functions

        // Acceptor side

            //- Return acceptor cell number
            label acceptorCell() const
            {
                return acceptorCell_;
            }

            //- Return access to acceptor cell number
            label& acceptorCell()
            {
                return acceptorCell_;
            }

            //- Return acceptor processor number
            label acceptorProcNo() const
            {
                return acceptorProcNo_;
            }

            //- Return access to acceptor processor number
            label& acceptorProcNo()
            {
                return acceptorProcNo_;
            }

            //- Return acceptor location
            const point& acceptorPoint() const
            {
                return acceptorPoint_;
            }


        // Donor side

            //- Has a donor been found?
            bool donorFound() const
            {
                return donorCell_ > -1;
            }

            //- Return donor cell number
            label donorCell() const
            {
                return donorCell_;
            }

            //- Return access to donor cell number
            label& donorCell()
            {
                return donorCell_;
            }

            //- Return donor processor number
            label donorProcNo() const
            {
                return donorProcNo_;
            }

            //- Return access to donor processor number
            label& donorProcNo()
            {
                return donorProcNo_;
            }


        // Edit

            //- Set hit: donor found
            void setDonor(const label& donorCell, const label& donorProcNo)
            {
                donorCell_ = donorCell;
                donorProcNo_ = donorProcNo;
            }


    // Member Operators

        void operator=(const donorAcceptor& rd)
        {
            acceptorCell_ = rd.acceptorCell_;
            acceptorProcNo_ = rd.acceptorProcNo_;
            acceptorPoint_ = rd.acceptorPoint_;

            donorCell_ = rd.donorCell_;
            donorProcNo_ = rd.donorProcNo_;
        }


    // Friend Operators

        friend bool operator==
        (
            const donorAcceptor& a,
            const donorAcceptor& b
        )
        {
            return
                a.acceptorCell_ == b.acceptorCell_
             && a.acceptorProcNo_ == b.acceptorProcNo_
             && a.acceptorPoint_ == b.acceptorPoint_
             && a.donorCell_ == b.donorCell_
             && a.donorProcNo_ == b.donorProcNo_;
        }

        friend bool operator!=
        (
            const donorAcceptor& a,
            const donorAcceptor& b
        )
        {
            return !(a == b);
        }


    // IOstream Operators

        //- Istream operator
        friend Istream& operator>>(Istream& is, donorAcceptor& rd)
        {
            is  >> rd.acceptorCell_
                >> rd.acceptorProcNo_
                >> rd.acceptorPoint_
                >> rd.donorCell_
                >> rd.donorProcNo_;

            return is;
        }

        //- Ostream operator
        friend Ostream& operator<<(Ostream& os, const donorAcceptor& rd)
        {
            os  << rd.acceptorCell_ << token::SPACE
                << rd.acceptorProcNo_ << token::SPACE
                << rd.acceptorPoint_ << token::SPACE
                << rd.donorCell_ << token::SPACE
                << rd.donorProcNo_;

            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
