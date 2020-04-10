// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Fri 10 Apr 2020 17:34:00

#ifndef NUHMSSMNoFVHimalaya_INPUT_PARAMETERS_H
#define NUHMSSMNoFVHimalaya_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct NUHMSSMNoFVHimalaya_input_parameters {
   double TanBeta{};
   double Qin{};
   double M1{};
   double M2{};
   double M3{};
   double AtIN{};
   double AbIN{};
   double AtauIN{};
   double AcIN{};
   double AsIN{};
   double AmuonIN{};
   double AuIN{};
   double AdIN{};
   double AeIN{};
   double MuIN{};
   double mA2IN{};
   double ml11IN{};
   double ml22IN{};
   double ml33IN{};
   double me11IN{};
   double me22IN{};
   double me33IN{};
   double mq11IN{};
   double mq22IN{};
   double mq33IN{};
   double mu11IN{};
   double mu22IN{};
   double mu33IN{};
   double md11IN{};
   double md22IN{};
   double md33IN{};
   double Mlow{};


   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const NUHMSSMNoFVHimalaya_input_parameters&);

} // namespace flexiblesusy

#endif
