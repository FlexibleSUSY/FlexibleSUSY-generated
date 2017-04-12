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

// File generated at Wed 12 Apr 2017 12:38:57

#ifndef lowNMSSM_INPUT_PARAMETERS_H
#define lowNMSSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct lowNMSSM_input_parameters {
   double Qin;
   double M1Input;
   double M2Input;
   double M3Input;
   double AtInput;
   double AbInput;
   double ATauInput;
   double TanBeta;
   double ml1Input;
   double ml2Input;
   double ml3Input;
   double me1Input;
   double me2Input;
   double me3Input;
   double mq1Input;
   double mq2Input;
   double mq3Input;
   double md1Input;
   double md2Input;
   double md3Input;
   double mu1Input;
   double mu2Input;
   double mu3Input;
   double LambdaInput;
   double KappaInput;
   double ALambdaInput;
   double AKappaInput;
   double MuEffInput;

   lowNMSSM_input_parameters()
      : Qin(0), M1Input(0), M2Input(0), M3Input(0), AtInput(0), AbInput(0),
   ATauInput(0), TanBeta(0), ml1Input(0), ml2Input(0), ml3Input(0), me1Input(0)
   , me2Input(0), me3Input(0), mq1Input(0), mq2Input(0), mq3Input(0), md1Input(
   0), md2Input(0), md3Input(0), mu1Input(0), mu2Input(0), mu3Input(0),
   LambdaInput(0), KappaInput(0), ALambdaInput(0), AKappaInput(0), MuEffInput(0
   )

   {}

   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const lowNMSSM_input_parameters&);

} // namespace flexiblesusy

#endif
