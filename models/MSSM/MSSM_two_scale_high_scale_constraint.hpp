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

// File generated at Mon 23 Feb 2015 14:23:10

#ifndef MSSM_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H
#define MSSM_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H

#include "MSSM_high_scale_constraint.hpp"
#include "MSSM_input_parameters.hpp"
#include "two_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class MSSM;

class Two_scale;

template<>
class MSSM_high_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   MSSM_high_scale_constraint();
   MSSM_high_scale_constraint(MSSM<Two_scale>*, const MSSM_input_parameters&);
   virtual ~MSSM_high_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   double get_initial_scale_guess() const;
   void initialize();
   void set_input_parameters(const MSSM_input_parameters&);
   void set_scale(double); ///< fix unification scale (0 = unfixed)

private:
   double scale;
   double initial_scale_guess;
   MSSM<Two_scale>* model;
   MSSM_input_parameters inputPars;

   void update_scale();
};

} // namespace flexiblesusy

#endif
