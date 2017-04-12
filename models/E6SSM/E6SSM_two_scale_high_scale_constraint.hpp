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

// File generated at Wed 12 Apr 2017 11:51:30

#ifndef E6SSM_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H
#define E6SSM_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H

#include "E6SSM_high_scale_constraint.hpp"
#include "E6SSM_input_parameters.hpp"
#include "two_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class E6SSM;

class Two_scale;

template<>
class E6SSM_high_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   E6SSM_high_scale_constraint();
   E6SSM_high_scale_constraint(E6SSM<Two_scale>*);
   virtual ~E6SSM_high_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   double get_initial_scale_guess() const;
   const E6SSM_input_parameters& get_input_parameters() const;
   E6SSM<Two_scale>* get_model() const;
   void initialize();
   void set_scale(double); ///< fix unification scale (0 = unfixed)

protected:
   void update_scale();
   bool check_non_perturbative();

private:
   double scale;
   double initial_scale_guess;
   E6SSM<Two_scale>* model;
};

} // namespace flexiblesusy

#endif
