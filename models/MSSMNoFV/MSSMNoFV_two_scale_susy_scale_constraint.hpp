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

// File generated at Tue 8 Sep 2015 13:58:55

#ifndef MSSMNoFV_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H
#define MSSMNoFV_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H

#include "MSSMNoFV_susy_scale_constraint.hpp"
#include "MSSMNoFV_input_parameters.hpp"
#include "two_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class MSSMNoFV;

class Two_scale;

template<>
class MSSMNoFV_susy_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   MSSMNoFV_susy_scale_constraint();
   MSSMNoFV_susy_scale_constraint(MSSMNoFV<Two_scale>*);
   virtual ~MSSMNoFV_susy_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   double get_initial_scale_guess() const;
   const MSSMNoFV_input_parameters& get_input_parameters() const;
   MSSMNoFV<Two_scale>* get_model() const;
   void initialize();

protected:
   void update_scale();

private:
   double scale;
   double initial_scale_guess;
   MSSMNoFV<Two_scale>* model;
};

} // namespace flexiblesusy

#endif
