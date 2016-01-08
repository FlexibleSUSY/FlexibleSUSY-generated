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

// File generated at Fri 8 Jan 2016 12:30:22

#ifndef UMSSM_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H
#define UMSSM_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H

#include "UMSSM_susy_scale_constraint.hpp"
#include "UMSSM_input_parameters.hpp"
#include "two_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class UMSSM;

class Two_scale;

template<>
class UMSSM_susy_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   UMSSM_susy_scale_constraint();
   UMSSM_susy_scale_constraint(UMSSM<Two_scale>*);
   virtual ~UMSSM_susy_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   double get_initial_scale_guess() const;
   const UMSSM_input_parameters& get_input_parameters() const;
   UMSSM<Two_scale>* get_model() const;
   void initialize();

protected:
   void update_scale();

private:
   double scale;
   double initial_scale_guess;
   UMSSM<Two_scale>* model;
};

} // namespace flexiblesusy

#endif
