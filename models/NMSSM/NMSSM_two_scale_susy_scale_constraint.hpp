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

// File generated at Tue 8 Sep 2015 13:12:33

#ifndef NMSSM_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H
#define NMSSM_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H

#include "NMSSM_susy_scale_constraint.hpp"
#include "NMSSM_input_parameters.hpp"
#include "two_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class NMSSM;

class Two_scale;

template<>
class NMSSM_susy_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   NMSSM_susy_scale_constraint();
   NMSSM_susy_scale_constraint(NMSSM<Two_scale>*);
   virtual ~NMSSM_susy_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   double get_initial_scale_guess() const;
   const NMSSM_input_parameters& get_input_parameters() const;
   NMSSM<Two_scale>* get_model() const;
   void initialize();

protected:
   void update_scale();

private:
   double scale;
   double initial_scale_guess;
   NMSSM<Two_scale>* model;
};

} // namespace flexiblesusy

#endif
