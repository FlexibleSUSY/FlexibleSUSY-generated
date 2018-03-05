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

// File generated at Mon 5 Mar 2018 18:57:18

#ifndef UMSSM_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H
#define UMSSM_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H

#include "UMSSM_susy_scale_constraint.hpp"
#include "UMSSM_input_parameters.hpp"
#include "single_scale_constraint.hpp"
#include "lowe.h"

namespace flexiblesusy {

template <class T>
class UMSSM;

class Two_scale;

template<>
class UMSSM_susy_scale_constraint<Two_scale> : public Single_scale_constraint {
public:
   UMSSM_susy_scale_constraint() = default;
   UMSSM_susy_scale_constraint(UMSSM<Two_scale>*, const softsusy::QedQcd&);
   virtual ~UMSSM_susy_scale_constraint() = default;
   virtual void apply() override;
   virtual double get_scale() const override;
   virtual std::string name() const override { return "UMSSM SUSY-scale constraint"; }
   virtual void set_model(Model*) override;

   void clear();
   double get_initial_scale_guess() const;
   const UMSSM_input_parameters& get_input_parameters() const;
   UMSSM<Two_scale>* get_model() const;
   void initialize();
   const softsusy::QedQcd& get_sm_parameters() const;
   void set_sm_parameters(const softsusy::QedQcd&);

protected:
   void update_scale();

private:
   double scale{0.};
   double initial_scale_guess{0.};
   UMSSM<Two_scale>* model{nullptr};
   softsusy::QedQcd qedqcd{};

   void check_model_ptr() const;
};

} // namespace flexiblesusy

#endif
