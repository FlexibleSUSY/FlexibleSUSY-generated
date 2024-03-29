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


#ifndef MRSSM2_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H
#define MRSSM2_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H

#include "MRSSM2_susy_scale_constraint.hpp"
#include "MRSSM2_input_parameters.hpp"
#include "single_scale_constraint.hpp"
#include "lowe.h"

namespace flexiblesusy {

template <class T>
class MRSSM2;

class Two_scale;

template<>
class MRSSM2_susy_scale_constraint<Two_scale> : public Single_scale_constraint {
public:
   MRSSM2_susy_scale_constraint() = default;
   MRSSM2_susy_scale_constraint(MRSSM2<Two_scale>*, const softsusy::QedQcd&);
   virtual ~MRSSM2_susy_scale_constraint() = default;

   virtual void apply() override;
   virtual double get_scale() const override;
   virtual std::string name() const override { return "MRSSM2 SUSY-scale constraint"; }
   virtual void set_model(Model*) override;

   void clear();
   double get_initial_scale_guess() const;
   const MRSSM2_input_parameters& get_input_parameters() const;
   MRSSM2<Two_scale>* get_model() const;
   void initialize();
   const softsusy::QedQcd& get_sm_parameters() const;
   void set_sm_parameters(const softsusy::QedQcd&);

protected:
   void update_scale();

private:
   double scale{0.};
   double initial_scale_guess{0.};
   MRSSM2<Two_scale>* model{nullptr};
   softsusy::QedQcd qedqcd{};

   void check_model_ptr() const;
};

} // namespace flexiblesusy

#endif
