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


#ifndef CNMSSM_SEMI_ANALYTIC_SUSY_SCALE_CONSTRAINT_H
#define CNMSSM_SEMI_ANALYTIC_SUSY_SCALE_CONSTRAINT_H

#include "CNMSSM_susy_scale_constraint.hpp"
#include "CNMSSM_input_parameters.hpp"
#include "single_scale_constraint.hpp"
#include "lowe.h"

namespace flexiblesusy {

template <class T>
class CNMSSM;

class Semi_analytic;

template<>
class CNMSSM_susy_scale_constraint<Semi_analytic> : public Single_scale_constraint {
public:
   using Scale_getter = std::function<double()>;

   CNMSSM_susy_scale_constraint() = default;
   CNMSSM_susy_scale_constraint(CNMSSM<Semi_analytic>*, const softsusy::QedQcd&);
   virtual ~CNMSSM_susy_scale_constraint() = default;
   virtual void apply() override;
   virtual double get_scale() const override;
   virtual std::string name() const override { return "CNMSSM SUSY-scale constraint"; }
   virtual void set_model(Model*) override;

   void clear();
   double get_initial_scale_guess() const { return initial_scale_guess; }
   const CNMSSM_input_parameters& get_input_parameters() const;
   CNMSSM<Semi_analytic>* get_model() const { return model; }
   void initialize();
   const softsusy::QedQcd& get_sm_parameters() const { return qedqcd; }
   void set_sm_parameters(const softsusy::QedQcd& qedqcd_) { qedqcd = qedqcd_; }
   void set_scale(const Scale_getter& s) { scale_getter = s; }
   void set_scale(Scale_getter&& s) { scale_getter = std::move(s); }

protected:
   void update_scale();

private:
   double scale{0.};
   double initial_scale_guess{0.};
   CNMSSM<Semi_analytic>* model{nullptr};
   softsusy::QedQcd qedqcd{};
   Scale_getter scale_getter{};

   void check_model_ptr() const;
};

} // namespace flexiblesusy

#endif
