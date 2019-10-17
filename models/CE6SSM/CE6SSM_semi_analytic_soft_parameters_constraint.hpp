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

// File generated at Wed 16 Oct 2019 20:12:13

#ifndef CE6SSM_SEMI_ANALYTIC_SOFT_PARAMETERS_CONSTRAINT_H
#define CE6SSM_SEMI_ANALYTIC_SOFT_PARAMETERS_CONSTRAINT_H

#include "CE6SSM_input_parameters.hpp"
#include "CE6SSM_soft_parameters_constraint.hpp"
#include "single_scale_constraint.hpp"
#include "lowe.h"

namespace flexiblesusy {

template <class T>
class CE6SSM;

class Semi_analytic;

template<>
class CE6SSM_soft_parameters_constraint<Semi_analytic> : public Single_scale_constraint {
public:
   using Scale_getter = std::function<double()>;

   CE6SSM_soft_parameters_constraint() = default;
   CE6SSM_soft_parameters_constraint(CE6SSM<Semi_analytic>*, const softsusy::QedQcd&);
   virtual ~CE6SSM_soft_parameters_constraint() = default;

   virtual void apply() override;
   virtual double get_scale() const override { return scale; }
   virtual std::string name() const override { return "CE6SSM soft parameters constraint"; }
   virtual void set_model(Model*) override;

   void clear();
   double get_initial_scale_guess() const { return initial_scale_guess; }
   double get_boundary_scale() const;
   void set_boundary_scale(const Scale_getter&);
   void set_boundary_scale(Scale_getter&&);
   const CE6SSM_input_parameters& get_input_parameters() const;
   CE6SSM<Semi_analytic>* get_model() const { return model; }
   void initialize();
   const softsusy::QedQcd& get_sm_parameters() const { return qedqcd; }
   void set_sm_parameters(const softsusy::QedQcd& qedqcd_) { qedqcd = qedqcd_; }

protected:
   void update_scale();

private:
   double scale{0.};
   double initial_scale_guess{0.};
   CE6SSM<Semi_analytic>* model{nullptr};
   softsusy::QedQcd qedqcd{};
   Scale_getter boundary_scale_getter{};

   void check_model_ptr() const;
};

} // namespace flexiblesusy

#endif
