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

// File generated at Wed 16 Oct 2019 22:56:03

#ifndef CMSSMSemiAnalytic_SEMI_ANALYTIC_HIGH_SCALE_CONSTRAINT_H
#define CMSSMSemiAnalytic_SEMI_ANALYTIC_HIGH_SCALE_CONSTRAINT_H

#include "CMSSMSemiAnalytic_high_scale_constraint.hpp"
#include "CMSSMSemiAnalytic_input_parameters.hpp"
#include "single_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class CMSSMSemiAnalytic;

class Semi_analytic;

template<>
class CMSSMSemiAnalytic_high_scale_constraint<Semi_analytic> : public Single_scale_constraint {
public:

   CMSSMSemiAnalytic_high_scale_constraint() = default;
   CMSSMSemiAnalytic_high_scale_constraint(CMSSMSemiAnalytic<Semi_analytic>*);
   virtual ~CMSSMSemiAnalytic_high_scale_constraint() = default;
   virtual void apply() override;
   virtual double get_scale() const override;
   virtual std::string name() const override { return "CMSSMSemiAnalytic high-scale constraint"; }
   virtual void set_model(Model*) override;

   void clear();
   double get_initial_scale_guess() const { return initial_scale_guess; }
   const CMSSMSemiAnalytic_input_parameters& get_input_parameters() const;
   CMSSMSemiAnalytic<Semi_analytic>* get_model() const { return model; }
   void initialize();
   void set_scale(double s) { scale = s; }

protected:
   void update_scale();
   bool check_non_perturbative();

private:
   double scale{0.};
   double initial_scale_guess{0.};
   CMSSMSemiAnalytic<Semi_analytic>* model{nullptr};

   void check_model_ptr() const;
};

} // namespace flexiblesusy

#endif
