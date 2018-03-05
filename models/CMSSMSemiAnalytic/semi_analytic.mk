#  ====================================================================
#  This file is part of FlexibleSUSY.
#
#  FlexibleSUSY is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation, either version 3 of the License,
#  or (at your option) any later version.
#
#  FlexibleSUSY is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FlexibleSUSY.  If not, see
#  <http://www.gnu.org/licenses/>.
#  ====================================================================

CMSSMSemiAnalytic_INCLUDE_MK += $(DIR)/semi_analytic.mk

LIBCMSSMSemiAnalytic_SRC += \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_convergence_tester.cpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_ewsb_solver.cpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_high_scale_constraint.cpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_initial_guesser.cpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_low_scale_constraint.cpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_model.cpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_soft_parameters_constraint.cpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_solutions.cpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_spectrum_generator.cpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_susy_convergence_tester.cpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_susy_scale_constraint.cpp
LIBCMSSMSemiAnalytic_HDR += \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_convergence_tester.hpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_ewsb_solver.hpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_high_scale_constraint.hpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_initial_guesser.hpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_low_scale_constraint.hpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_model.hpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_soft_parameters_constraint.hpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_solutions.hpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_spectrum_generator.hpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_susy_convergence_tester.hpp \
		$(DIR)/CMSSMSemiAnalytic_semi_analytic_susy_scale_constraint.hpp \
		$(DIR)/CMSSMSemiAnalytic_soft_parameters_constraint.hpp \
		$(DIR)/CMSSMSemiAnalytic_susy_convergence_tester.hpp
