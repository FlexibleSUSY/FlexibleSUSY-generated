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

E6SSMEFTHiggs_INCLUDE_MK += $(DIR)/two_scale.mk

LIBE6SSMEFTHiggs_SRC += \
		$(DIR)/E6SSMEFTHiggs_two_scale_convergence_tester.cpp \
		$(DIR)/E6SSMEFTHiggs_two_scale_ewsb_solver.cpp \
		$(DIR)/E6SSMEFTHiggs_two_scale_high_scale_constraint.cpp \
		$(DIR)/E6SSMEFTHiggs_two_scale_initial_guesser.cpp \
		$(DIR)/E6SSMEFTHiggs_two_scale_low_scale_constraint.cpp \
		$(DIR)/E6SSMEFTHiggs_two_scale_model.cpp \
		$(DIR)/E6SSMEFTHiggs_two_scale_spectrum_generator.cpp \
		$(DIR)/E6SSMEFTHiggs_two_scale_susy_scale_constraint.cpp
LIBE6SSMEFTHiggs_HDR += \
		$(DIR)/E6SSMEFTHiggs_two_scale_convergence_tester.hpp \
		$(DIR)/E6SSMEFTHiggs_two_scale_ewsb_solver.hpp \
		$(DIR)/E6SSMEFTHiggs_two_scale_high_scale_constraint.hpp \
		$(DIR)/E6SSMEFTHiggs_two_scale_initial_guesser.hpp \
		$(DIR)/E6SSMEFTHiggs_two_scale_low_scale_constraint.hpp \
		$(DIR)/E6SSMEFTHiggs_two_scale_model.hpp \
		$(DIR)/E6SSMEFTHiggs_two_scale_spectrum_generator.hpp \
		$(DIR)/E6SSMEFTHiggs_two_scale_susy_scale_constraint.hpp

LIBE6SSMEFTHiggs_SRC += \
		models/E6SSMEFTHiggs/E6SSMEFTHiggs_standard_model_two_scale_matching.cpp

LIBE6SSMEFTHiggs_HDR += \
		models/E6SSMEFTHiggs/E6SSMEFTHiggs_standard_model_two_scale_matching.hpp
