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

MRSSMEFTHiggs_INCLUDE_MK += $(DIR)/two_scale.mk

LIBMRSSMEFTHiggs_SRC += \
		$(DIR)/MRSSMEFTHiggs_two_scale_convergence_tester.cpp \
		$(DIR)/MRSSMEFTHiggs_two_scale_ewsb_solver.cpp \
		$(DIR)/MRSSMEFTHiggs_two_scale_high_scale_constraint.cpp \
		$(DIR)/MRSSMEFTHiggs_two_scale_initial_guesser.cpp \
		$(DIR)/MRSSMEFTHiggs_two_scale_low_scale_constraint.cpp \
		$(DIR)/MRSSMEFTHiggs_two_scale_model.cpp \
		$(DIR)/MRSSMEFTHiggs_two_scale_spectrum_generator.cpp \
		$(DIR)/MRSSMEFTHiggs_two_scale_susy_scale_constraint.cpp
LIBMRSSMEFTHiggs_HDR += \
		$(DIR)/MRSSMEFTHiggs_two_scale_convergence_tester.hpp \
		$(DIR)/MRSSMEFTHiggs_two_scale_ewsb_solver.hpp \
		$(DIR)/MRSSMEFTHiggs_two_scale_high_scale_constraint.hpp \
		$(DIR)/MRSSMEFTHiggs_two_scale_initial_guesser.hpp \
		$(DIR)/MRSSMEFTHiggs_two_scale_low_scale_constraint.hpp \
		$(DIR)/MRSSMEFTHiggs_two_scale_model.hpp \
		$(DIR)/MRSSMEFTHiggs_two_scale_spectrum_generator.hpp \
		$(DIR)/MRSSMEFTHiggs_two_scale_susy_scale_constraint.hpp

LIBMRSSMEFTHiggs_SRC += \
		models/MRSSMEFTHiggs/MRSSMEFTHiggs_standard_model_two_scale_matching.cpp

LIBMRSSMEFTHiggs_HDR += \
		models/MRSSMEFTHiggs/MRSSMEFTHiggs_standard_model_two_scale_matching.hpp
