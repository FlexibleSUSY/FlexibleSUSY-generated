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

NMSSMEFTHiggs_INCLUDE_MK += $(DIR)/two_scale.mk

LIBNMSSMEFTHiggs_SRC += \
		$(DIR)/NMSSMEFTHiggs_two_scale_convergence_tester.cpp \
		$(DIR)/NMSSMEFTHiggs_two_scale_ewsb_solver.cpp \
		$(DIR)/NMSSMEFTHiggs_two_scale_high_scale_constraint.cpp \
		$(DIR)/NMSSMEFTHiggs_two_scale_initial_guesser.cpp \
		$(DIR)/NMSSMEFTHiggs_two_scale_low_scale_constraint.cpp \
		$(DIR)/NMSSMEFTHiggs_two_scale_model.cpp \
		$(DIR)/NMSSMEFTHiggs_two_scale_spectrum_generator.cpp \
		$(DIR)/NMSSMEFTHiggs_two_scale_susy_scale_constraint.cpp
LIBNMSSMEFTHiggs_HDR += \
		$(DIR)/NMSSMEFTHiggs_two_scale_convergence_tester.hpp \
		$(DIR)/NMSSMEFTHiggs_two_scale_ewsb_solver.hpp \
		$(DIR)/NMSSMEFTHiggs_two_scale_high_scale_constraint.hpp \
		$(DIR)/NMSSMEFTHiggs_two_scale_initial_guesser.hpp \
		$(DIR)/NMSSMEFTHiggs_two_scale_low_scale_constraint.hpp \
		$(DIR)/NMSSMEFTHiggs_two_scale_model.hpp \
		$(DIR)/NMSSMEFTHiggs_two_scale_spectrum_generator.hpp \
		$(DIR)/NMSSMEFTHiggs_two_scale_susy_scale_constraint.hpp

LIBNMSSMEFTHiggs_SRC += \
		models/NMSSMEFTHiggs/NMSSMEFTHiggs_standard_model_two_scale_matching.cpp

LIBNMSSMEFTHiggs_HDR += \
		models/NMSSMEFTHiggs/NMSSMEFTHiggs_standard_model_two_scale_matching.hpp
