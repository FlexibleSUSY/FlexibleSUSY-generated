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

E6SSMEFTHiggs_INCLUDE_MK += $(DIR)/FlexibleEFTHiggs.mk

LIBE6SSMEFTHiggs_SRC += \
		models/E6SSMEFTHiggs/E6SSMEFTHiggs_standard_model_matching.cpp

LIBE6SSMEFTHiggs_HDR += \
		models/E6SSMEFTHiggs/E6SSMEFTHiggs_standard_model_matching.hpp
