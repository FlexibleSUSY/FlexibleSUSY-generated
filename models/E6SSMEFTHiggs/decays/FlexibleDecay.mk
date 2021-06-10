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

LIBE6SSMEFTHiggs_SRC += \
		$(DIR)/decays/E6SSMEFTHiggs_decay_table.cpp \
		$(DIR)/decays/E6SSMEFTHiggs_decays.cpp

LIBE6SSMEFTHiggs_HDR += \
		$(DIR)/decays/E6SSMEFTHiggs_decay_table.hpp \
		$(DIR)/decays/E6SSMEFTHiggs_decays.hpp \
		$(DIR)/decays/E6SSMEFTHiggs_decay_amplitudes.hpp

EXEE6SSMEFTHiggs_SRC += \
		$(DIR)/run_decays_E6SSMEFTHiggs.cpp
