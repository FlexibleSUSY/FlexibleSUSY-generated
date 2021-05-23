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

LIBMSSMEFTHiggs_SRC += \
		$(DIR)/decays/MSSMEFTHiggs_decay_table.cpp \
		$(DIR)/decays/MSSMEFTHiggs_decays.cpp

LIBMSSMEFTHiggs_HDR += \
		$(DIR)/decays/MSSMEFTHiggs_decay_table.hpp \
		$(DIR)/decays/MSSMEFTHiggs_decays.hpp \
		$(DIR)/decays/MSSMEFTHiggs_decay_amplitudes.hpp

EXEMSSMEFTHiggs_SRC += \
		$(DIR)/run_decays_MSSMEFTHiggs.cpp
