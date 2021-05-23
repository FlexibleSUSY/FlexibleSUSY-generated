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


#include "MSSMRHN_decay_table.hpp"

#include "error.hpp"

#include <iostream>
#include <sstream>

namespace flexiblesusy {

MSSMRHN_decay_table::MSSMRHN_decay_table()
   : decay_table({Decays_list(25), Decays_list(35), Decays_list(-37), Decays_list(
      36)})
{
}

void MSSMRHN_decay_table::clear()
{
   for (auto& particle: decay_table) {
      particle.clear();
   }
}

void MSSMRHN_decay_table::print(std::ostream& ostr) const
{
   for (const auto& particle: decay_table) {
      const auto pdg = particle.get_particle_id();
      for (const auto& decay: particle) {
         const auto& final_state = decay.second.get_final_state_particle_ids();
         const std::size_t final_state_size = final_state.size();
         ostr << "Partial width["
              << MSSMRHN_info::get_particle_name_from_pdg(pdg)
              << " -> ";
         for (std::size_t i = 0; i < final_state_size; ++i) {
            ostr << MSSMRHN_info::get_particle_name_from_pdg(final_state[i]);
            if (i < final_state_size - 1) {
               ostr << ", ";
            }
         }
         ostr << "] = " << decay.second.get_width() << '\n';
      }
   }
}

std::ostream& operator<<(std::ostream& ostr,
                         const MSSMRHN_decay_table& table)
{
   table.print(ostr);
   return ostr;
}

Decays_list& MSSMRHN_decay_table::get_hh_decays(int gI1)
{
   switch (gI1) {
   case 0: return decay_table[0]; break;
   case 1: return decay_table[1]; break;
   }

   std::ostringstream sstr;
   sstr << "invalid particle index " << std::to_string(gI1) << '\n';

   throw OutOfBoundsError(sstr.str());
}

const Decays_list& MSSMRHN_decay_table::get_hh_decays(int gI1) const
{
   switch (gI1) {
   case 0: return decay_table[0]; break;
   case 1: return decay_table[1]; break;
   }

   std::ostringstream sstr;
   sstr << "invalid particle index " << std::to_string(gI1) << '\n';

   throw OutOfBoundsError(sstr.str());
}
Decays_list& MSSMRHN_decay_table::get_Hpm_decays(int gI1)
{
   switch (gI1) {
   case 1: return decay_table[2]; break;
   }

   std::ostringstream sstr;
   sstr << "invalid particle index " << std::to_string(gI1) << '\n';

   throw OutOfBoundsError(sstr.str());
}

const Decays_list& MSSMRHN_decay_table::get_Hpm_decays(int gI1) const
{
   switch (gI1) {
   case 1: return decay_table[2]; break;
   }

   std::ostringstream sstr;
   sstr << "invalid particle index " << std::to_string(gI1) << '\n';

   throw OutOfBoundsError(sstr.str());
}
Decays_list& MSSMRHN_decay_table::get_Ah_decays(int gI1)
{
   switch (gI1) {
   case 1: return decay_table[3]; break;
   }

   std::ostringstream sstr;
   sstr << "invalid particle index " << std::to_string(gI1) << '\n';

   throw OutOfBoundsError(sstr.str());
}

const Decays_list& MSSMRHN_decay_table::get_Ah_decays(int gI1) const
{
   switch (gI1) {
   case 1: return decay_table[3]; break;
   }

   std::ostringstream sstr;
   sstr << "invalid particle index " << std::to_string(gI1) << '\n';

   throw OutOfBoundsError(sstr.str());
}

} // namespace flexiblesusy
