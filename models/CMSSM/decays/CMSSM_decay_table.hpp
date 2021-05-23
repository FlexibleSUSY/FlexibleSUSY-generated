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


#ifndef CMSSM_DECAY_TABLE_H
#define CMSSM_DECAY_TABLE_H

#include "CMSSM_info.hpp"

#include "decays/decay.hpp"

#include <array>
#include <iosfwd>

namespace flexiblesusy {

class CMSSM_decay_table {
private:
   static const int number_of_decay_particles = 4;
   using Table_type = std::array<Decays_list, number_of_decay_particles>;
public:
   using iterator = Table_type::iterator;
   using const_iterator = Table_type::const_iterator;

   CMSSM_decay_table();
   ~CMSSM_decay_table() = default;
   CMSSM_decay_table(const CMSSM_decay_table&) = default;
   CMSSM_decay_table(CMSSM_decay_table&&) = default;
   CMSSM_decay_table& operator=(const CMSSM_decay_table&) = default;
   CMSSM_decay_table& operator=(CMSSM_decay_table&&) = default;

   iterator begin() noexcept { return decay_table.begin(); }
   const_iterator begin() const noexcept { return decay_table.begin(); }
   const_iterator cbegin() const noexcept { return decay_table.cbegin(); }
   iterator end() noexcept { return decay_table.end(); }
   const_iterator end() const noexcept { return decay_table.end(); }
   const_iterator cend() const noexcept { return decay_table.end(); }

   std::size_t size() const noexcept { return decay_table.size(); }

   void clear();
   void print(std::ostream&) const;

   Decays_list& get_hh_decays(int);
   const Decays_list& get_hh_decays(int) const;
   Decays_list& get_Hpm_decays(int);
   const Decays_list& get_Hpm_decays(int) const;
   Decays_list& get_Ah_decays(int);
   const Decays_list& get_Ah_decays(int) const;
private:
   Table_type decay_table;
};

std::ostream& operator<<(std::ostream&, const CMSSM_decay_table&);

} // namespace flexiblesusy

#endif
