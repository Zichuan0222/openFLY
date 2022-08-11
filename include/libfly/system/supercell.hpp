#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <cstddef>
#include <type_traits>

#include "libfly/system/SoA.hpp"
#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/utility/core.hpp"
/**
 * \file supercell.hpp
 *
 * @brief Classes to represent entire systems of atoms.
 */

namespace fly::system {

  template <typename...>
  class TypeMap;

  namespace detail {

    template <typename...>
    struct same_properties : std::false_type {};

    template <typename... T>
    struct same_properties<TypeMap<T...>, T...> : std::true_type {};

  }  // namespace detail

  /**
   * @brief Utility for mapping an atoms TypeID to atomic properties.
   *
   * This is used for properties that are the same for many atoms e.g. atomic numbers.
   *
   * @tparam Mems Tags derived from ``MemTag``, to describe each property.
   */
  template <typename... Mems>
  class TypeMap : private SoA<Type, Mems...> {
  private:
    using SOA = SoA<Type, Mems...>;

    static_assert(SOA::owns_all, "TypeMap must own all its data,");

  public:
    /**
     * @brief Construct a TypeMap to hold ``num_types`` types.
     *
     * \rst
     * .. warning::
     *    All the new types properties are uninitialised.
     * \endrst
     */
    explicit TypeMap(Eigen::Index num_types) : SOA(num_types) {}

    /**
     * @brief Construct a new TypeMap by slicing a different kind of TypeMap.
     */
    template <typename... T, typename = std::enable_if_t<!detail::same_properties<TypeMap, T...>::value>>
    explicit TypeMap(TypeMap<T...> map) : SOA(static_cast<SoA<Type, T...>&&>(map)) {}

    /**
     * @brief Construct a new Type Map object
     */
    TypeMap(TypeMap const&) = default;

    /**
     * @brief Construct a new Type Map object
     */
    TypeMap(TypeMap&&) = default;

    /**
     * @brief Fetch the number of types stored in the TypeMap
     */
    Eigen::Index num_types() const { return SOA::size(); }

    /**
     * @brief Get a property corresponding to the id ``id``.
     */
    template <typename T>
    decltype(auto) get(T, std::uint32_t id) const {
      return SOA::operator()(T{}, safe_cast<Eigen::Index>(id));
    }

    /**
     * @brief Set a property corresponding to the id ``id``.
     */
    template <typename T, typename U = typename remove_cref_t<T>::matrix_t>
    void set(T, std::uint32_t id, U&& value) {
      SOA::operator()(T{}, safe_cast<Eigen::Index>(id)) = std::forward<U>(value);
    }

  private:
    /* clang-format off */ template <typename...>  friend class TypeMap; /* clang-format on */
  };

  /**
   * @brief LibFLY's representation of a system of atoms
   *
   * The Supercell is libFLY's amalgamation of all the data required for a simulation. It **is a** SoA containing all the atoms in the
   * system **has a** ``Box`` and a ``TypeMap``. A Supercell always has ``Position`` and ``TypeID`` members for each atom and accepts
   * the rest as template arguments.
   *
   * @tparam Mems Tags derived from ``MemTag``, to describe each member.
   */
  template <typename Map = TypeMap<>, typename... Mems>
  class Supercell : public SoA<TypeID, Position, Mems...> {
  private:
    using SOA = SoA<TypeID, Position, Mems...>;

    static_assert(SOA::owns_all, "Supercells must own all their data");

  public:
    /**
     * @brief Construct a new Supercell object to store `num_atoms` atoms.
     */
    explicit Supercell(Box const& box, Map const& map, Eigen::Index num_atoms) : SOA(num_atoms), m_map(map), m_box(box) {}

    /**
     * @brief Simulation box getter.
     */
    Box const& box() const { return m_box; }

    /**
     * @brief Simulation box getter.
     */
    Box& box() { return m_box; }

    /**
     * @brief TypeMap getter, read-only.
     */
    Map const& map() const { return m_map; }

  private:
    Box m_box;
    Map m_map;
  };

}  // namespace fly::system