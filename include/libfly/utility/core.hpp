#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-2.0

// This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <type_traits>
#include <utility>

/**
 * \file core.hpp
 *
 * @brief Miscellaneous utilities.
 */

/**
 * @brief Specify the number of spatial dimensions at compile time, defaults to 3.
 *
 * \rst
 *
 * .. _`configure FLY_SPATIAL_DIMS`:
 *
 * To customize, at the compilation :ref:`configuration <Compiling openFLY>` step append:
 *
 * .. code:: console
 *
 *     -DCMAKE_CXX_FLAGS='-DFLY_SPATIAL_DIMS=<number>'
 *
 * replacing ``<number>`` with the number of dimensions.
 *
 * \endrst
 */
#ifndef FLY_SPATIAL_DIMS
#  define FLY_SPATIAL_DIMS 3
#endif  // !FLY_SPATIAL_DIMS

namespace fly {

  /**
   * @brief The number of spatial dimensions.
   *
   * \rst
   *
   * Configurable using the :ref:`FLY_SPATIAL_DIMS <configure FLY_SPATIAL_DIMS>` macro.
   *
   * \endrst
   */
  inline constexpr int spatial_dims = FLY_SPATIAL_DIMS;

  static_assert(spatial_dims >= 2, "libFLY is not optimal for 1D simulations");

  /**
   * @brief The maximum atomic number that any atom can have.
   */
  inline constexpr int max_atomic_num = 111;

  /**
   * @brief Shorthand for creating an \c Eigen::Vector of length \c spatial_dims
   *
   * @tparam T The scalar type of the \c Eigen::Vector
   */
  template <typename T>
  using Vec = Eigen::Vector<T, spatial_dims>;

  /**
   * @brief Shorthand for creating an \c spatial_dims x 1 \c Eigen::Array.
   *
   * @tparam T The scalar type of the \c Eigen::Array
   */
  template <typename T>
  using Arr = Eigen::Array<T, spatial_dims, 1>;

  /**
   * @brief Strongly typed +/- sign.
   */
  enum class Sign : int {
    plus = 1,
    minus = -1,
  };

  /**
   * @brief Shorthand for creating an \c spatial_dims x \c spatial_dims \c Eigen::Matrix.
   *
   * @tparam T The scalar type of the \c Eigen::Matrix
   */
  template <typename T>
  using Mat = Eigen::Matrix<T, spatial_dims, spatial_dims>;

  namespace detail {

    template <typename T, typename...>
    struct First {
      using type = T;
    };

  }  // namespace detail

  /**
   * @brief Extracts the first type from a parameter pack.
   */
  template <typename... Ts>
  using first_t = typename detail::First<Ts...>::type;

  /**
   * @brief Non-deducable \c false for use in \c static_assert.
   *
   * \rst
   *
   * Example:
   *
   * .. include:: ../../examples/utility/static_assert_false.cpp
   *    :code:
   *
   * \endrst
   */
  template <typename...>
  inline constexpr bool always_false = false;

  /**
   * @brief Strip all reference and const qualifications from \c T.
   *
   * @tparam T The type to strip ref/const qualifications from.
   */
  template <typename T>
  using remove_cref_t = std::remove_const_t<std::remove_reference_t<T>>;

  /**
   * @brief Test if two floating point numbers are within 0.01% of each other.
   */
  template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
  constexpr bool near(T a, T b) {
    return std::abs(a - b) <= 0.0001 * std::max(std::abs(a), std::abs(b));
  }

  /**
   * @brief Generalised dot-product between two \c Eigen::Array objects.
   *
   * \rst
   *
   * Computes:
   *
   * .. math::
   *
   *    \sum_{i} \sum_{j} a_{ij} b_{ij}
   *
   * \endrst
   */
  template <typename E1, typename E2>
  auto gdot(Eigen::ArrayBase<E1> const& a, Eigen::ArrayBase<E2> const& b) {
    return (a * b).sum();
  }

  /**
   * @brief Generic squared Frobenius norm of an \c Eigen::Array.
   *
   * \rst
   *
   * Computes:
   *
   * .. math::
   *
   *    \sum_{i} \sum_{j} r^2_{ij}
   *
   * \endrst
   */
  template <typename E>
  auto norm_sq(E const& r) {
    return (r.array() * r.array()).sum();
  }

  /**
   * @brief Generic Frobenius norm of an \c Eigen::Array.
   *
   * \rst
   *
   * Computes:
   *
   * .. math::
   *
   *    \sqrt{\sum_{i} \sum_{j} r^2_{ij}}
   *
   * \endrst
   */
  template <typename E>
  auto norm(E&& expr) {
    return std::sqrt(norm_sq(std::forward<E>(expr)));
  }

  /**
   * @brief Compute integer powers of arithmetic types at compile time.
   *
   * Only SFINE if base is an arithmetic type.
   *
   * \rst
   *
   * Computes:
   *
   * .. math::
   *
   *    \text{x}^{\text{Exp}}
   *
   * \endrst
   */
  template <std::size_t Exp, typename T>
  constexpr std::enable_if_t<std::is_arithmetic_v<T>, T> ipow(T x) {
    if constexpr (Exp == 0) {
      return T(1);
    }
    if constexpr (Exp == 1) {
      return x;
    }
    if constexpr (Exp % 2 == 0) {
      return ipow<Exp / 2>(x) * ipow<Exp / 2>(x);
    } else {
      return ipow<Exp - 1>(x) * x;
    }
  }

  /**
   * @brief Basic implementation of a Golang like defer.
   *
   * Only SFINE if callable is noexcept.
   *
   * \rst
   *
   * Example:
   *
   * .. include:: ../../examples/utility/defer.cpp
   *    :code:
   *
   * \endrst
   */
  template <class F, typename = std::enable_if_t<std::is_nothrow_invocable_v<F&&>>>
  class [[nodiscard]] Defer {
  public:
    /**
     * @brief Construct a new Defer object
     *
     * @param f Forwarded into object and invoked by destructor.
     */
    Defer(F&& f) : m_f(std::forward<F>(f)) {}

    Defer(const Defer&) = delete;
    Defer(Defer&& other) = delete;
    Defer& operator=(const Defer&) = delete;
    Defer& operator=(Defer&&) = delete;

    ~Defer() noexcept { std::invoke(std::forward<F>(m_f)); }

  private:
    F m_f;
  };

  /**
   * @brief Forwarding deduction guide.
   */
  template <typename F>
  Defer(F&&) -> Defer<F>;

}  // namespace fly
