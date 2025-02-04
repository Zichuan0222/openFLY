// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/neigh/list.hpp"

#include <Eigen/src/Core/util/Meta.h>
#include <fmt/core.h>
#include <omp.h>

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <limits>
#include <random>

#include "libfly/neigh/sort.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/random.hpp"

using namespace fly;

struct Neigh {
  Eigen::Index i;
  fly::Vec dr;
};

void slow_neigh_list(Vector<Vector<Neigh>>& nl, system::Box const& box, system::SoA<Position const&> atoms, double r_cut) {
  //
  nl.resize(atoms.size());

  // The displacement vectors move +/- 1 or zero supercells in each direction.
  std::vector<Vec> disp_vectors;

  Mat basis = box.basis();

  template_for<int>(Arr<int>::Constant(-1), Arr<int>::Constant(2), [&](auto... args) {
    //
    Arr<int> signs{args...};

    Vec offset = Vec::Zero();

    for (int i = 0; i < spatial_dims; i++) {
      if (box.periodic(i)) {  // Only in periodic axis
        offset += basis.col(i) * signs(i);
      }
    }

    if ((signs == 0).all() || offset != Vec::Zero()) {
      disp_vectors.push_back(offset);
    }
  });

  for (Eigen::Index i = 0; i < atoms.size(); i++) {
    //
    nl[i].clear();

    for (Eigen::Index j = 0; j < atoms.size(); j++) {
      if (i != j) {
        Neigh best{
            i,
            Vec::Constant(std::numeric_limits<double>::max()),
        };
        // Brute for min-image
        for (auto const& d : disp_vectors) {
          Vec dr = (atoms(r_, j) + d) - atoms(r_, i);

          if (gnorm(dr) < gnorm(best.dr)) {
            best.i = j;
            best.dr = dr;
          }
        }

        ASSERT(best.i != i, "should not have found self", 0);

        if (gnorm_sq(best.dr) < r_cut * r_cut) {
          nl[i].push_back(best);
        }
      }
    }
    std::sort(nl[i].begin(), nl[i].end(), [](auto const& a, auto const& b) { return a.i < b.i; });
  }
}

void test(Vector<Vector<Neigh>> const& nl, neigh::List const& neigh, double r_cut) {
  //
  //   neigh::List neigh(box, r_cut);

  //   timeit(fmt::format("\trebuild {:<2}", num_threads), [&] { neigh.rebuild(atoms, num_threads); });

  for (Eigen::Index i = 0; i < neigh.size(); i++) {
    //
    Vector<Neigh> nl2;

    neigh.for_neighbours(i, r_cut, [&](Eigen::Index n, double, Vec const& dr) {
      //
      nl2.push_back({n, dr});
    });

    std::sort(nl2.begin(), nl2.end(), [](auto const& a, auto const& b) { return a.i < b.i; });

    // // Test same number of neighbours
    REQUIRE(nl2.size() == nl[i].size());

    for (Eigen::Index j = 0; j < nl2.size(); j++) {
      // Test all neighbours have the same index
      REQUIRE(nl2[j].i == nl[i][j].i);
      // Test all neighbours have the same minimum image positions
      REQUIRE(near(0., gnorm(nl2[j].dr - nl[i][j].dr)));
    }
  }
}

std::pair<system::SoA<Position>, system::Box> gen_rand(fly::Xoshiro& gen, bool full_periodic = false) {
  //

  std::uniform_real_distribution<double> dis(0, 1);
  std::uniform_int_distribution<Eigen::Index> idis(100, 1000);

  auto rand = [&dis, &gen]() { return dis(gen); };

  // Random simulation box
  auto box = [&]() -> fly::system::Box {
    //
    Mat basis = Mat::NullaryExpr(rand);

    basis.diagonal() += Vec::Ones() + 9 * Vec::NullaryExpr(rand);

    basis = basis.triangularView<Eigen::Upper>();

    if (dis(gen) < 0.5) {
      basis = basis.triangularView<Eigen::Lower>();
    }

    if (full_periodic) {
      return {basis, fly::Arr<bool>::Constant(true)};
    } else {
      return {basis, fly::Arr<double>::NullaryExpr(rand) < 0.5};
    }
  }();

  system::SoA<Position> cell(idis(gen));

  Mat basis = box.basis();

  for (int j = 0; j < cell.size(); j++) {
    cell(r_, j).noalias() = basis * fly::Vec::NullaryExpr(rand);
  }

  return {cell, box};
}

TEST_CASE("List fuzz-testing", "[neigh]") {
  //

  //
  fly::Xoshiro gen({1, 3, 3, 4});

  for (int i = 0; i < 100; i++) {
    //
    auto [cell, box] = gen_rand(gen);

    double r_cut_max = visit(box.get(), [](auto const& g) { return g.min_width(); });

    double r_cut = r_cut_max / 3;

    cell = neigh::sort(box, r_cut, cell);

    static Vector<Vector<Neigh>> nl;

    Eigen::Index sum = 0;

    for (auto&& elem : nl) {
      sum += elem.size();
    }

    fmt::print("{} atoms and shaped={}, avg_num_neigh={}\n",
               cell.size(),
               visit(box.make_grid(r_cut), [](auto const g) { return g.shape(); }) - 2,
               static_cast<double>(sum) / static_cast<double>(nl.size()));

    timeit("\tslow build", [cell = cell, box = box, r_cut] { slow_neigh_list(nl, box, cell, r_cut); });

    neigh::List neigh(box, r_cut);

    timeit(fmt::format("\trebuild {:<2}", 1), [&, cell = cell] { neigh.rebuild(cell, 1); });

    test(nl, neigh, r_cut);

    timeit(fmt::format("\trebuild {:<2}", omp_get_max_threads()), [&, cell = cell] { neigh.rebuild(cell, omp_get_max_threads()); });

    test(nl, neigh, r_cut);
  }
}

TEST_CASE("List::update()", "[neigh]") {
  //
  fly::Xoshiro gen({1, 1, 1, 1});

  for (int anon = 0; anon < 10; anon++) {
    //
    auto [cell, box] = gen_rand(gen, true);

    double r_cut_max = visit(box.get(), [](auto const& g) { return g.min_width(); });

    double r_cut = r_cut_max / 3;

    double skin = 0.1;

    neigh::List neigh(box, r_cut + 2 * skin);

    //////////////////////////

    neigh.rebuild(cell, omp_get_max_threads());

    //////////////////////////

    system::SoA<Delta> deltas(cell.size());

    std::uniform_real_distribution<double> dis(0, 1);

    // Set up a displacement less than the skin distance.
    for (int i = 0; i < deltas.size(); i++) {
      deltas(del_, i) = Vec::NullaryExpr([&] { return dis(gen); }).normalized() * (skin * 0.99);
    }

    neigh.update(deltas);

    Vector<Vector<Neigh>> nl;

    cell[r_] -= deltas[del_];

    slow_neigh_list(nl, box, cell, r_cut);

    test(nl, neigh, r_cut);
  }
}
