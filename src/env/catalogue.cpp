// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/env/catalogue.hpp"

#include <fmt/core.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "libfly/env/geometry.hpp"
#include "libfly/env/heuristics.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

namespace fly::env {

  void Catalogue::rebuild_env(int i,
                              Eigen::Index num_types,
                              Geometry<Index> &scratch,
                              system::SoA<Position const &, TypeID const &, Frozen const &> const &info) {
    //
    RelEnv &env = m_real[std::size_t(i)];
    // Build geometries
    rebuild_geo_from_nl(i, env.geo, info, *m_nl, m_opt.r_env);
    // Make fingerprint.
    env.f.rebuild(env.geo);
    // Compute graph-hash and canonise
    env.hash = canon_hash(env.geo, m_opt.r_edge, safe_cast<std::size_t>(2 * num_types), &scratch);
  }

  std::vector<int> Catalogue::rebuild_impl(system::SoA<Position const &, TypeID const &, Frozen const &> const &info,
                                           Eigen::Index num_types,
                                           int num_threads) {
    // Prepare memory.
    m_real.resize(std::size_t(info.size()));

    // Prep neigh list.
    m_nl->rebuild(info, num_threads);

    Geometry<Index> scratch;
    // Rebuild the RelEnv's
#pragma omp parallel for num_threads(num_threads) firstprivate(scratch) schedule(static)
    for (int i = 0; i < info.size(); i++) {
      rebuild_env(i, num_types, scratch, info);
    }

    // Make sure every hash is in the map.
    for (auto const &elem : m_real) {
      m_cat.try_emplace(elem.hash);
    }

    bool flag = false;

    // Find in catalogue, optimising for found case
#pragma omp parallel for num_threads(num_threads) schedule(static)
    for (std::size_t i = 0; i < m_real.size(); i++) {
      //   if (!(m_real[i].ptr = timeit("canon", [&] { return canon_find(m_real[i]); }))) {
      if (!(m_real[i].ptr = canon_find(m_real[i]))) {
#pragma omp atomic write
        flag = true;
        if (m_opt.debug) {
          fmt::print("CAT: environment around {} is new\n", i);
        }
      }
    }

    if (!flag) {
      if (m_opt.debug) {
        fmt::print("CAT: No new environments\n");
      }
      return {};
    }

    // Now must operate single threaded as we modify buckets
    std::vector<int> new_idx;

    struct Pair {
      Pointer ptr;
      std::size_t hash;
    };

    std::vector<Pair> new_ptrs;

    // Insert missing
    for (std::size_t i = 0; i < m_real.size(); i++) {
      if (!m_real[i].ptr) {
        //
        auto it = std::find_if(new_ptrs.begin(), new_ptrs.end(), [&](Pair const &ref) {
          return m_real[i].hash == ref.hash && canon_equiv(m_real[i], *(ref.ptr));
          //
        });

        if (it == new_ptrs.end()) {
          //
          m_real[i].ptr = insert(m_real[i]);
          new_idx.push_back(int(i));
          new_ptrs.push_back({*m_real[i].ptr, m_real[i].hash});

          if (m_opt.debug) {
            fmt::print("CAT: Unknown at {} is new\n", i);
          }

        } else {
          // New environment equivalent to other new environment.
          m_real[i].ptr = it->ptr;

          if (m_opt.debug) {
            fmt::print("CAT: Unknown at {} is a duplicate new\n", i);
          }
        }
      }
    }

    if (m_opt.debug) {
      fmt::print("CAT: found {} new environments\n", new_idx.size());
    }

    // Update frequencies
    for (auto const &elem : m_real) {
      ASSERT(elem.ptr, "Null pointer in catalogue\n", 0);
      (**elem.ptr).m_freq++;
    }

    return new_idx;
  }

  std::optional<Catalogue::Pointer> Catalogue::canon_find(Catalogue::RelEnv &env) {
    // "it" always points to valid bucket (possibly empty)
    auto it = m_cat.find(env.hash);

    ASSERT(it != m_cat.end(), "Catalogue missing key!", 0);

    // Existing key, must search bucket for explicit match
    auto match = std::find_if(it->second.begin(), it->second.end(), [&](Env const &ref) { return canon_equiv(env, ref); });

    // If found a match, return it
    if (match != it->second.end()) {
      return Pointer(it, match - it->second.begin());
    } else {
      return {};
    }
  }

  Catalogue::Pointer Catalogue::insert(Catalogue::RelEnv &env) {
    //

    ASSERT(!canon_find(env), "Environment already in catalogue.", 0);

    auto it = m_cat.find(env.hash);

    ASSERT(it != m_cat.end(), "Catalogue missing key!", 0);

    it->second.emplace_back(Env{env.geo, env.f, m_size++, std::min(env.f.r_min() * 0.4, m_opt.delta_max)});

    return {it, xise(it->second) - 1};
  }

  bool Catalogue::canon_equiv(Catalogue::RelEnv &mut, Env const &ref) const {
    //
    double delta = calc_delta(mut.f, ref);

    // Test if fuzzy keys match (fast)
    if (!ref.m_finger.equiv(mut.f, delta * m_opt.overfuzz)) {
      return false;
    }

    if (mut.geo.permute_onto(ref, delta)) {
      return true;
    } else {
      //   fmt::print("False positive\n");
      return false;
    }
  }

  auto Catalogue::set_mechs(int i, std::vector<Mechanism> const &m) -> void {
    //
    Env &env = **(m_real[std::size_t(i)].ptr);

    verify(env.m_mechs.empty(), "We already have {} mechanisms, set_mech() should only be called once.", env.m_mechs.size());

    verify(std::all_of(m.begin(),
                       m.end(),
                       [s = env.size()](Mechanism const &x) {
                         //
                         return x.delta_sp.size() == s && x.delta_fwd.size() == s;
                       }),
           "Wrong number of atoms.");

    env.m_mechs = m;
  }

  double Catalogue::refine_tol(int i, double min_delta) {
    //
    std::size_t si = safe_cast<std::size_t>(i);

    double delta = calc_delta(m_real[si].f, get_ref(i));

    std::optional res = m_real[si].geo.permute_onto(get_ref(i).ref_geo(), delta);

    verify(bool(res), "While tightening @{} perm failed with delta={}", i, delta);

    double new_delta_max = std::max(min_delta, res->rmsd / 1.5);

    if (m_opt.debug) {
      fmt::print("Refining delta_max @{} from {} to {}\n", i, get_ref(i).m_delta_max, new_delta_max);
    }

    (**(m_real[si].ptr)).m_delta_max = new_delta_max;

    return new_delta_max;
  }

  void Catalogue::reconstruct_impl(Mechanism const &mech,
                                   int i,
                                   system::SoA<Position const &, TypeID const &, Frozen const &> in,
                                   system::SoA<Position &> out,
                                   bool in_ready_state,
                                   Eigen::Index num_types,
                                   int num_threads) {
    //

    verify(out.size() == in.size(), "Output system is a different size from input system {}!={}", out.size(), in.size());

    out[r_] = in[r_];

    std::size_t si = safe_cast<std::size_t>(i);

    if (!in_ready_state) {
      // Prepare memory.
      m_real.resize(std::size_t(in.size()));
      // Prep neigh list.
      m_nl->rebuild(in, num_threads);

      Geometry<Index> scratch;

      rebuild_env(i, num_types, scratch, in);

      m_real[si].ptr = canon_find(m_real[si]);

      if (!m_real[si].ptr) {
        throw error("Rebuilding the geo centred on {} cat failed to find a matching Env", i);
      }
    } else {
      verify(i >= 0 && size_t(i) < m_real.size(), "Atom with index {} is not in catalogue", i);
    }

    double delta = calc_delta(m_real[si].f, get_ref(i));

    if (std::optional res = m_real[si].geo.permute_onto(get_ref(i), delta)) {
      //
      auto const &geo = get_geo(i);

      Mat O = res->O.transpose();

      for (int j = 0; j < geo.size(); ++j) {
        out(r_, geo[j][i_]).noalias() += O * mech.delta_fwd[j][del_];
      }
    } else {
      throw error("Unable to align cell perm geo {} onto reference, ready={}", i, in_ready_state);
    }
  }

}  // namespace fly::env