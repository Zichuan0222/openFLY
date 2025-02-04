// Copyright © 2020-2022 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see https :
// // www.gnu.org/licenses/>.

#include "libfly/saddle/find.hpp"

#include <fmt/core.h>
#include <fmt/format.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <optional>
#include <random>
#include <thread>
#include <utility>
#include <vector>

#include "libfly/env/catalogue.hpp"
#include "libfly/env/geometry.hpp"
#include "libfly/env/mechanisms.hpp"
#include "libfly/io/gsd.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle/dimer.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/VoS.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/hessian.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/lattice.hpp"
#include "libfly/utility/random.hpp"

/**
 * \file find.hpp
 *
 * @brief Class to coordinate a group of saddle-point searches.
 */

namespace fly::saddle {

  /////////// static functions ///////////

  // Reconstruct sp/minima according to reference geometry
  static system::SoA<Position> reconstruct(env::Geometry<Index> const& ref,
                                           system::VoS<Delta> const& m,
                                           system::SoA<Position const&> in) {
    system::SoA<Position> sp(in);

    for (int j = 0; j < m.size(); ++j) {
      sp(r_, ref[j][i_]) += m[j][del_];
    }

    return sp;
  }

  /**
   * @brief Reconstruct and relax a mechanism's saddle point and minima.
   */
  Master::Recon Master::recon_relax(env::Geometry<Index> const& geo,
                                    env::Mechanism const& m,
                                    system::SoA<Position const&, TypeID const&, Frozen const&> in) {
    //
    struct Recon res {
      reconstruct(geo, m.delta_fwd, in), reconstruct(geo, m.delta_sp, in), {}, {},
    };

    auto& thr = thread();

    {  // Minima
      system::SoA<Position, PotentialGradient, TypeID const&, Frozen const&> minima(in.size());
      minima.rebind(id_, in);
      minima.rebind(fzn_, in);
      minima[r_] = res.min[r_];

      if (!thr.min.minimise(minima, minima, thr.pot, 1)) {
        res.rel_min = std::move(minima);
      }
    }

    {  // SP

      system::SoA<Position, Axis, TypeID const&, Frozen const&> dim(in.size());
      dim.rebind(id_, in);
      dim.rebind(fzn_, in);

      // Ensures Axis transforms with the symmetry
      dim[ax_] = res.sp[r_] - in[r_];
      // ... and normalise.
      dim[ax_] /= gnorm(dim[ax_]);
      // Align res.sp post-use for constructing axis!
      if (m_count_frozen == 0) {
        centroid_align(res.sp, in);
      }

      dim[r_] = res.sp[r_];

      if (!thr.dimer.find_sp(dim, dim, in, thr.pot, {}, m_count_frozen, 1)) {
        res.rel_sp = std::move(dim);
      }
    }

    return res;
  }

  // Compute angle between VoSs (in degrees)
  static double theta_mech(system::VoS<Delta> const& a, system::VoS<Delta> const& b) {
    //
    verify(a.size() == b.size(), "");

    auto n = a.size();

    double ct = 0;
    double sa = 0;
    double sb = 0;

    for (int i = 0; i < n; i++) {
      ct += gdot(a[i][del_], b[i][del_]);
      sa += gnorm_sq(a[i][del_]);
      sb += gnorm_sq(b[i][del_]);
    }

    return std::acos(std::clamp(ct / std::sqrt(sa * sb), -1., 1.)) / 2 / M_PI * 360.;
  }

  // Find the indices of minimally and maximally separated atoms.
  static std::pair<int, int> min_max(system::SoA<Position const&> a, system::SoA<Position const&> b) {
    //
    ASSERT(a.size() > 0, "Input has only {} atoms?", a.size());
    ASSERT(a.size() == b.size(), "min_max()'s inputs are different lengths {} != {}", a.size(), b.size());

    int min = 0;
    int max = 0;

    double r_min = std::numeric_limits<double>::max();
    double r_max = 0;

    for (int i = 0; i < a.size(); i++) {
      //
      double r = gnorm_sq(a(r_, i) - b(r_, i));

      if (r < r_min) {
        min = i;
        r_min = r;
      }
      if (r > r_max) {
        max = i;
        r_max = r;
      }
    }

    return {min, max};
  }

  ////////////////////////////////////////

  Master::Master(Master::Options const& opt,
                 system::Box const& box,
                 potential::Generic const& pot,
                 minimise::LBFGS const& min,
                 saddle::Dimer const& dimer)
      : m_opt{opt}, m_box{box} {
    //
    std::random_device rd;

    Xoshiro prng(rd);

    for (int i = 0; i < opt.num_threads; i++) {
      m_data.push_back({pot, min, dimer, prng, {box, pot.r_cut()}, {}});
      prng.long_jump();
    }
  }

  std::vector<Master::LocalisedGeo> Master::package(std::vector<int> const& ix,
                                                    env::Catalogue const& cat,
                                                    int num_threads) {
    //
    std::vector<LocalisedGeo> out_data(ix.size());

#pragma omp parallel for num_threads(num_threads) schedule(static)
    for (std::size_t i = 0; i < ix.size(); i++) {
      //
      out_data[i].geo = cat.get_geo(ix[i]);
      out_data[i].centre = out_data[i].geo[0][i_];
      out_data[i].syms = cat.calc_self_syms(ix[i]);
    }

    return out_data;
  }

  std::vector<Master::Found> Master::find_mechs(std::vector<LocalisedGeo> const& geo_data,
                                                SoA in,
                                                std::optional<Hint> const& hint) {
    // Early exit!
    if (geo_data.empty()) {
      return {};
    }

    neigh::List nl_pert{m_box, m_opt.r_pert};

    nl_pert.rebuild(in, m_opt.num_threads);

    calc_minima_hess(in);

    {
      m_sep_list.resize(safe_cast<std::size_t>(in.size()));

      auto mi = m_box.slow_min_image_computer();

#pragma omp parallel for num_threads(m_opt.num_threads)
      for (Eigen::Index i = 0; i < in.size(); i++) {
        //

        Eigen::Index k = 0;
        double max = 0;

        for (Eigen::Index j = 0; j < in.size(); j++) {
          if (double dr = mi(in(r_, i), in(r_, j)); dr > max) {
            k = j;
            max = dr;
          }
        }

        m_sep_list[safe_cast<std::size_t>(i)] = k;
      }
    }

    std::vector<Found> out(geo_data.size());

//
#pragma omp parallel
#pragma omp single nowait
    {
      for (std::size_t j = 0; j < geo_data.size(); ++j) {
#pragma omp task untied default(none) firstprivate(j) shared(out, in, nl_pert, geo_data, hint)
        {
          find_n(out[j], geo_data[j], in, nl_pert, hint);

          if (out[j]) {
            std::vector<double> barriers;

            for (auto const& m : out[j].mechs()) {
              barriers.push_back(m.barrier);
            }

            std::sort(barriers.begin(), barriers.end());

            auto last = std::unique(barriers.begin(), barriers.end());

            barriers.erase(last, barriers.end());

            fmt::print("FINDER: @{:<4} found {:>3} mech(s) {::.3f}\n",
                       geo_data[j].centre,
                       out[j].mechs().size(),
                       barriers);

          } else {
            fmt::print("FINDER: @{:<4} symmetry-break!\n", geo_data[j].centre);
          }
        }
      }
    }

    // Compute pre factors
    for (auto& f : out) {
      if (f) {
        for (auto& m : f.m_mechs) {
          m.kinetic_pre
              = std::sqrt(std::exp(m_log_prod_eigen - m.kinetic_pre) / (2 * M_PI * 1.6605390666050e-27));
        }
      }
    }

    if (m_opt.fout || m_opt.debug) {
      auto c = m_opt.fout ? m_opt.fout->n_frames() : 0;

      for (std::size_t i = 0; i < out.size(); ++i) {
        if (out[i]) {
          for (auto& mech : out[i].m_mechs) {
            fmt::print(stderr,
                       "FINDER: @{} frame={}, Delta={}eV, A={:e}Hz\n",
                       geo_data[i].centre,
                       c,
                       mech.barrier,
                       mech.kinetic_pre);

            if (m_opt.fout) {
              m_opt.fout->commit(
                  [&] { m_opt.fout->write(r_, reconstruct(geo_data[i].geo, mech.delta_sp, in)); });
              m_opt.fout->commit(
                  [&] { m_opt.fout->write(r_, reconstruct(geo_data[i].geo, mech.delta_fwd, in)); });
            }

            c += m_opt.fout ? 2 : 1;
          }
        }
      }
    }

    return out;
  }

  void Master::process_hint(Found& out,
                            LocalisedGeo const& geo_data,
                            SoA in,
                            Hint const& hint,
                            std::vector<system::SoA<Position>>& cache) {
    // This is copypasta from recon_relax as we need the axis.

    system::SoA<Position> re_sp = reconstruct(hint.geo, hint.delta_sp, hint.prev_state);

    system::SoA<Position, Axis, TypeID const&, Frozen const&> dimer(in.size());

    dimer.rebind(id_, in);
    dimer.rebind(fzn_, in);

    // Ensures Axis transforms with the symmetry
    dimer[ax_] = re_sp[r_] - hint.prev_state[r_];
    // ... and normalise.
    dimer[ax_] /= gnorm(dimer[ax_]);
    // Align post-use for constructing axis!
    if (m_count_frozen == 0) {
      centroid_align(re_sp, hint.prev_state);
    }
    dimer[r_] = re_sp[r_];

    if (thread().dimer.find_sp(dimer, dimer, hint.prev_state, thread().pot, {}, m_count_frozen, 1)) {
      fmt::print("FINDER: Hint's find_sp failed!\n");
      return;
    }

    if (std::optional mech = saddle_2_mech(in, dimer, geo_data.geo)) {
      double b = mech->barrier;

      verify(add_mech(out, in, std::move(*mech), dimer, cache, geo_data), "Hinted mech was not added");

      if (out.m_fail) {
        fmt::print("FINDER: Hint reveals asym!\n");
        return;
      }

      fmt::print("FINDER: @{:<4} hints {:>3} mech(s) [{:.3f}]\n", geo_data.centre, out.mechs().size(), b);
      return;
    }

    fmt::print("FINDER: Hint was useless!\n");

    // Debugging output

    static int count = 0;

    fly::io::BinaryFile file(fmt::format("build/gsd/hint.{}.gsd", count++), fly::io::create);

    file.commit([&] {
      file.write("particles/N", fly::safe_cast<std::uint32_t>(in.size()));
      file.write(m_box);
      file.write(id_, in);
      file.write(r_, hint.prev_state);
    });

    file.commit([&] { file.write(r_, re_sp); });

    file.commit([&] { file.write(r_, dimer); });

    file.commit([&] { file.write(r_, in); });
  }

  ///////////////////////////////////////////////////////

  void Master::find_n(Found& out,
                      LocalisedGeo const& geo_data,
                      SoA in,
                      neigh::List const& nl_pert,
                      std::optional<Hint> const& hint) {
    // Cache saddle-points.
    std::vector<system::SoA<Position>> cache;

    if (hint && hint->centre == geo_data.centre) {
      process_hint(out, geo_data, in, *hint, cache);

      if (out.m_fail) {
        return;
      }
    }

    std::vector<Batch> batch(safe_cast<std::size_t>(m_opt.batch_size), Batch(in.size()));

    int tot = 0;
    int fail = 0;

    int mod = in(id_, geo_data.centre) == 0 ? 1 : 4;

    while (tot < m_opt.max_searches * mod && fail < m_opt.max_failed_searches * mod) {
      //
      if (find_batch(tot, out, batch, geo_data, in, nl_pert, cache)) {
        fail = 0;
      } else {
        fail += m_opt.batch_size;
      }

      if (out.m_fail) {
        return;
      }

      tot += m_opt.batch_size;

      if (m_opt.debug) {
        //
        double t_min = 360;

        for (auto const& a : out.m_mechs) {
          for (auto const& b : out.m_mechs) {
            if (&a != &b) {
              t_min = std::min(t_min, theta_mech(a.delta_sp, b.delta_sp));
            }
          }
        }

        fmt::print("FINDER: {} mech(s) @{}: fail={}, tot={}, t_min={}deg, n-cache={}\n",
                   out.m_mechs.size(),
                   geo_data.centre,
                   fail,
                   tot,
                   t_min,
                   cache.size());
      }
    }
  }

  void Master::dump_recon(system::SoA<Position const&, TypeID const&> in,
                          Index::scalar_t centre,
                          Recon const& recon,
                          system::SoA<Position const&> dimer) const {
    //
    constexpr char const* fname = "crash.check.gsd";

    fly::io::BinaryFile file(fname, fly::io::create);

    fmt::print(stderr,
               "ERROR: First symmetry (identity) should be guaranteed to reconstruct at atom {}. Failure was "
               "written to \"{}\" in the current working directory whose frames are: initial configuration, "
               "dimer final configuration, reconstructed saddle-point, reconstructed final-minima.\n",
               centre,
               fname);

    file.commit([&] {
      file.write("particles/N", fly::safe_cast<std::uint32_t>(in.size()));
      file.write(m_box);
      file.write(id_, in);
      file.write(r_, in);
    });

    file.commit([&] { file.write(r_, dimer); });

    file.commit([&] { file.write(r_, recon.sp); });

    file.commit([&] { file.write(r_, recon.min); });

    if (recon.rel_sp) {
      fmt::print(stderr, "ERROR: and relaxed saddle, \n");
      file.commit([&] { file.write(r_, *recon.rel_sp); });
    }

    if (recon.rel_min) {
      fmt::print(stderr, "ERROR: and relaxed min, \n");
      file.commit([&] { file.write(r_, *recon.rel_min); });
    }
  }

  std::optional<system::SoA<Position>> Master::check_mech(Found& out,
                                                          system::SoA<Position const&> dimer,
                                                          env::Mechanism const& mech,
                                                          std::size_t sym_indx,
                                                          LocalisedGeo const& geo_data,
                                                          SoA in) {
    bool fail;

#pragma omp atomic read
    fail = out.m_fail;

    if (fail) {
      return {};
    }

    //
    auto recon = recon_relax(geo_data.geo, mech, in);

    auto set_fail_flag = [&, sym_indx] {
      if (sym_indx == 0) {
#pragma omp critical
        dump_recon(in, geo_data.centre, recon, dimer);

        throw error("First symmetry (identity) should be guaranteed to reconstruct");
      }

#pragma omp atomic write
      out.m_fail = true;
    };

    // Common

    ThreadData& thr = thread();

    thr.pot_nl.rebuild(in);
    double E0 = thr.pot.energy(in, thr.pot_nl, 1);

    // Check minima reconstructs.

    if (!recon.rel_min) {
      dprint(m_opt.debug, "FINDER: @ symmetry #{} recon_relax-min failied\n", sym_indx);
      set_fail_flag();
      return {};
    }

    centroid_align(*recon.rel_min, recon.min);

    thr.pot_nl.rebuild(*recon.rel_min);
    double Ef = thr.pot.energy(in, thr.pot_nl, 1);

    double re_delta = Ef - E0;

    double d_delta = std::abs(re_delta - mech.delta);

    double frac_delta = d_delta / std::abs(mech.delta);

    double err_fwd = gnorm((*recon.rel_min)[r_] - recon.min[r_]);

    double d_fwd = std::abs(mech.err_fwd - err_fwd);

    double frac_fwd = d_fwd / mech.err_fwd;

    dprint(m_opt.debug,
           "FINDER: Recon @{:>4} sym={:<2}  Δ(ΔE)={:.3f} [{:.4f}] ΔR_min={:.3f} [{:.4f}]\n",
           geo_data.centre,
           sym_indx,
           d_delta,
           frac_delta,
           d_fwd,
           frac_fwd);

    if ((d_delta > m_opt.recon_e_tol_abs && frac_delta > m_opt.recon_e_tol_frac)
        || (frac_fwd > m_opt.recon_norm_frac_tol && d_fwd > m_opt.recon_norm_abs_tol)) {
      set_fail_flag();
      return {};
    }

    // Check SP reconstructs if expecting it to.

    if (mech.poison_sp) {
      dprint(m_opt.debug, "FINDER: @ symmetry #{} exiting from poison-sp\n", sym_indx);
      return {};
    }

    if (!recon.rel_sp) {
      dprint(m_opt.debug, "FINDER: @ symmetry #{} recon_relax-sp\n", sym_indx);
      set_fail_flag();
      return {};
    }

    centroid_align(*recon.rel_sp, recon.sp);

    thr.pot_nl.rebuild(*recon.rel_sp);
    double Esp = thr.pot.energy(in, thr.pot_nl, 1);

    double re_barrier = Esp - E0;

    double d_barrier = std::abs(re_barrier - mech.barrier);

    double frac_barrier = d_barrier / std::abs(mech.barrier);

    double err_sp = gnorm((*recon.rel_sp)[r_] - recon.sp[r_]);

    double d_sp = std::abs(mech.err_sp - err_sp);

    double frac_sp = d_sp / mech.err_sp + err_sp;

    dprint(m_opt.debug,
           "FINDER: Recon @{:>4} sym={:<2} Δ(ΔE#)={:.3f} [{:.4f}]  ΔR_sp={:.3f} [{:.4f}]\n",
           geo_data.centre,
           sym_indx,
           d_barrier,
           frac_barrier,
           d_sp,
           frac_sp);

    if ((d_barrier > m_opt.recon_e_tol_abs && frac_barrier > m_opt.recon_e_tol_frac)
        || (frac_sp > m_opt.recon_norm_frac_tol && d_sp > m_opt.recon_norm_abs_tol)) {
      set_fail_flag();
      return {};
    } else {
      if (m_count_frozen == 0) {
        centroid_align(*recon.rel_sp, in);
      }
      return std::move(recon.rel_sp);
    }
  }

  // Add mech should: add one or more mechs and return true or return false if no mech was added due to a mech
  // collision
  bool Master::add_mech(Found& out,
                        SoA in,
                        env::Mechanism&& mech,
                        system::viewSoA<Position, Axis> dimer,
                        std::vector<system::SoA<Position>>& cache,
                        LocalisedGeo const& geo_data) {
    // Found a mechanism.

    if (!is_new_mech(mech, out.m_mechs)) {
      dprint(m_opt.debug, "FINDER: Duplicate mech, caching\n");
      cache.emplace_back(dimer);
      return false;
    }

    // Found a new mechanism.

    if (mech.poison_fwd) {
      dprint(m_opt.debug, "FINDER: caching poisoned\n");
      out.m_mechs.push_back(std::move(mech));
      cache.emplace_back(dimer);
      return true;
    }

    // Found a new non-poisoned forward mechanism (sp may be poisoned).

    std::size_t num_new;

    try {
      // Adds mechanism + all symmetries
      num_new = append_syms(geo_data.syms, mech, out.m_mechs);
    } catch (...) {
#pragma omp critical
      {
        fly::io::BinaryFile file("crash.gsd", fly::io::create);

        fmt::print(
            stderr,
            "Append symmetries threw @{}, this implies a symmetrical saddle-point but sym-breaking minima, "
            "perhaps the minimiser reached the wrong minima i.e. failed to follow the minimum mode or the "
            "sp was higher order. Failure was written to \"crash.gsd\" in the current working directory\n",
            geo_data.centre);

        file.commit([&] {
          file.write("particles/N", fly::safe_cast<std::uint32_t>(in.size()));
          file.write(r_, in);
          file.write(m_box);
          file.write(id_, in);
        });
        file.commit([&] { file.write(r_, reconstruct(geo_data.geo, mech.delta_sp, in)); });
        file.commit([&] { file.write(r_, dimer); });
        file.commit([&] { file.write(r_, reconstruct(geo_data.geo, mech.delta_fwd, in)); });
      }
      throw;
    }

    verify(num_new > 0, "Found none @{} but should have had at least 1", geo_data.centre);

    std::size_t n = cache.size();

    if (mech.poison_sp) {
      cache.emplace_back(dimer);
    } else {
      cache.resize(n + num_new);
    }

    for (std::size_t k = 0; k < num_new; k++) {
#pragma omp task untied default(none) firstprivate(n, k, num_new) \
    shared(out, in, cache, geo_data, dimer, mech)
      {
        std::optional recon_sp
            = check_mech(out, dimer, out.m_mechs[out.m_mechs.size() - num_new + k], k, geo_data, in);

        if (mech.poison_sp) {
          // This check may be too strong
          ASSERT(!recon_sp, "Logic error, check should have returned false", 0);

        } else {
          if (recon_sp) {
            cache[n + k] = std::move(*recon_sp);
          } else {
            bool fail;
#pragma omp atomic read
            fail = out.m_fail;
            // double check that check_mech has set fail flag
            ASSERT(fail, "Logic error if no recon_sp then fail should be set", 0);
          }
        }
      }
    }
#pragma omp taskwait

    return true;
  }

  bool Master::find_batch(int tot,
                          Found& out,
                          std::vector<Batch>& batch,
                          LocalisedGeo const& geo_data,
                          SoA in,
                          neigh::List const& nl_pert,
                          std::vector<system::SoA<Position>>& cache) {
    // Abort SPS if the cosine of the angle between the dimer and a known SP is greater than this.
    double theta_tol = ((30 - 2.5) * std::exp(-0.02 * tot) + 2.5) / 360. * 2. * M_PI;

    dprint(m_opt.debug, "FINDER: Theta tolerance = {}\n", theta_tol / 2. / M_PI * 360.);

    //  Do batch_size SP searches
    for (Batch& elem : batch) {
#pragma omp task untied default(none) firstprivate(theta_tol) \
    shared(elem, in, nl_pert, batch, cache, geo_data)
      {
        perturb(elem.dimer, in, geo_data.centre, nl_pert);

        if (m_count_frozen == 0) {
          centroid_align(elem.dimer, in);
        }

        system::SoA<Position&, Axis&, Frozen const&, TypeID const&> dimer(in.size());

        dimer.rebind(r_, elem.dimer);
        dimer.rebind(ax_, elem.dimer);
        dimer.rebind(fzn_, in);
        dimer.rebind(id_, in);

        elem.exit = find_sp(dimer, in, cache, theta_tol);

        if (elem.exit == Dimer::success) {
          elem.mech = saddle_2_mech(in, dimer, geo_data.geo);
        } else {
          elem.mech = std::nullopt;
        }
      }
    }
#pragma omp taskwait

    bool found_one_or_more = false;

    // Process batch
    for (auto&& elem : batch) {
      if (!elem.mech) {
        //
        ASSERT(elem.exit != Dimer::Exit::uninit, "Uninitialised exit/return code exit={}!", elem.exit);

        if (elem.exit == Dimer::Exit::success) {
          dprint(m_opt.debug, "FINDER: Caching failure not at this atom\n");
          cache.emplace_back(elem.dimer);
        } else {
          dprint(m_opt.debug, "FINDER: Not caching SPS collision/fail\n");
        }

      } else if (add_mech(out, in, std::move(*elem.mech), elem.dimer, cache, geo_data)) {
        found_one_or_more = true;
      }
    }

    return found_one_or_more;
  }

  std::optional<env::Mechanism> Master::saddle_2_mech(system::viewSoA<Position, Frozen, TypeID> in,
                                                      system::viewSoA<Position, Axis, Frozen, TypeID> dimer,
                                                      env::Geometry<Index> const& geo) {
    // Saddle search

    std::optional path = do_adj_min(dimer, in, geo[0][i_]);

    if (!path) {
      return {};
    }

    auto& [rev, sp, fwd] = *path;

    //  CHECK pathway is min->sp->min

    double i2r = gnorm(rev[r_] - in[r_]);
    double r2f = gnorm(fwd[r_] - rev[r_]);
    double r2w = gnorm(rev[r_] - sp[r_]);
    double w2f = gnorm(fwd[r_] - sp[r_]);

    dprint(m_opt.debug, "FINDER: i2r={:.4f}, r2f={:.4f} r2sp={:.4f} sp2f={:.4f}\n", i2r, r2f, r2w, w2f);

    if (i2r > m_opt.basin_tol) {
      dprint(m_opt.debug, "FINDER: Mech starting {}A from initial is unlinked\n", i2r);
      return {};
    }
    if (r2f < m_opt.basin_tol) {
      dprint(m_opt.debug, "FINDER: Mech total displacement={} => converged back to initial\n", r2f);
      return {};
    }
    if (r2w < m_opt.stationary_tol) {
      dprint(m_opt.debug, "FINDER: Min->Sp displacement={}\n", r2w);
      return {};
    }
    if (w2f < m_opt.stationary_tol) {
      dprint(m_opt.debug, "FINDER: Sp->Min displacement={}\n", w2f);
      return {};
    }

    // We now have a mech via min->sp->min, time to build a mechanism.

    ThreadData& thr = thread();

    thr.pot_nl.rebuild(rev);
    double E0 = thr.pot.energy(in, thr.pot_nl, 1);

    thr.pot_nl.rebuild(sp);
    double Esp = thr.pot.energy(in, thr.pot_nl, 1);

    thr.pot_nl.rebuild(fwd);
    double Ef = thr.pot.energy(in, thr.pot_nl, 1);

    env::Mechanism mech;

    mech.barrier = Esp - E0;
    mech.delta = Ef - E0;

    for (auto const& atom : geo) {
      //
      auto i = atom[i_];

      mech.delta_sp.emplace_back(sp(r_, i) - rev(r_, i));
      mech.delta_fwd.emplace_back(fwd(r_, i) - rev(r_, i));
    }

    //////////////// Partial hessian compute. ////////////////

    thr.pot_nl.rebuild(sp, 1);

    thr.pot.hessian(thr.hess, in, thr.pot_nl, 1);

    thr.pot.mw_hessian(thr.hess, in, 1);

    system::Hessian::Vector const& freq = thr.hess.eigenvalues();

    // Must have at least one neg
    verify(freq[0] < -m_opt.hessian_eigen_zero_tol, "Saddle-point with minimum mode={}", freq[0]);

    int count_zeros = 0;
    mech.kinetic_pre = 0;

    for (int i = 1; i < freq.size(); i++) {
      //
      if (freq[i] < -m_opt.hessian_eigen_zero_tol) {
        dprint(m_opt.debug, "FINDER: Second order SP or higher, modes {} = {}\n", i, freq.head(i));
        return {};
      } else if (freq[i] > m_opt.hessian_eigen_zero_tol) {
        mech.kinetic_pre += std::log(freq[i]);
      } else {
        ++count_zeros;
      }
    }

    if (count_zeros != m_num_zero_modes) {
      for (int j = 0; j < std::max(m_num_zero_modes, count_zeros) + 3; j++) {
        fmt::print(stderr, "FINDER: [ERR] min mode {} = {}\n", j, freq[j]);
      }
      throw error("Expecting {} zero modes but hessian has {}", m_num_zero_modes, count_zeros);
    }

    // Test reconstruction of the mechanism, mark poisoning and calc recon errors

    system::SoA<Position const&, TypeID const&, Frozen const&> ref(rev.size());
    ref.rebind(r_, rev);
    ref.rebind(fzn_, in);
    ref.rebind(id_, in);

    auto [re_min, re_sp, rel_min, rel_sp] = recon_relax(geo, mech, ref);

    double const r_tol = m_opt.capture_r_tol;
    double const e_tol = m_opt.capture_E_tol;

    if (!rel_min) {
      dprint(m_opt.debug, "FINDER: Mech min-poisoned, could not relax!\n");
      mech.poison_fwd = true;
    } else {
      mech.err_fwd = gnorm((*rel_min)[r_] - re_min[r_]);

      dprint(m_opt.debug, "FINDER: In-place reconstruct, unaligned: err_fwd={}\n", mech.err_fwd);

      centroid_align(*rel_min, fwd);

      thr.pot_nl.rebuild(*rel_min);
      double re_Ef = thr.pot.energy(in, thr.pot_nl, 1);

      double err_fwd = gnorm((*rel_min)[r_] - fwd[r_]);

      double del_Efwd = std::abs(re_Ef - Ef);

      dprint(m_opt.debug, "FINDER: In-place reconstruct, aligned: err_fwd={} dEfwd={}\n", err_fwd, del_Efwd);

      if (err_fwd > r_tol || del_Efwd > e_tol) {
        dprint(m_opt.debug, "FINDER: Mech min-poisoned!\n");
        mech.poison_fwd = true;
      } else {
        mech.poison_fwd = false;
      }
    }

    if (!rel_sp) {
      dprint(m_opt.debug, "FINDER: Mech sp-poisoned, could not relax!\n");
      mech.poison_sp = true;
    } else {
      mech.err_sp = gnorm((*rel_sp)[r_] - re_sp[r_]);

      dprint(m_opt.debug, "FINDER: In-place reconstruct, unaligned: err_sp={}\n", mech.err_sp);

      centroid_align(*rel_sp, sp);

      thr.pot_nl.rebuild(*rel_sp);
      double re_Esp = thr.pot.energy(in, thr.pot_nl, 1);

      double err_sp = gnorm((*rel_sp)[r_] - sp[r_]);

      double del_Esp = std::abs(re_Esp - Esp);

      dprint(m_opt.debug, "FINDER: In-place reconstruct, aligned: err_sp={} dEsp={}\n", err_sp, del_Esp);

      if (err_sp > r_tol || del_Esp > e_tol) {
        dprint(m_opt.debug, "FINDER: Mech sp-poisoned!\n");
        mech.poison_sp = true;
      } else {
        mech.poison_sp = false;
      }
    }

    // Debugging extreme poisoning

    if (mech.poison_fwd && mech.barrier < 1.) {
#pragma omp critical
      {
        fly::io::BinaryFile file("poison.gsd", fly::io::create);

        file.commit([&] {
          file.write("particles/N", fly::safe_cast<std::uint32_t>(in.size()));
          file.write(m_box);
          file.write(r_, in);
          file.write(id_, in);
        });

        file.commit([&, rev = rev] { file.write(r_, rev); });
        file.commit([&, sp = sp] { file.write(r_, sp); });
        file.commit([&, fwd = fwd] { file.write(r_, fwd); });

        file.commit([&] { file.write(id_, in); });
        file.commit([&, re_sp = re_sp] { file.write(r_, re_sp); });
        file.commit([&, re_min = re_min] { file.write(r_, re_min); });

        file.commit([&] { file.write(id_, in); });

        if (rel_sp) {
          file.commit([&, rel_sp = rel_sp] { file.write(r_, *rel_sp); });
        }

        if (rel_min) {
          file.commit([&, rel_min = rel_min] { file.write(r_, *rel_min); });
        }
      }

      throw error("Mechanism @{} with energy barrier = {}eV is poisoned! err_sp={}, err_min={}",
                  geo[0][i_],
                  mech.barrier,
                  mech.err_sp,
                  mech.err_fwd);
    }

    if (m_opt.debug) {
      fmt::print(
          "FINDER: Mech ΔEsp={:.3e}, ΔEfwd={:.3e}, N_zeros={}\n", mech.barrier, mech.delta, count_zeros);

      for (int i = 0; i < m_num_zero_modes + 3; i++) {
        fmt::print("FINDER: sp mode {} = {}\n", i, freq[i]);
      }
    }

    return mech;
  }  // namespace fly::saddle

  Dimer::Exit Master::find_sp(system::SoA<Position&, Axis&, Frozen const&, TypeID const&> dimer,
                              system::SoA<Position const&> in,
                              std::vector<system::SoA<Position>> const& hist_sp,
                              double theta_tol) {
    //
    ThreadData& thr = thread();

    auto err = thr.dimer.find_sp(dimer, dimer, in, thr.pot, hist_sp, theta_tol, m_count_frozen, 1);

    if (err && m_opt.debug) {
      switch (err) {
        case Dimer::collision:
          fmt::print("FINDER: SPS fail - collision\n");
          break;
        case Dimer::iter_max:
          fmt::print("FINDER: SPS fail - iter_max\n");
          break;
        case Dimer::convex:
          fmt::print("FINDER: SPS fail - convex\n");
          break;
        default:
          ASSERT(false, "Unknown error = {}", err);
      }
    }

    return err;
  }

  /**
   * @brief Given a saddle point produce a min->sp->min pathway
   */
  std::optional<Master::Pathway> Master::do_adj_min(
      system::SoA<Position const&, Axis const&, Frozen const&, TypeID const&> dimer,
      system::SoA<Position const&> in,
      Index::scalar_t centre) {
    //   Minimisations

    double disp = gnorm(dimer[r_] - in[r_]);

    system::SoA<Position, PotentialGradient, Frozen, TypeID const&> relax{dimer.size()};
    relax[r_] = dimer[r_] + dimer[ax_] * disp * m_opt.nudge_frac;
    relax[fzn_] = dimer[fzn_];
    relax.rebind(id_, dimer);

    if (m_count_frozen == 0) {
      // Freeze furthest atom to remove translational degrees of freedom.

      auto min = m_sep_list[safe_cast<size_t>(centre)];

      dprint(m_opt.debug, "FINDER: Freezing atom #{}\n", min);

      relax(fzn_, min) = true;
    }

    ThreadData& thr = thread();

    if (thr.min.minimise(relax, relax, thr.pot, 1)) {
      dprint(m_opt.debug, "FINDER: minimisation failed\n");
      return {};
    }

    system::SoA<Position> fwd{relax};

    relax[r_] = dimer[r_] - dimer[ax_] * disp * m_opt.nudge_frac;

    if (thr.min.minimise(relax, relax, thr.pot, 1)) {
      dprint(m_opt.debug, "FINDER: minimisation failed\n");
      return {};
    }

    system::SoA<Position> rev{relax};

    // Swap if forward went to final.
    if (gnorm(rev[r_] - in[r_]) > gnorm(fwd[r_] - in[r_])) {
      using std::swap;
      swap(fwd, rev);
    }

    // We now have a min->sp->min pathway with no translation (due to freeze). Need to correct for COM drift
    // during SPS.

    Vec in_com = centroid(in);
    Vec rev_com = centroid(rev);

    Vec drift = rev_com - in_com;

    dprint(m_opt.debug, "FINDER: Centre of mass drift={}\n", gnorm(drift));

    system::SoA<Position> sp_{dimer};

    for (int i = 0; i < in.size(); i++) {
      fwd(r_, i) -= drift;
      sp_(r_, i) -= drift;
      rev(r_, i) -= drift;
    }

    ASSERT(gnorm(centroid(rev) - centroid(in)) < 1e-8, "Drift correction failed", 0);

    // Check mech centred on centre.

    auto [_, max] = min_max(rev, fwd);

    if (max != centre) {
      dprint(m_opt.debug, "FINDER: found mech at {} but wanted {}\n", max, centre);
      return {};
    }

    return Pathway{std::move(rev), std::move(sp_), std::move(fwd)};
  }

  ///////////////////////////////////////////////////////////

  void Master::calc_minima_hess(system::SoA<Position const&, Frozen const&, TypeID const&> in) {
    //

    m_count_frozen = 0;

    for (int i = 0; i < in.size(); i++) {
      if (in(fzn_, i)) {
        m_count_frozen++;
      }
    }
    // If all unfrozen then zero modes due to translational degrees of freedom.
    int exp_zero_modes = m_count_frozen > 0 ? spatial_dims * m_count_frozen : spatial_dims;

    m_data[0].pot_nl.rebuild(in);

    m_data[0].pot.hessian(m_data[0].hess, in, m_data[0].pot_nl, m_opt.num_threads);

    m_data[0].pot.mw_hessian(m_data[0].hess, in, m_opt.num_threads);

    system::Hessian::Vector const& freq = m_data[0].hess.eigenvalues();

    m_log_prod_eigen = 0;
    m_num_zero_modes = 0;

    for (int i = 0; i < freq.size(); i++) {
      // Negative modes = not a minima.
      verify(freq[i] > -m_opt.hessian_eigen_zero_tol, "Mode {} has eigenvalue = {}", i, freq[i]);

      // Sum non-zero modes.
      if (freq[i] > m_opt.hessian_eigen_zero_tol) {
        m_log_prod_eigen += std::log(freq[i]);
      } else {
        ++m_num_zero_modes;
      }
    }

    if (m_num_zero_modes != exp_zero_modes) {
      for (int j = 0; j < std::max(m_num_zero_modes, exp_zero_modes) + 5; j++) {
        fmt::print(stderr, "FINDER: [ERR ]min mode {} = {}\n", j, freq[j]);
      }
      throw error("Expecting {} zero modes but fount {}", exp_zero_modes, m_num_zero_modes);
    }

    if (m_opt.debug) {
      fmt::print(
          "FINDER: Null space of Hessian has dimension  = {}, exp = {}\n", m_num_zero_modes, exp_zero_modes);

      for (int i = 0; i < m_num_zero_modes + 5; i++) {
        fmt::print("FINDER: Initial min mode {} = {}\n", i, freq[i]);
      }
    }
  }

  // Perturb in-place positions around centre and write axis,
  auto Master::perturb(system::SoA<Position&, Axis&> out,
                       system::SoA<Position const&, Frozen const&> in,
                       Index::scalar_t centre,
                       neigh::List const& nl) -> double {
    //
    Xoshiro& prng = thread().prng;

    std::normal_distribution normal(0., 1.);

    std::normal_distribution prime(m_opt.stddev, m_opt.stddev / 3.0);

    double stddev = -1;

    while (stddev < 0) {
      stddev = prime(prng);
    }

    dprint(m_opt.debug, "FINDER: standard deviation={}\n", stddev);

    std::normal_distribution<double> gauss(0, stddev);

    out[r_] = in[r_];
    out[ax_] = 0;

    nl.for_neighbours(centre, m_opt.r_pert, [&](auto n, double r, auto const&) {
      if (!in(Frozen{}, n)) {
        out(r_, n) += Vec::NullaryExpr([&] { return gauss(prng); }) * (1. - r / m_opt.r_pert);
        out(ax_, n) += Vec::NullaryExpr([&] { return normal(prng); });
      }
    });

    /* Experimental quasi-random centre

    constexpr double g = 1.32471795724474602596;  // Golden ratio
    constexpr double a1 = 1.0 / g;
    constexpr double a2 = 1.0 / (g * g);

    static int n = 0;

    double u = 0.5 + a1 * n - std::floor(0.5 + a1 * n);
    double v = 0.5 + a2 * n - std::floor(0.5 + a2 * n);

    n++;

    ASSERT(u >= 0 && u < 1, "u={}", u);
    ASSERT(u >= 0 && u < 1, "v={}", v);

    //    Lambert projection

    double phi = 2 * M_PI * v;
    double sin_lam = 2 * u - 1;
    double cos_lam = std::cos(std::asin(sin_lam));

    Vec p = {
        cos_lam * std::cos(phi),
        cos_lam * std::sin(phi),
        sin_lam,
    };

    */

    ASSERT(!in(Frozen{}, centre), "perturbation centred on a frozen atom {}", centre);

    out(r_, centre)
        += Vec::NullaryExpr([&] { return gauss(prng); });  // p * std::abs(gauss(prng)) * std::sqrt(3);
    out(ax_, centre) += Vec::NullaryExpr([&] { return normal(prng); });

    out[ax_] /= gnorm(out[ax_]);  // normalize

    return stddev;
  }

  bool Master::is_new_mech(env::Mechanism const& maybe, std::vector<env::Mechanism> const& hist) const {
    // All searches have been around same atom hence orientation is fixed.

    if (hist.empty()) {
      dprint(m_opt.debug, "FINDER: First mech must be new\n");
      return true;
    }

    double min_d_sp = std::numeric_limits<double>::max();
    double min_d_fwd = std::numeric_limits<double>::max();
    double min_d_theta = std::numeric_limits<double>::max();

    for (env::Mechanism const& m : hist) {
      //
      double d_sp = env::rmsd<Delta>(maybe.delta_sp, m.delta_sp);
      double d_fwd = env::rmsd<Delta>(maybe.delta_fwd, m.delta_fwd);

      double dot = 0;
      double n1 = 0;
      double n2 = 0;

      for (Eigen::Index i = 0; i < m.delta_sp.size(); ++i) {
        n1 += gnorm_sq(maybe.delta_sp[i][del_]);
        n2 += gnorm_sq(m.delta_sp[i][del_]);
        dot += gdot(maybe.delta_sp[i][del_], m.delta_sp[i][del_]);
      }

      n1 = std::sqrt(n1);
      n2 = std::sqrt(n2);

      double theta = std::acos(std::clamp(dot / (n1 * n2), -1., 1.)) / (2 * M_PI) * 360;

      min_d_sp = std::min(min_d_sp, d_sp);
      min_d_fwd = std::min(min_d_fwd, d_fwd);
      min_d_theta = std::min(min_d_theta, theta);

      if (d_sp < m_opt.mech_tol && d_fwd < m_opt.mech_tol && theta < 10) {
        dprint(m_opt.debug,
               "FINDER: Mech is NOT new, distance: old-sp={:.3f} old-min={:.3f}, theta={:.1f}, n1={:.3f}, "
               "n2={:.3f}\n",
               d_sp,
               d_fwd,
               theta,
               n1,
               n2);
        return false;
      }
    }

    dprint(m_opt.debug,
           "FINDER: Mech is new, min distance: old-sp={:.3e} old-min={:.3e}, theta-min={:.1f}\n",
           min_d_sp,
           min_d_fwd,
           min_d_theta);

    return true;
  }

  std::size_t Master::append_syms(std::vector<env::Catalogue::SelfSymetry> const& syms,
                                  env::Mechanism const& new_mech,
                                  std::vector<env::Mechanism>& mechs) const {
    //
    env::Mechanism mech = new_mech;

    std::size_t count = 0;
    //
    for (auto const& [O, perm] : syms) {
      for (std::size_t i = 0; i < perm.size(); ++i) {
        mech.delta_sp[Eigen::Index(i)][del_].noalias() = O * new_mech.delta_sp[perm[i]][del_];
        mech.delta_fwd[Eigen::Index(i)][del_].noalias() = O * new_mech.delta_fwd[perm[i]][del_];
      }

      if (is_new_mech(mech, mechs)) {
        mechs.push_back(mech);
        ++count;
      }
    }

    dprint(m_opt.debug, "FINDER: Added {} mechs related by symmetry\n", count);

    return count;
  }

}  // namespace fly::saddle