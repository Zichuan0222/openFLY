

#include <fmt/core.h>
#include <omp.h>

//

#include <algorithm>
#include <array>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "libfly/env/catalogue.hpp"
#include "libfly/io/gsd.hpp"
#include "libfly/kinetic/basin.hpp"
#include "libfly/kinetic/skmc.hpp"
#include "libfly/kinetic/superbasin.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/neigh/sort.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle//find.hpp"
#include "libfly/saddle/dimer.hpp"
#include "libfly/saddle/perturb.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/lattice.hpp"
#include "libfly/utility/random.hpp"

using namespace fly;

template <typename... T>
system::Supercell<system::TypeMap<>, Position, Frozen, T...> bcc_iron_motif(double a = 2.855300) {
  //
  system::TypeMap<> FeH(3);

  FeH.set(0, tp_, "Fe");
  FeH.set(1, tp_, "H");
  FeH.set(2, tp_, "V");

  Mat basis{
      {a, 0, 0},
      {0, a, 0},
      {0, 0, a},
  };

  system::Supercell motif
      = system::make_supercell<Position, Frozen, T...>({basis, Arr<bool>::Constant(true)}, FeH, 2);

  motif[fzn_] = false;
  motif[id_] = 0;

  motif(r_, 0) = Vec::Zero();
  motif(r_, 1) = Vec::Constant(0.5);

  return motif;
}

fly::system::SoA<TypeID, Position> explicit_V(std::vector<Vec> const &vac,
                                              system::viewSoA<TypeID, Position> cell) {
  //
  fly::system::SoA<TypeID, Position> special(cell.size() + fly::ssize(vac));

  special[id_].head(cell.size()) = cell[id_];  // Copy cells types
  special[id_].tail(fly::ssize(vac)) = 2;      // TypeID for vacancies == 2

  special[r_].head(cell[r_].size()) = cell[r_];

  Eigen::Index x = cell.size();

  for (auto const &v : vac) {
    special(r_, x++) = v;
  }

  return special;
}

int main() {
  //

  system::Supercell perfect = motif_to_lattice(bcc_iron_motif(), {6, 6, 6});

  DetectVacancies detect(0.75, perfect.box(), perfect);

  system::Supercell cell = remove_atoms(perfect, {1, 3});

  Vec r_H = {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857 / 4 + 3.14};

  cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen>(1, r_H, false)});

  //   cell = add_atoms(
  //       cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857
  //       / 4 * 2 + 3.14}, false, 0)});

  //   cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {0.992116, 6.01736, 4.56979},
  //   false, 0)});

  //   cell(r_, 431) = Vec{4.5896, 4.5892, 5.91909};
  //   cell(r_, 0) = Vec{4.45211, 4.45172, 4.2526};

  fly::io::BinaryFile file("build/gsd/sim.gsd", fly::io::create);

  file.commit([&] {
    file.write(cell.box());
    file.write(cell.map());

    file.write("particles/N", fly::safe_cast<std::uint32_t>(cell.size()));

    file.write(id_, cell);
    file.write(r_, cell);

    file.write("log/time", -1.);
    file.write("log/energy", -1.);
    file.write("log/barrier", -1.);
    file.write("log/kinetic", -1.);
  });

  kinetic::SKMC runner = {
      {
          .debug = true,
          .fread = "build/gsd/cat.v2.2.bin",
          .opt_cache = {
              .barrier_tol = 0.45,
              .debug = true,
              .opt_basin = {
                  .debug = true,
                  .temp = 500,
              },
              .opt_sb = {
                  .debug = true,
              },
          },
          .opt_master = {
              .num_threads = omp_get_max_threads(),
              .max_searches = 200,
              .max_failed_searches = 75,
              .debug = false,
          }
      },
      cell.box(),
      {
          {},
          cell.box(),
      },
      potential::Generic{
          potential::EAM{
              cell.map(),
              std::make_shared<potential::DataEAM>(
                potential::DataEAM{
                  {
                    .debug = false, 
                    .symmetric = false,
                  }, 
                  std::ifstream{"data/wen.eam.fs"}
                }
              ),
          },
      },
      {
          {},
          {},
          cell.box(),
      },
  };

  auto const min_image = cell.box().slow_min_image_computer();

  double run_time = 0;
  bool v_diss = false;
  bool h_escaped = false;

  runner.skmc(cell,
              omp_get_max_threads(),
              [&](double time,                         ///< Total time just after system at post
                  system::SoA<Position const &> pre,   ///< State just before mech applied
                  double E0,                           ///< Energy of the system in state pre.
                  int,                                 ///< Index of central atom of mechanism
                  env::Mechanism const &mech,          ///< Chosen mechanism
                  system::SoA<Position const &> post,  ///< Final state of system after this iteration / mech
                  double Ef                            ///< Energy of system in state post.
              ) {
                //
                run_time = time;

                fly::system::SoA<TypeID const &, Position const &> tmp(cell.size());

                tmp.rebind(r_, post);
                tmp.rebind(id_, cell);

                std::vector vac = detect.detect_vacancies(tmp);

                fmt::print("Found {} vacancies @{:::.2f}\n", vac.size(), vac);

                // V-dis testing

                double const VV = kruskal_max(vac, min_image);

                fmt::print("MST max V-V = {:.3e}\n", VV);

                v_diss = VV > 5;  // Vacancy-cluster dissociation criterion

                // H-escape testing
                if (auto last = post.size() - 1; cell(id_, last) == 1) {
                  double vh = std::numeric_limits<double>::max();

                  for (auto const &v : vac) {
                    vh = std::min(vh, min_image(post(r_, last), v));
                  }

                  fmt::print("Min V-H = {:.3e}\n", vh);

                  h_escaped = vh > 6;  // H-escape criterion set here.
                }

                // Write to GSD

                file.commit([&] {
                  file.write("particles/N", fly::safe_cast<std::uint32_t>(cell.size()));

                  file.write(r_, pre);

                  file.write("log/energy", E0);
                });

                file.commit([&] {
                  //
                  auto vpost = explicit_V(vac, tmp);

                  file.write("particles/N", fly::safe_cast<std::uint32_t>(vpost.size()));

                  file.write(id_, vpost);
                  file.write(r_, vpost);

                  file.write("log/time", time);
                  file.write("log/energy", Ef);
                  file.write("log/barrier", mech.barrier);
                  file.write("log/kinetic", mech.kinetic_pre);
                });

                fmt::print("Just wrote frame index No. {}\n", file.n_frames() - 1);

                return v_diss;
              });

  fmt::print("It took {:.3e}s to terminate\n", run_time);
  fmt::print("H escaped = {}, cluster dissociated = {}\n", h_escaped, v_diss);

  return 0;
}
