#include <fmt/core.h>
#include <omp.h>

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
#include <filesystem>
#include <chrono>

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
  special[id_].tail(fly::ssize(vac)) = 2;      // TypeID for vacacies == 2

  special[r_].head(cell[r_].size()) = cell[r_];

  Eigen::Index x = cell.size();

  for (auto const &v : vac) {
    special(r_, x++) = v;
  }

  return special;
}


struct Result {
  double v_v;  ///< The maximum of the V-V neighrest-neighbour distances. E.G for each vacancy compute the
               ///< distance to its closest neighbour then v_v is the maximum of these.
  double v_h;  ///< The minimum V-H distance.
};

/**
 * @brief
 */
Result distances(system::Box const &box, std::vector<Vec> const &vac, Vec hy) {
  //
  auto mi = box.slow_min_image_computer();

  Result r{
      0,
      std::numeric_limits<double>::max(),
  };

  for (auto const &v : vac) {
    r.v_h = std::min(r.v_h, mi(hy, v));

    for (auto const &n : vac) {
      r.v_v = std::max(r.v_v, mi(n, v));
    }
  }

  return r;
}


// struct 

struct cluster {
  std::string name;
  std::string binname;
  std::vector<long int> atom_index;
  bool adding_H {};
  bool forVdiss {};
  double diss_criterion {};
};

//// ---------------------------------------------- Clusters ------------------------------------------------////

const cluster V1H_Hdis = {"V1H_Hdis", "v1h", {1}, true, false, 6.0};

const cluster V2_Vdis = {"V2_Vdis", "v2", {1, 3}, false, true, 4.85};
const cluster V2H_Hdis = {"V2H_Hdis", "v2h_v3h", {1, 3}, true, false, 6.0};
const cluster V2H_Vdis = {"V2H_Vdis", "v2h_v3h", {1, 3}, true, true, 4.85};

const cluster V3_Vdis = {"V3_Vdis", "v3", {1, 3, 86}, false, true, 4.5};
const cluster V3H_Hdis = {"V3_Vdis", "v2h_v3h", {1, 3, 86}, false, true, 6.0};
const cluster V3H_Vdis = {"V3H_Vdis", "v2h_v3h", {1, 3, 86}, true, true, 4.5};

//// ---------------------------------------------- -------- ------------------------------------------------////

struct TempVdis {  
  int temp {};        //< Temperature of simulation
  int totalcount {};  //< Total no. of sim at this temp
  int vcount {};      //< The no. of times V-dis before H-dis
};




// This function output current system time
std::string findsystemtime(){

    auto systemtime =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  
    struct tm * now = std::localtime(&systemtime);
    char buffer[80];
    strftime(buffer, 80, "%H-%M-%S", now);

    return buffer;
}

// This function output current system date
std::string findsystemdate(){

    auto systemtime =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  
    struct tm * now = std::localtime(&systemtime);
    char buffer[80];
    strftime(buffer, 80, "%Y-%m-%d", now);

    return buffer;
}



// This function output 3 vectors to a txt file under the name of filename.txt
int exporttotxt(std::string filename,
                std::vector<double> exportvector1,
                std::vector<double> exportvector2,
                std::vector<double> exportvector3)
{
  std::ofstream myfile(filename);
    myfile.precision(20);

  for(int i=0; i < exportvector1.size(); ++i){
    myfile << exportvector1[i] << " " << exportvector2[i] << " " << exportvector3[i] << "\n";
  }

  return 0;
}





// Input: temperature, indices of vacancies 
// Output: lifetime of V-H complex or V cluster

double complex_lifetime(double const& temp, 
                        cluster cluster,
                        std::string gsdfiledirectory,
                        int loopcounter){

  system::Supercell perfect = motif_to_lattice(bcc_iron_motif(), {6, 6, 6});

   DetectVacancies detect(0.75, perfect.box(), perfect);


  // ------------------- delete atoms ------------------- //
  system::Supercell cell = remove_atoms(perfect, cluster.atom_index);               // remove atoms by indices
  
  // Indices of atoms to be removed:
  // 1V - 1
  // 2V - 1, 3
  // 3V - 1, 3, 86

  Vec r_H = {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857 / 4 + 3.14};

  // ------------------- adding H ------------------- //
  if (cluster.adding_H){
    cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen>(1, r_H, false)});       //  add H into lattice
  } 



  //   cell = add_atoms(
  //       cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857
  //       / 4 * 2 + 3.14}, false, 0)});

  //   cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {0.992116, 6.01736, 4.56979},
  //   false, 0)});

  //   cell(r_, 431) = Vec{4.5896, 4.5892, 5.91909};
  //   cell(r_, 0) = Vec{4.45211, 4.45172, 4.2526};


    std::string gsdfilename                                 //<<<<<<<<<<<<<<<<< Use this line to save all .gsd file
        = gsdfiledirectory + "/" 
          + std::to_string(int(temp)) + "K_"            
          + findsystemtime() + ".gsd";                   

  // std::string gsdfilename = gsdfiledirectory + "sim.gsd";    //<<<<<<<<<<<<<< Use this line to save only current .gsd file

  fly::io::BinaryFile file(gsdfilename, fly::io::create);
  
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

          //---------------------------- bin file directory --------------------------//
          .fread = "build/gsd/cat." + cluster.binname + ".bin",      
            // Use a different bin file for with/without H
            // Names of bin files: v1h, v2, v2h_v3h, v3
          .opt_cache = {                            
              .barrier_tol = 0.45,
              .debug = true,
              .opt_basin = {
                  .debug = true,
                  .temp = temp,
              },
              .opt_sb = {
                  .debug = true,
              },
          },
          .opt_master = {
              .hessian_eigen_zero_tol = 1e-3,
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
  bool dissociation = false;
  // bool h_escaped = false;

  runner.skmc(cell,
              omp_get_max_threads(),
              [&](double time,                        ///< Total time just after system at post
                  system::SoA<Position const &> pre,  ///< State just before mech applied
                  double E0,                    // Enegy of pre
                  int atom,                           ///< Index of central atom of mechanism
                  env::Mechanism const &mech,         ///< Chosen mechanism
                  system::SoA<Position const &> post,  ///< Final state of system after this iteration / mech
                  double Ef                     // energy of post
              ) {
                
                  
                run_time = time;

                fly::system::SoA<TypeID const &, Position const &> tmp(cell.size());

                tmp.rebind(r_, post);
                tmp.rebind(id_, cell);

                std::vector vac = detect.detect_vacancies(tmp);

                fmt::print("Found {} vacancies @{:::.2f}\n", vac.size(), vac);


                // V-dis testing by minimum spanning tree (MST)

                double const VV = kruskal_max(vac, min_image);

                fmt::print("MST max V-V = {:.3e}\n", VV);


                if (cluster.forVdiss){
                  dissociation = VV > cluster.diss_criterion;  // Vacancy-cluster dissociation criterion
                }
                else {
                  // H-escape testing
                  if (auto last = post.size() - 1; cell(id_, last) == 1) {
                    double vh = std::numeric_limits<double>::max();

                    for (auto const &v : vac) {
                      vh = std::min(vh, min_image(post(r_, last), v));
                    }

                    fmt::print("Min V-H = {:.3e}\n", vh);



                    dissociation = vh > cluster.diss_criterion;  // H-escape criterion set here.
                }
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

                fmt::print("Just wrote frame index No. {} (iteration no. {}, {} K) \n\n", 
                            file.n_frames() - 1, 
                            loopcounter + 1,          // a counter to print the current no. of iteration
                            temp);

                return dissociation;
              });

  // fmt::print("It took {:.3e}s for H to detrap\n", d_time);
  // fmt::print("It took {:.3e}s for cluster to dissociate\n", d_time);

  return run_time;
}




int main(){

  //-------------------------------- Type of cluster -------------------------------//
  std::string clustertype = V3H_Vdis.name;
    // Use: (+.name)
    //    V1H_Hdis; 
    //    V2H_Hdis, V2H_Vdis, V2_Vdis; 
    //    V3H_Vids; V3H_Hdis, V3_Vdis;

  std::string gsddirectory = "build/gsd/sim_output/" + clustertype + "_" + findsystemdate();   // gsd directory named after date of simulation
  std::filesystem::create_directory(gsddirectory);                                            


  int no_iteration = 5;                                          // Setting no. of iterations for each temperature
  std::vector<double> temperature;                                 // Setting the range of temperatures  
  for (int i=3; i<8; i++){
    temperature.push_back(i * 100.0);
  }


  for (int i=0; i< temperature.size(); ++i){  
    std::string datadirectory = "/home/zichuan/openFLY/build/data/";
    std::filesystem::create_directory(datadirectory + clustertype); 
    std::string taufiletemp = std::to_string(int(temperature[i]));                   // Files named after temperature for txt output

    std::ofstream taufile (datadirectory + clustertype + "/" + taufiletemp + "K.txt", std::ios::app);     // using append mode

    TempVdis diss_data {int(temperature[i]), no_iteration, 0};                       // Using struct defined above to count v-dis at each temperature

    for (int k=0; k<no_iteration; ++k){   


    // ------------------- Calculate lifetime --------------------- //
      double tauvalue = complex_lifetime(temperature[i], V3H_Vdis, gsddirectory, k);       


      // if (tauvalue == 0){                     // If dissociated - don't write anything
      //   diss_data.vcount++;
      // } 
      // else {
        taufile << tauvalue << "\n";          // Write lifetime to txt file
    //   }
    }

    std::ifstream countinfile (datadirectory + clustertype + "/dissociation_count.txt");
    int a, b, c, counter;
      a = b = c = counter = 0;
    std::vector<int> transfer_temp;
    std::vector<int> transfer_totalcount;
    std::vector<int> transfer_vcount;

    while (countinfile >> a >> b >> c ){
      if(int(temperature[i]) == a){                                     // detect if the same temperature exists
        transfer_temp.push_back(a); 
        transfer_totalcount.push_back(b + diss_data.totalcount);        // adding counts to existing data
        transfer_vcount.push_back(c + diss_data.vcount);
        counter++;                                                      
      }
      else{
        transfer_temp.push_back(a);                                     // transfer existing data
        transfer_totalcount.push_back(b);
        transfer_vcount.push_back(c);
      }
    }

    if(counter == 0){                                                   // Activated when there's no matching temp
        transfer_temp.push_back(diss_data.temp);
        transfer_totalcount.push_back(diss_data.totalcount);
        transfer_vcount.push_back(diss_data.vcount);
    }

    countinfile.close();


    std::ofstream countoutfile (datadirectory + clustertype + "/dissociation_count.txt"); 
      for (int j=0; j<transfer_temp.size(); j++){
        countoutfile << transfer_temp[j] << " " << transfer_totalcount[j] << " " << transfer_vcount[j] << "\n";
      }
  }
  

  return 0;
}

//------------------------- KEY PAEAMETERS -----------------------//
//  Deleted atoms / adding H or not
//  dissociation criterion: v-dis or h-dis
//  clustertype string name
//  bin file
//------------------------- Check parameters every time running the code -------------------------------//
