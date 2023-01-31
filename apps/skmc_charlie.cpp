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



struct TempVdis {  
  int temp {};        //< Temperature of simulation
  int totalcount {};  //< Total no. of sim at this temp
  int vcount {};      //< The no. of times V-dis before H-dis
};


// // Input cell and positions (in armstrongs) of atoms to be deleted, output the indices of those atoms.
// int position_to_index(system::SoA<Position const&> pos, Vec const& coor){
    
//     for (int i=0; i<pos.size(); ++i) {
//         Vec coori = pos(r_, i);
//         double disti = gnorm (coori - coor);
//         if (disti < 0.2){
//             return i;
//         }}

// throw error("Could not find such atoms");
// }



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
                        std::string gsdfiledirectory,
                        int loopcounter){

  system::Supercell perfect = motif_to_lattice(bcc_iron_motif(), {6, 6, 6});

  DetectVacancies detect(4, perfect.box(), perfect);


  // ------------------- delete atoms ------------------- //
  system::Supercell cell = remove_atoms(perfect, {1, 3});               // remove atoms by indices
  
  // Indices of atoms to be removed:
  // 1V - 1
  // 2V - 1, 3
  // 3V

  Vec r_H = {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857 / 4 + 3.14};

  // ------------------- adding H ------------------- //
  cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen>(1, r_H, false)});       //  add H into lattice


  //   cell = add_atoms(
  //       cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857
  //       / 4 * 2 + 3.14}, false, 0)});

  //   cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {0.992116, 6.01736, 4.56979},
  //   false, 0)});

  //   cell(r_, 431) = Vec{4.5896, 4.5892, 5.91909};
  //   cell(r_, 0) = Vec{4.45211, 4.45172, 4.2526};


  // std::string gsdfilename                                 //<<<<<<<<<<<<<<<<< Use this line to save all .gsd file
  //     = gsdfiledirectory + "/" 
  //       + std::to_string(int(temp)) + "K_"            
  //       + findsystemtime() + ".gsd";                   

  std::string gsdfilename = gsdfiledirectory + "sim.gsd";    //<<<<<<<<<<<<<< Use this line to save only current .gsd file

  fly::io::BinaryFile file(gsdfilename, fly::io::create);

  auto vac = detect.detect_vacancies(cell);

  fmt::print("Found {} vacancies @{:::.2f}\n", vac.size(), vac);

  auto const N = vac.size();

  file.commit([&] {
    file.write(cell.box());
    file.write(cell.map());

    auto special = explicit_V(vac, cell);

    file.write("particles/N", fly::safe_cast<std::uint32_t>(special.size()));

    file.write("log/time", -1.);                // write -1 to initialise
    file.write("log/barrier", -1.);
    file.write("log/energy", -1.);

    file.write(id_, special);
    file.write(r_, special);
  });

  kinetic::SKMC runner = {
      {
          .debug = true,

          //---------------------------- bin file directory --------------------------//
          .fread = "build/gsd/cat.v2h.bin",      
            // Use a different bin file for with/without H
            // Names of bin files: v1h, v2, v2h,
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

  double d_time = 0;
  int count = 0;


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
                
                d_time = time;

                fly::system::SoA<TypeID const &, Position const &> tmp(cell.size());

                tmp.rebind(r_, post);
                tmp.rebind(id_, cell);

                auto v2 = detect.detect_vacancies(tmp);

                verify(v2.size() == N, "Num v changed");

                fmt::print("Found {} vacancies @{:::.2f}\n", v2.size(), v2);

                fmt::print("E0={:.8e}, Ef={:.8e}\n", E0, Ef);

                auto dist = distances(cell.box(), v2, post(r_, post.size() - 1));

                fmt::print("Max V-V = {:.3e}, min V-H = {:.3e}\n", dist.v_v, dist.v_h);

                auto vpost = explicit_V(v2, tmp);

                file.commit([&] {
                  auto copy = vpost;
                  copy[r_].head(pre[r_].size()) = pre[r_];
                  file.write(r_, copy);
                  file.write("log/time", -1.);
                  file.write("log/barrier", -1.);
                  file.write("log/energy", E0);
                });

                file.commit([&] { 
                  file.write(r_, vpost);  
                  file.write("log/time", time);
                  file.write("log/barrier", mech.barrier);
                  file.write("log/energy", Ef);
                });

                fmt::print("Just wrote frame index No. {} (iteration no. {}, {} K) \n", 
                            file.n_frames() - 1, 
                            loopcounter + 1,          // a counter to print the current no. of iteration
                            temp);



                
                //-------------------------- Dissociation criterion ----------------------------//

                //// ---- H-diss criterion ---- ////
                // if (dist.v_v > 4.85){           // detect if V diss occurs before H diss
                //   ++count;
                //   fmt::print("The cluster dissociates before H escapes after {:.3e}s \n", d_time);
                //   d_time = 0.0;
                // }
                
                // if (dist.v_h > 6){
                //   ++count;
                //   fmt::print("It took {:.3e}s for H to detrap\n", d_time);
                // }

                //// ---- V-diss criterion ----- ////
                if (dist.v_v > 4.85){
                  ++count;
                  fmt::print("The cluster dissociates after {:.3e}s \n", d_time);
                }

                return count >= 1;           // Hdiss 
                // return dist.v_v > 4.85;      // Vdiss (  5NN - 4.94A, 4NN - 4.74A)
              });

  // fmt::print("It took {:.3e}s for H to detrap\n", d_time);
  // fmt::print("It took {:.3e}s for cluster to dissociate\n", d_time);

  return d_time;
}



int main(){

  //-------------------------------- Type of cluster -------------------------------//
  std::string clustertype = "V2H_Vdis";
    
    // Use: V2H_Hdis, 

  std::string gsddirectory = "build/gsd/sim_output/" + clustertype + "_" + findsystemdate();   // gsd directory named after date of simulation
  std::filesystem::create_directory(gsddirectory);                                            


  // int countvdis = 0;                        // counting cases when v-dis occurs before h-diss
  // std::vector<double> countvdis
  // int counttotal = 0;                       

  // std::vector<double> tau_mean;             // Declare mean value for tau here, each entry is calculated from a number of iteration at a temperature
  // std::vector<double> tau_stdev;            // Same as above, std dev of tau


  int no_iteration = 4;                                            // Setting no. of iterations for each temperature
  std::vector<double> temperature;                                 // Setting the range of temperatures  
  for (int i=3; i<11; i++){
    temperature.push_back(i * 100.0);
  }


  for (int i=0; i< temperature.size(); ++i){  
    std::filesystem::create_directory("/home/zichuan/openFLY/build/data/" + clustertype); 
    std::string taufiletemp = std::to_string(int(temperature[i]));                   // Files named after temperature for txt output

    std::ofstream taufile ("/home/zichuan/openFLY/build/data/" + clustertype + "/" + taufiletemp + "K.txt", std::ios::app);     // using append mode

    TempVdis diss_data {int(temperature[i]), no_iteration, 0};                       // Using struct defined above to count v-dis at each temperature

    for (int k=0; k<no_iteration; ++k){   
      double tauvalue = complex_lifetime(temperature[i], gsddirectory, k);           // Call the above function to calculate lifetime

      if (tauvalue == 0){                     // If dissociated - don't write anything
        diss_data.vcount++;
      } 
      else {
        taufile << tauvalue << "\n";          // Write lifetime to txt file
      }
    }

    std::ifstream countinfile ("/home/zichuan/openFLY/build/data/" + clustertype + "/dissociation_count.txt");
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


    std::ofstream countoutfile ("/home/zichuan/openFLY/build/data/" + clustertype + "/dissociation_count.txt"); 
      for (int j=0; j<transfer_temp.size(); j++){
        countoutfile << transfer_temp[j] << " " << transfer_totalcount[j] << " " << transfer_vcount[j] << "\n";
      }
      
    // std::sort(tau.begin(), tau.end());                // Sort and print liftimes in ascending order
    // std::cout << "\nLiftimes: \n";
    // for (int k=0; k<tau.size(); ++k){
    //   std::cout << tau[k] << ", ";
    // }     

    // // Calculating mean and std dev of tau
    // double sum = std::accumulate(tau.begin(), tau.end(), 0.0);
    // double mean = sum / tau.size(); 

    // std::vector<double> squareterms(tau.size());
    // for (int k=0; k<tau.size(); ++k){
    //   squareterms[k] = ((tau[k] - mean)) * ((tau[k] - mean));                                              
    // }
    
    // double squaresum = std::accumulate(squareterms.begin(), squareterms.end(), 0.0);
    // double stdev = std::sqrt(squaresum / tau.size());                                        


    // tau_mean.push_back(mean);
    // tau_stdev.push_back(stdev);

  }
  

// // export data to txt
//   std::string outputdirectory = "/home/zichuan/openFLY/build/data/" + clustertype + "/";     //< export directory
//   std::filesystem::create_directory(outputdirectory);  
//   std::string outputfile                                                                     //< V: v-dissociation, H: h-dissociation
//     = outputdirectory                                               
//       + findsystemdate() + "_" + findsystemtime() + ".txt";                                  //< and date + time of simulation
//   exporttotxt(outputfile, temperature, tau_mean, tau_stdev);

  return 0;
}

//------------------------- KEY PAEAMETERS -----------------------//
//  Deleted atoms / adding H or not
//  dissociation criterion: v-dis or h-dis
//  clustertype string name
//  bin file
//------------------------- Check parameters every time running the code -------------------------------//
