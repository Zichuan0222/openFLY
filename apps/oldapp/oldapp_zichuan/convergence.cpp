#include <fmt/core.h>
#include <math.h>
#include <cstdint>
#include <vector>
#include "omp.h"
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>


#include "libfly/utility/core.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/box.hpp"
#include "libfly/utility/lattice.hpp"
#include "libfly/io/gsd.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/potential/EAM/eam.hpp"


using namespace fly;                            // fly can be ignored from now



// Questions
  //  connecting to the department cluster / vpn
  //  plan for next two meetings


// This function will minimise its input (overwriting positions) and output the energy of the minimised system.
double min_and_cal_energy(system::SoA<Position &, TypeID const&, Frozen const&> in, minimise::LBFGS &minimiser, potential::Generic & pot, neigh::List & nl){
    
    system::SoA<Position, PotentialGradient> out(in.size());                                

    minimiser.minimise(out, in, pot, omp_get_max_threads());

    nl.rebuild(out, omp_get_max_threads());
    
    double min_E = pot.energy(in, nl, omp_get_max_threads());                        // output: minimised energy

    in[r_] = out[r_];

    return min_E;
}



// Input cell and positions (in armstrongs) of atoms to be deleted, output the indices of those atoms.
int position_to_index(system::SoA<Position const&> pos, Vec const& coor){
    
    for (int i=0; i<pos.size(); ++i) {
        Vec coori = pos(r_, i);
        double disti = gnorm (coori - coor);
        if (disti < 0.2){
            return i;
        }}

throw error("Could not find such atoms");
}


// positions of atoms to be deleted for v2 to v5
  std::vector<Vec> del_for_v2 = {
    {4.5, 4.5, 4.5},
    {4.5, 7.4, 4.5},
  };

  std::vector<Vec> del_for_v3 = {
    {4.5, 4.5, 4.5},
    {3.1, 6.0, 6.0},
    {6.0, 6.0, 6.0},
  };

  std::vector<Vec> del_for_v4 = {
    {4.5, 4.5, 4.5},
    {4.5, 4.5, 7.4},
    {3.1, 6.0, 6.0},
    {6.0, 6.0, 6.0},
  };

  std::vector<Vec> del_for_v5 = {
    {4.5, 4.5, 4.5},
    {4.5, 4.5, 7.4},
    {3.1, 6.0, 6.0},
    {4.5, 7.4, 7.4},
    {6.0, 6.0, 6.0},
  };



// This function exports two vectors to txt files in ./data
int exporttotxt(std::vector<double> exportvector1, std::vector<double> exportvector2, int const& nvac, std::vector<int> const& cellsize){

  std::ostringstream filenamenvac;
  filenamenvac << "/home/zichuan/openFLY/data/formation_energy_" << nvac << "V.txt";
  std::string filename = filenamenvac.str();
  std::ofstream MyFile(filename);
    MyFile.precision(15);


  for(int i; i < exportvector1.size(); ++i){
    MyFile << cellsize[i] << " " << exportvector1[i] << " " << exportvector2[i] << "\n";
  }

  return 0;
}


// input: number of vacancies, output Ek(N)
int formationenergy(int nvac)
{

  // Create a list of cluster formation energy
  std::vector<double> Ef_cluster;
  std::vector<double> Ef_clusterperfect;
  std::vector<int> cellsize;

for (int n = 14; n < 15; ++n) {

  system::TypeMap<> FeH(2);                     // type of atoms

  FeH.set(0, tp_, "Fe");
  FeH.set(1, tp_, "H");

  
  Mat basis{                            
      {2.855300, 0.000000, 0.000000},           // imput basis and initial lattice parameters
      {0.000000, 2.855300, 0.000000},
      {0.000000, 0.000000, 2.855300},
  };

  system::Supercell motif = system::make_supercell<Position, Frozen>({basis, Arr<bool>::Constant(true)}, FeH, 2);

  // cell contains position, typeID, and frozen elements


  motif[fzn_] = false;                                                  /// No frozen
  motif[id_] = 0;                                                       /// Iron set to motif

  motif(r_, 0) = Vec::Zero();                   
  motif(r_, 1) = Vec::Constant(0.5);



/// Set simulation scale
  int cubex = int((n+2)/3);                                             /// Take the integer part (round down) to achieve increase of 1 in one dimension each time
  int cubey = int((n+1)/3);
  int cubez = int(n/3);
      
    fmt::print("--------- supercell dimension: {}, {}, {}\n", cubex, cubey, cubez);                                                     



// Construct supercell
  system::Supercell cell = motif_to_lattice(motif, {cubex, cubey, cubez});     
    int natoms = cell.size();
    cellsize.push_back(natoms);

/// Create .gsd file for visualiser
  fly::io::BinaryFile file("build/test.gsd", fly::io::create);                      // To create a gsd file for output with "file name.gsd"
  
   file.commit([&] {                                                           
    file.write(cell.box());                                                         // Write the box to frame 0
    file.write(cell.map());                                                         // Write the map to frame 0
    file.write("particles/N", fly::safe_cast<std::uint32_t>(cell.size()));          // Write the number of atoms to frame 0
    file.write(id_, cell);                                                          // Write the typeID of atoms (id_) to frame 0
    file.write(r_, cell);                                                           // Write the position of atoms (r_) to frame 0 
  });



// Set minimiser with EAM potential
minimise::LBFGS minimiser({}, cell.box());                                         

  potential::Generic pot{                                                            
      potential::EAM{
          cell.map(),
          std::make_shared<potential::DataEAM>(std::ifstream{"data/wen.eam.fs"}),
      },
  };

  neigh::List nl(cell.box(), pot.r_cut());

nl.rebuild(cell, omp_get_max_threads());


/// Initial relaxation of supercell (no vacancy)
  double E0 = min_and_cal_energy(cell, minimiser, pot, nl);  

    file.commit([&] {      
      file.write("particles/N", fly::safe_cast<std::uint32_t>(cell.size()));          // Write the number of atoms to frame 0
      file.write(id_, cell);                                                          // Write the typeID of atoms (id_) to frame 0                                                     
      file.write(r_, cell);                                                           // Write the position of atoms (r_) to frame 0 
    });

    fmt::print("E0 = {} (no vacancy)\n", E0);



/// Create vacancies

  auto cell_v1 = remove_atoms(cell, {1});                                                  // atom 1 is always deleted 

  double E1 = min_and_cal_energy(cell_v1, minimiser, pot, nl);                             // lattice energy with a monovacancy at the position of atom 1
  fmt::print("E1 = {} (monovacancy)\n", E1);

      file.commit([&] {            
        file.write("particles/N", fly::safe_cast<std::uint32_t>(cell_v1.size()));          // Write the number of atoms to frame 0
        file.write(id_, cell_v1);                                                          // Write the typeID of atoms (id_) to frame 0                                               
        file.write(r_, cell_v1);                                                           // Write the position of atoms (r_) to frame 0 
    });

  


  std::vector<std::vector<Vec>> del_atoms;
    del_atoms = {del_for_v2, del_for_v3, del_for_v4, del_for_v5};
    

  std::vector<Eigen::Index> atom_to_remove;                                               // delete more atoms for cluster
  
    for (int i=0; i < nvac; ++i){
        atom_to_remove.push_back(position_to_index(cell, del_atoms[nvac-2][i]));
    }

    fmt::print("Removing:{}\n", atom_to_remove);




  auto cell_vk = timeit("remove", [&]{return remove_atoms(cell, atom_to_remove); });                                      // remove more atoms to form the cluster
  
  double Ek = timeit("min calc", [&]{return  min_and_cal_energy(cell_vk, minimiser, pot, nl);});                     // lattice energy with a k-vacancy cluster

  fmt::print("Ek = {} ({}-vacancy cluster)\n", Ek, nvac);

      file.commit([&] {            
        file.write("particles/N", fly::safe_cast<std::uint32_t>(cell_vk.size()));          // Write the number of atoms to frame 0
        file.write(id_, cell_vk);                                                          // Write the typeID of atoms (id_) to frame 0                                               
        file.write(r_, cell_vk);                                                           // Write the position of atoms (r_) to frame 0 
    });




  double Ef1 = (E1 - (natoms - 1) * E0 / natoms);                                                 // Ef for monovacancy
  double Efk = (Ek - (natoms - nvac) * E0 / natoms);                                              // Ef for k-vac cluster from perfect xtal
  double Ef = (Efk - nvac * Ef1);                                                                 // Ef for k-vac cluster from k independent vacs
  
    fmt::print("Ef1 = {}\n", Ef1);
    fmt::print("Efk = {}\n", Efk);
    fmt::print("cluster Ef = {}\n", Ef);

// Make a list of lattice energies at different N values

    Ef_cluster.push_back(Ef);
    fmt::print("Ef = {}\n", Ef_cluster);

    Ef_clusterperfect.push_back(Efk);
    fmt::print("Efk = {}\n", Ef_clusterperfect);

    
} // end of for loop in main
  

// Export Ef, Efk to csv files
  exporttotxt(Ef_cluster, Ef_clusterperfect, nvac, cellsize);

  return 0;
}



int main(){

  for(int nvac = 2; nvac < 6; ++nvac){

  formationenergy(nvac);

  }
  return 0;
}

