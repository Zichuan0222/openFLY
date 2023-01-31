#include <iostream>
#include <chrono>
#include <ctime>    
#include <typeinfo>
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


// std::string findsystemtime()
// {
//     auto systemtime =
//       std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    
//     std::cout << std::ctime(&systemtime) << "\n";
//     struct tm * now = std::localtime(&systemtime);

//     char buffer[80];
//     strftime(buffer, 80, "%Y-%m-%d_%H-%M-%S", now);
    
//     std::cout << buffer << "\n";

//     std::cout << typeid(buffer).name() << "\n";

//     return buffer;
// }



// int  exporttotxt(std::string filename,
//                 std::vector<double> exportvector1,
//                 std::vector<double> exportvector2,
//                 std::vector<double> exportvector3)
// {
//   std::string outputfile = "/home/zichuan/openFLY/build/data/test/" + filename + ".txt";
//   std::ofstream myfile(outputfile);
//     myfile.precision(20);

//   for(int i=0; i < exportvector1.size(); ++i){
//     myfile << exportvector1[i] << " " << exportvector2[i] << " " << exportvector3[i] << "\n";
//   }

//   return 0;
// }


struct TempVdis {  
  int temp {};     //< Temperature of simulation
  int totalcount {};
  int vcount {};    //< The no. of times V-dis before H-dis
};


int main(){

  int no_iteration = 5;

  std::vector<double> temperature = {300.0, 400.0, 500.0};

  std::vector<double> tau300K = {1.0, 1.2}; 
  std::vector<double> tau400K = {3.14, 2.22};
  std::vector<double> tau500K = {4.4, 4.26}; 




  for(int i=0; i<temperature.size(); ++i){
    
    // std::string taufiletemp = std::to_string(int(temperature[i]));
    // std::string outputfile 
    //   = "/home/zichuan/openFLY/build/data/test/" 
    //     + taufiletemp + "K.txt";
    
    // std::ofstream myfile (outputfile, std::ios::app);
    // myfile.precision(20);

    std::ifstream infile ("/home/zichuan/openFLY/build/data/test/dissociation_count.txt");
    int a, b, c, counter;
      a = b = c = counter = 0;
    std::vector<int> transfer_temp;
    std::vector<int> transfer_totalcount;
    std::vector<int> transfer_vcount;
    
    TempVdis diss_data {int(temperature[i]), 100, 10};

    while (infile >> a >> b >> c ){
      if(int(temperature[i]) == a){
        transfer_temp.push_back(a); 
        transfer_totalcount.push_back(b + diss_data.totalcount);
        transfer_vcount.push_back(c + diss_data.vcount);
        counter++;
      }
      else{
        transfer_temp.push_back(a);
        transfer_totalcount.push_back(b);
        transfer_vcount.push_back(c);
      }
    }

    if(counter == 0){
        transfer_temp.push_back(diss_data.temp);
        transfer_totalcount.push_back(diss_data.totalcount);
        transfer_vcount.push_back(diss_data.vcount);
    }

    infile.close();

    std::ofstream outfile ("/home/zichuan/openFLY/build/data/test/dissociation_count.txt");
      for (int j=0; j<transfer_temp.size(); j++){
        outfile << transfer_temp[j] << " " << transfer_totalcount[j] << " " << transfer_vcount[j] << "\n";
      }

    // for (int k=0; k< no_iteration; ++k){
    //   myfile <<  tau300K[k]  << "\n";
    // }
  }




  // for(int i=0; i<10; ++i){
  //   V11.push_back(i*1.0);
  //   V22.push_back(i*1.5);
  //   V33.push_back(i*3.14);
  //     std::cout << "test " << i <<"\n";
  // }

  // std::cout << "test \n";

  // for(int i=0; i<V11.size(); ++i){
  //   std::cout << V11[i] << " \n";
  // }
 
  // exporttotxt(filename, V11, V22, V33);

  return 0;
}