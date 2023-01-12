#include <iostream>
#include <chrono>
#include <ctime>    



int main()
{

    auto systemtime =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    
    std::cout << std::ctime(&systemtime) << "\n";
    struct tm * now = std::localtime(&systemtime);

    char buffer[80];
    strftime(buffer, 80, "%Y-%m-%d_%H-%M-%S", now);
    
    std::cout << buffer << "\n";

    // auto year = 1900 + now->tm_year;
    // auto month = 1 + now->tm_mon;
    // auto day = now->tm_mday;
    // auto hhmmss = now->hh_mm_ss;    

    // std::cout << year << "\n";
    // std::cout << year << "\n";
    // std::cout << year << "\n";
    // std::cout << year << "\n";

    return 0;

}