#include "Exception.hpp"
#include <iostream>



Exception::Exception(std::string sum, std::string prob)
{
  problem = prob;
  summary = sum;
}



void Exception::DebugPrint()
{
  std::cerr << "**  Exception ("<<summary<<") **\n";
  std::cerr << "Problem: " << problem << "\n\n";
}
