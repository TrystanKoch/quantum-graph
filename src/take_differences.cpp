/*----------------------------------------------------------------------
   take_differences.cpp 
    - Finds differences between numbers in a data file

   Copyright (C) 2015  Trystan Koch

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------*/

#include <limits>
#include <numeric>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cstring>



////////////////////////////////////////////////////////////////////////
/// Prototypes

bool ParseInputForNormalization(int, char*[]);
void ReadInDataFile(std::vector<std::string>&, std::vector<double>&);
void PrintPreviousHeader(std::vector<std::string>);
void OutputDifferenceInformation(double, bool);
void PrintDifferences(std::vector<double>);





////////////////////////////////////////////////////////////////////////
/// Main Function

int main(int argc, char *argv[])
{
  bool normalization_flag = ParseInputForNormalization(argc, argv);

  std::vector<std::string> header;
  std::vector<double> data;

  ReadInDataFile(header, data);
 
  std::sort(data.begin(), data.end());
  
  std::vector<double> differences(data.size());
  
  data.erase(std::remove(data.begin(),data.end(),0.0),data.end());

  // Take the differences between roots using an iterator method.
  std::adjacent_difference(data.begin(), data.end(), 
                           differences.begin());
      
  differences.erase(differences.begin());
  

  double difference_average 
      = std::accumulate(differences.begin(), differences.end(), 0.0);

  difference_average = difference_average / differences.size();
  
  if (normalization_flag)
  {
    for (unsigned int i=0; i<differences.size(); i++)
    {
      differences[i] /= difference_average;
    }
  }

  // Make sure program puts whole double into the output files.
  // Use scientific notation and constant length for numbers.
  // Print "true" or "false" for booleans
  std::cout.precision( std::numeric_limits<double>::max_digits10 );
  std::cout << std::fixed << std::scientific << std::boolalpha;

  
  // Output what we received from the input file
  PrintPreviousHeader(header);

  // Add some Information we acquired here so that the differences file
  // is more easily understood when we come back to it.
  OutputDifferenceInformation(difference_average, normalization_flag);

  // Print what we found.
  PrintDifferences(differences);

  return 0;
}





////////////////////////////////////////////////////////////////////////
/// Helper Functions

// There may or may not be a commandline argument for this function. If
// there is an argument that indicates that we want a list of normalized
// differences, we normalize the differences so that the average 
// difference is one. Otherwise, we just return a raw list of 
// differences.
bool ParseInputForNormalization(int argc, char* argv[])
{
  if (argc > 1)
  {
    for (int i=1; i<argc ; i++)
    {
      if ( !std::strcmp(argv[i],"-n") 
           or !std::strcmp(argv[i],"--normalized")
           or !std::strcmp(argv[i],"--normed") )
      {
        return true;
      }
    }
  }

  return false;
}





// Reads the data in the file into a header, containing the input file's
// comments and headmatter, and a vecator of doubles, which we use later
// to calculate the differences. Note that this function takes two
// output variables.
void ReadInDataFile(std::vector<std::string>& header,
                    std::vector<double>& data)
{
  // Read the input file. Pass the header along to the output file and
  // then send the data to a vector for further processing.
  std::string line; 
  while (std::getline(std::cin, line))
  {
    // Check if there is a header.
    if (line[0]=='#')
    {
      header.push_back(line);
      continue;
    }
    // After header, read the root values and store them.
    data.push_back(std::stod(line));
  }
}





// Simple Macro function that just makes the main code easier to read.
void PrintPreviousHeader(std::vector<std::string> header)
{
  for (unsigned int i=0; i<header.size(); i++)
  {
    std::cout << header[i] << std::endl;
  }
}





// Passes information about the differences we took to the output
// file. The meta-information about which run was what is kept so that
// it's easier to see what a file corresponds to.
void OutputDifferenceInformation(double difference_average,
                                 bool normalization_flag)
{
  // Print new information to the header. This will be the third heading
  // begining with "# BEGIN" and ending with "#END" in the file.
  // Whether we normalized the data or not, pass on what the old spacing
  // was and then say if we normalized the spacings before placing them
  // in the file, so it's clear regardless of the real average.
  std::cout << "# BEGIN"                             << std::endl;
  std::cout << "#  Spacing Information"              << std::endl;
  std::cout << "#  "                                 << std::endl;
  std::cout << "#    Average Spacing: "              << std::endl;
  std::cout << "#      " << difference_average       << std::endl;
  std::cout << "#    Spacings Normalized?"           << std::endl;
  std::cout << "#      " << normalization_flag       << std::endl;
  std::cout << "#  "                                 << std::endl;
  std::cout << "# END"                               << std::endl;
}





// Simple Macro function that just makes the main code easier to read.
void PrintDifferences( std::vector<double> differences)
{
  for (unsigned int i=0; i<differences.size(); i++)
  {
    std::cout << differences[i] << std::endl;
  }
}



