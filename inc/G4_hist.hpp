// Copyright © 2024 GSI Helmholtzzentrum für Schwerionenforschung gmbh
// <https://gsi.de/>
//
// Everyone is permitted to copy and distribute verbatim copies of this license
// document, but changing it is not allowed.
//
// Preamble
// The GNU General Public License is a free, copyleft license for software and
// other kinds of works. The licenses for most software and other practical
// works are designed to take away your freedom to share and change the works.
// By contrast, the GNU General Public License is intended to guarantee your
// freedom to share and change all versions of a program--to make sure it
// remains free software for all its users. We, GSI Helmholtzzentrum für
// Schwerionenforschung, use the GNU General Public License for most of our
// software; it applies also to any other work released this way by its authors.
// You can apply it to your programs, too.

#pragma once

#include "G4NucleiProperties.hh"
#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <numeric>
#include <vector>
/* class containing the G4 histogram
this class can host a G4 histogram
 */

class G4_hist {
public:
  G4_hist(std::filesystem::path path_to_file_to_write, std::string hist_type,
          std::string name, std::string description, std::string comments,
          int Z, std::vector<double> right_bin_edge,
          std::vector<double> bin_count, double relative_intensity);

  /* Print the data in the class as csv
      if success: return 1
  */
  int print_data();

  /* Save the data of the class to file*/
  int save_data_as_csv();

  /* Write G4 phase space format*/
  /** @brief create a file by writing the data to a text file with the same as a
   * Macro.
   *
   * Write the class content to a file compatible with geant4,
   * the file contains the required information to create an histogram
   * as Phase Space in geant4
   * If the file exists, the file will be overwritten
   *
   * It creates a file with the name as first parameter,
   * and description as second parameter
   * Then it continues linke this:
   * // # name
   * // # description
   * // # comments
   *  ##########################
   * gps/hist/point	1.0023	202.26
   * gps/hist/point	1.0162	207.42
   * gps/hist/point	1.0304	212.75
   * gps/hist/point	1.0447	218.26
   * for each of the keys: Ek, X, Y, theta and phi.
   * creates
    @return 0 if no error
  */
  int write_g4_ps_file();

  std::map<int, int> associated_A = {
      {1, 1},    {2, 4},    {3, 7},
      {4, 9},    {5, 11},   {6, 12},
      {7, 14},   {8, 16},   {9, 19},
      {10, 20},  {11, 23},  {12, 24},
      {13, 27},  {14, 28},  {15, 31},
      {16, 32},  {17, 35},  {18, 40},
      {19, 39},  {20, 40},  {21, 45},
      {22, 48},  {23, 51},  {24, 52},
      {25, 55},  {26, 56},  {27, 59},
      {28, 59},  {29, 63},  {30, 65},
      {31, 70},  {32, 73},  {33, 75},
      {34, 79},  {35, 80},  {36, 84},
      {37, 85},  {38, 88},  {39, 89},
      {40, 91},  {41, 93},  {42, 96},
      {43, 79},  {44, 101}, {45, 103},
      {46, 106}, {47, 108}, {48, 112},
      {49, 115}, {50, 119}, {51, 123},
      {52, 128}, {53, 126}, {54, 132},
      {55, 133}, {56, 137}, {57, 139},
      {58, 140}, {59, 141}, {60, 144},
      {61, 145}, {62, 150}, {63, 152},
      {64, 157}, {65, 159}, {66, 162} /* This is wrong */,
      {67, 165}, {68, 167}, {69, 169},
      {70, 173}, {71, 175}, {72, 178},
      {73, 181}, {74, 184}, {75, 186},
      {76, 190}, {77, 192}, {78, 195},
      {79, 197}, {80, 201}};

  std::map<int, int> source_intensity = {
      {0, 1},  {1, 1},  {2, 1},  {3, 1},  {4, 1},  {5, 1},  {6, 1},  {7, 1},
      {8, 1},  {9, 1},  {10, 1}, {11, 1}, {12, 1}, {13, 1}, {14, 1}, {15, 1},
      {16, 1}, {17, 1}, {18, 1}, {19, 1}, {20, 1}, {21, 1}, {22, 1}, {23, 1},
      {24, 1}, {25, 1}, {26, 1}, {27, 1}, {28, 1}};

  double dimHalfX = 200; // mm
  double dimHalfY = 200; // mm
private:
  std::filesystem::path path_to_file_to_write;
  std::string hist_type;
  std::string name;
  std::string description;
  std::string comments;
  int Z;
  std::vector<double> right_bin_edge;
  std::vector<double> bin_count;
  double relative_intensity;

  const std::string beam_offset =
      "200"; // mm of offset from the origin
             // this shift the beams at -200 mm from the origin
             // in the generated phase space file

  // TODO: define getter and setter
  /* actual dimension of the beam generator
   * NOTE: must match with the scoring surface dimensions
   * */
};

/*
 *  @brief Generate the sim.mac file containing all the sourcing.
 * This fill will source all the macro files thas are generated.
 * By sourcing this file in geant all the geneate gps macro are included.
 * @params commandFileName: file name to open and write on
 * @params outputPath: no used
 * @params configuration: defines which configuration is used and
 * consequentially set up the file name of the sourced files
 * @param maxZ: maximum Z gerated in the macro. So the sourcing will continue up
 * to this defined Z.
 * @return void
 *
 * */
class CommandGenerator {
public:
  static void generateCommandFile(const std::string &commandFileName,
                                  const double outputPath,
                                  const std::string &configuration,
                                  const std::vector<int> Z,
                                  const std::vector<double> counts);
};
