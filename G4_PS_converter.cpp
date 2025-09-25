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

/**
 * @file G4_PS_converter.cc
 *
 * @brief main program to create PS file
 *
 *
 * @author enrico.pierobon
 * Contact: e.pierobon@gsi.de
 *
 */

#include "G4_PS_hist.hpp"
#include "G4_hist.hpp"
#include <iostream>
#include <stdexcept>

int printUsage() {
  std::cout << "Name" << std::endl;
  std::cout << "\tG4_PS_creator - create a phase space file with macros for "
               "geant4 given an input root file."
            << std::endl;
  std::cout << "SYNOPSIS" << std::endl;
  std::cout << "\tG4_PS_creator [input] [output] [options]" << std::endl;
  std::cout << "DESCRIPTION" << std::endl;
  std::cout << "\tCreate multiple .mac files. The macro file are to be used in "
               "a geant4 simulation with GPS. "
            << std::endl;
  std::cout << "\t[input]" << std::endl;
  std::cout << "\t\tinput path of the root file. The root file should contain "
               "the following keys:"
            << std::endl;
  std::cout << "\t\t\tEkinU : kietick energy distribution in MeV/u"
            << std::endl;
  std::cout << "\t\t\tX: position of the x coordinate in mm" << std::endl;
  std::cout << "\t\t\tY: position of the y coordinate in mm" << std::endl;
  std::cout << "\t\t\t\t X and Y indentify the plane orthogonal to the beam "
               "direction (z)"
            << std::endl;
  std::cout << "\t\t\tTheta: theta coordinate in spherical coordinates"
            << std::endl;
  std::cout << "\t\t\tPhi: phi coordinate in spherical coordinates"
            << std::endl;
  std::cout << "\t[output]" << std::endl;
  std::cout << "\t\toutput path for the macro that will be generated."
            << std::endl;

  std::cout << "\t[options]" << std::endl;
  std::cout << "\t\tTwo integers selectiong the user configuraiton and the "
               "user dimensions respectivly."
            << std::endl;
  std::cout << "USAGE EXAMPLE" << std::endl;
  std::cout << "\t./G4_PS_converter ../root_files/join.root /tmp/" << std::endl;
  std::cout << "\t./G4_PS_converter ../root_files/join.root /tmp/ 1 2"
            << std::endl;
  return 1;
}

int main(int argc, char *argv[]) {
  /*   G4_PS_hist H1(1024, // Npoints
                  1,    // min value
                  10,   // max value
                  "/home/questo/data/sim/G4_phase_space_example/"
                  "gcr_p_10gcm2.root",  // path to root input file
                  "/tmp/schifezza.mac", //
                  "Title!",             //
                  "Description long "   //
    ); */

  std::string inputFile;
  std::string outputFile;
  std::string outputPath;
  uint user_configuration = 0;
  uint user_dimensions = 0;
  if (argc == 3) {
    /*
     * two arguments were passed: input and output file, process them instead
     */
    inputFile = argv[1];
    outputPath = argv[2];
    outputFile = outputPath + "/sim.mac";
    user_configuration = 1;
    user_dimensions = 1;

  } else if (argc == 5) {
    inputFile = argv[1];
    outputPath = argv[2];
    outputFile = outputPath + "/sim.mac";
    try {
      user_configuration = std::stoi(argv[3]);
    } catch (std::invalid_argument e) {
      std::cerr << "Cannot convert user configuration " << e.what()
                << std::endl;
    }
    try {
      user_dimensions = std::stoi(argv[4]);
    } catch (std::invalid_argument e) {
      std::cerr << "Cannot convert user dimensions " << e.what() << std::endl;
    }

  } else {
    std::cerr << "Error: incorrect parameter usage detected" << std::endl;
    printUsage();
    inputFile = "../root_files/join.root";
    outputFile = "/tmp/sim.mac";
    outputPath = "/tmp";
    return 1;
  }

#ifdef DEBUG
  std::cout << "Warning: running with debug flags" << std::endl;
#endif

  std::string configuration =
      G4_PS_hist::AVAILABLE_CONFIGURATION[user_configuration];
  std::string dimensions = G4_PS_hist::AVAILABLE_DIMENSIONS[user_dimensions];
  std::cout << "configuration = " << configuration << std::endl;
  std::cout << "dimensions = " << dimensions << std::endl;
  std::vector<int> availableZ;

  G4_PS_hist H1(inputFile, outputPath, configuration, dimensions);

  availableZ = H1.process_file();

  if (availableZ[0] == -1) {
    std::cerr << "Error: no valid input provided" << std::endl;
    return 1;
  }
#ifdef DEBUG
  std::cout << "Got an array of processed Z = " << std::endl;
  for (auto a : availableZ) {
    std::cout << a << " ,";
  }
  std::cout << std::endl;
#endif
  CommandGenerator::generateCommandFile(
      outputFile, G4_PS_hist::SOURCE_DIMENSIONS[user_dimensions], configuration,
      availableZ, H1.relative_abundance.second);

  return 0;
}
