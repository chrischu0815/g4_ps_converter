// Copyright © 2026 GSI Helmholtzzentrum für Schwerionenforschung gmbh
// <https://gsi.de/>
//
// Everyone is permitted to copy and distribute verbatim copies of this license
// document, but changing it is not allowed.
//
// Preamble
// The GNU Affero General Public License (version 3) is a free, copyleft
// license for software and other kinds of works. The licenses for most
// software and other practical works are designed to take away your freedom
// to share and change the works. By contrast, the GNU Affero General Public
// License is intended to guarantee your freedom to share and change all
// versions of a program--to make sure it remains free software for all its
// users. We, GSI Helmholtzzentrum für Schwerionenforschung, use the GNU
// Affero General Public License for most of our software; it applies also to
// any other work released this way by its authors. You can apply it to your
// programs, too.

/**
 * @file G4_PS_converter.cpp
 *
 * @brief Main program of the G4 phase space converter (G4_PS_creator).
 *
 * Parses the command line, instantiates G4_PS_hist to convert the input
 * ROOT file into Geant4 GPS macro files, and finally generates the main
 * sim.mac file with CommandGenerator::generateCommandFile().
 *
 * @author enrico.pierobon
 * Contact: e.pierobon@gsi.de
 *
 */

#include "G4_PS_hist.hpp"
#include "G4_hist.hpp"
#include <iostream>
#include <stdexcept>

/**
 * @brief Print the command line usage and the available selections.
 *
 * The configuration and dimension selections are 1-based indices into
 * G4_PS_hist::AVAILABLE_CONFIGURATION and G4_PS_hist::AVAILABLE_DIMENSIONS.
 *
 * @return 1 (the program exits with an error status when usage is printed)
 * */
int printUsage() {
  std::cout << "Name" << std::endl;
  std::cout << "\tG4_PS_creator - create a phase space file with macros for "
               "geant4 given an input root file."
            << std::endl;
  std::cout << "SYNOPSIS" << std::endl;
  std::cout << "\tG4_PS_creator [input] [output] [configuration] [dimensions]"
            << std::endl;
  std::cout << "DESCRIPTION" << std::endl;
  std::cout << "\tCreate multiple .mac files. The macro file are to be used in "
               "a geant4 simulation with GPS. "
            << std::endl;
  std::cout << "\t[input]" << std::endl;
  std::cout << "\t\tinput path of the root file. The root file should contain "
               "the following keys:"
            << std::endl;
  std::cout << "\t\t\tEkinU : kinetic energy distribution in MeV/u"
            << std::endl;
  std::cout << "\t\t\tX: position of the x coordinate in mm" << std::endl;
  std::cout << "\t\t\tY: position of the y coordinate in mm" << std::endl;
  std::cout << "\t\t\t\t X and Y indentify the plane orthogonal to the beam "
               "direction (z)"
            << std::endl;
  std::cout << "\t\t\tTheta: theta coordinate in spherical coordinates "
               "(radians)"
            << std::endl;
  std::cout << "\t\t\tPhi: phi coordinate in spherical coordinates (radians)"
            << std::endl;
  std::cout << "\t[output]" << std::endl;
  std::cout << "\t\toutput path for the macro that will be generated. The "
               "simulation has to be run from this directory, using sim.mac "
               "as the Geant4 macro."
            << std::endl;

  std::cout << "\t[options]" << std::endl;
  std::cout << "\t\tTwo integers selecting the user configuration and the "
               "user dimensions respectively."
            << std::endl;
  std::cout << "AVAILABLE SELECTIONS" << std::endl;
  std::cout << "\tconfiguration:";
  for (size_t i = 0; i < G4_PS_hist::AVAILABLE_CONFIGURATION.size(); i++) {
    std::cout << " " << (i + 1) << "="
              << G4_PS_hist::AVAILABLE_CONFIGURATION[i];
  }
  std::cout << std::endl;
  std::cout << "\tdimensions:";
  for (size_t i = 0; i < G4_PS_hist::AVAILABLE_DIMENSIONS.size(); i++) {
    std::cout << " " << (i + 1) << "=\"" << G4_PS_hist::AVAILABLE_DIMENSIONS[i]
              << "\"";
  }
  std::cout << std::endl;
  std::cout << "USAGE EXAMPLE" << std::endl;
  std::cout << "\t./G4_PS_converter ../root_files/join.root /tmp/" << std::endl;
  std::cout << "\t./G4_PS_converter ../root_files/join.root /tmp/ 4 1  # "
               "configuration \"1\", normal scoring plane"
            << std::endl;
  return 1;
}

/**
 * @brief Program entry point.
 *
 * Usage:
 * @code
 * G4_PS_creator [input] [output] [configuration] [dimensions]
 * @endcode
 *
 *  - [input]: path of the ROOT file containing the phase space histograms
 *    (Ekin in MeV/u, X/Y in mm, Theta/Phi in radians).
 *  - [output]: directory where the generated GPS macros and sim.mac are
 *    written (created if missing).
 *  - [configuration] / [dimensions]: optional 1-based selections into
 *    G4_PS_hist::AVAILABLE_CONFIGURATION / G4_PS_hist::AVAILABLE_DIMENSIONS.
 *    If omitted, the first entry of each list (configuration "100", no
 *    dimension) is used.
 *
 * The program builds a G4_PS_hist, runs process_file() and finally
 * CommandGenerator::generateCommandFile() to produce sim.mac.
 *
 * @param argc argument count
 * @param argv argument values
 * @return 0 on success, 1 on any error (bad usage, invalid selection, or a
 * failure in the conversion)
 * */
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
  int user_configuration = 0;
  int user_dimensions = 0;
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
    } catch (const std::exception &e) {
      std::cerr << "Error: cannot convert user configuration \"" << argv[3]
                << "\": " << e.what() << std::endl;
      printUsage();
      return 1;
    }
    try {
      user_dimensions = std::stoi(argv[4]);
    } catch (const std::exception &e) {
      std::cerr << "Error: cannot convert user dimensions \"" << argv[4]
                << "\": " << e.what() << std::endl;
      printUsage();
      return 1;
    }

  } else {
    std::cerr << "Error: incorrect parameter usage detected" << std::endl;
    printUsage();
    return 1;
  }

  // validate the user selections before indexing the available vectors
  // (user selections are 1-based indices)
  if (user_configuration < 0 ||
      user_configuration >
          static_cast<int>(G4_PS_hist::AVAILABLE_CONFIGURATION.size() - 1) ||
      user_dimensions < 0 ||
      user_dimensions >
          static_cast<int>(G4_PS_hist::AVAILABLE_DIMENSIONS.size()) - 1) {
    std::cerr << "Error: invalid user configuration/dimensions: "
              << user_configuration << "/" << user_dimensions << std::endl;
    std::cerr << "Available configurations:";
    for (size_t i = 0; i < G4_PS_hist::AVAILABLE_CONFIGURATION.size(); i++) {
      std::cerr << " " << i << "=" << G4_PS_hist::AVAILABLE_CONFIGURATION[i];
    }
    std::cerr << std::endl;
    std::cerr << "Available dimensions:";
    for (size_t i = 0; i < G4_PS_hist::AVAILABLE_DIMENSIONS.size(); i++) {
      std::cerr << " " << i << "=\"" << G4_PS_hist::AVAILABLE_DIMENSIONS[i]
                << "\"";
    }
    std::cerr << std::endl;
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
  std::vector<std::string> available_particles;

  G4_PS_hist H1(inputFile, outputPath, configuration, dimensions);

  try {
    available_particles = H1.process_file();
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  if (available_particles.empty() || available_particles[0] == "-1") {
    std::cerr << "Error: no valid data found in the input file for the "
                 "selected configuration/dimensions."
              << std::endl;
    return 1;
  }
#ifdef DEBUG
  std::cout << "Got an array of processed Z = " << std::endl;
  for (auto a : available_particles) {
    std::cout << a << " ,";
  }
  std::cout << std::endl;
#endif
  try {
    CommandGenerator::generateCommandFile(
        outputFile, G4_PS_hist::SOURCE_DIMENSIONS[user_dimensions],
        configuration, available_particles);
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
