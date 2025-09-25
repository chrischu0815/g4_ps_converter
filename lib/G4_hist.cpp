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

#include "G4_hist.hpp"
#include <algorithm>
G4_hist::G4_hist(std::filesystem::path path_to_file_to_write,
                 std::string hist_type, std::string name,
                 std::string description, std::string comments, int Z,
                 std::vector<double> right_bin_edge,
                 std::vector<double> bin_count, double relative_intensity)
    : path_to_file_to_write(path_to_file_to_write), hist_type(hist_type),
      name(name), description(description), comments(comments), Z(Z),
      right_bin_edge(right_bin_edge), bin_count(bin_count),
      relative_intensity(relative_intensity) {}

int G4_hist::print_data() {
  for (size_t i = 0; i < this->right_bin_edge.size(); i++) {
    std::cout << i << ", " << right_bin_edge[i] << ", " << bin_count[i]
              << std::endl;
  }
  return 1;
}

int G4_hist::save_data_as_csv() {
  std::ofstream hist_geant(path_to_file_to_write);
  for (size_t i = 0; i < this->right_bin_edge.size(); i++) {
    hist_geant << this->right_bin_edge[i] << ", " << this->bin_count[i]
               << std::endl;
  }
  return 1;
}

int G4_hist::write_g4_ps_file() {
  double imass = 1;

  std::ofstream file_to_write(this->path_to_file_to_write, std::ofstream::out);

  file_to_write << "# " << this->name << std::endl;
  file_to_write << "# " << this->description << std::endl;
  file_to_write << "# " << this->comments << std::endl;

  file_to_write << "##########################" << std::endl;

  if (this->hist_type == "Ek") {

    file_to_write << "/gps/source/add " << this->relative_intensity
                  << std::endl;
    if (this->Z == 0) {
      file_to_write << "/gps/particle neutron" << std::endl;
    }

    else if (this->Z == 1) {
      file_to_write << "/gps/particle proton" << std::endl;
    } else {
      // TODO: use G4IonTable instead
      imass = G4NucleiProperties::GetNuclearMass(associated_A[Z], Z) /
              CLHEP::amu_c2;
      file_to_write << "/gps/particle ion" << std::endl;
      file_to_write << "/gps/ion " << Z << " " << this->associated_A[Z]
                    << " # with u = " << imass << std::endl;
    }

    file_to_write << "/gps/pos/type Plane" << std::endl;
    file_to_write << "/gps/pos/shape Square" << std::endl;
    file_to_write << "/gps/pos/centre 0 0 -" + this->beam_offset + " mm"
                  << std::endl; // FIXME: make this sign-depentend
    file_to_write << "/gps/direction 0 0 1" << std::endl;

    file_to_write << "/gps/pos/halfx " << this->dimHalfX << " mm" << std::endl;
    file_to_write << "/gps/pos/halfy " << this->dimHalfY << " mm" << std::endl;
    file_to_write << "/gps/ene/type Arb" << std::endl;
    file_to_write << "/gps/ene/diffspec 1" << std::endl;
    file_to_write << "/gps/hist/type arb" << std::endl;

    const double normalization =
        std::accumulate(this->bin_count.begin(), this->bin_count.end(), 0.);

    for (size_t i = 0; i < this->right_bin_edge.size(); i++) {

      file_to_write << "/gps/hist/point " << this->right_bin_edge[i] * imass
                    << " "
                    << (double)this->bin_count[i] / std::abs(normalization)
                    << std::endl;
    }
    file_to_write << "/gps/hist/inter Lin" << std::endl;

    file_to_write << "# EOF " << std::endl;

    file_to_write.close();
    return 0;

  } else if (this->hist_type == "X") {

    file_to_write << "/gps/hist/type biasx" << std::endl;

    const double norm = std::abs(
        *std::max_element(this->bin_count.begin(), this->bin_count.end()));

    const double offset = this->dimHalfX;

    const double max_element = 2 * this->dimHalfX;
    for (size_t i = 0; i < this->right_bin_edge.size(); i++) {
      file_to_write << "/gps/hist/point "
                    << (this->right_bin_edge[i] + offset) / max_element << " "
                    << this->bin_count[i] / norm << std::endl;
    }

    file_to_write << "# EOF " << std::endl;

    file_to_write.close();
    return 0;

  } else if (this->hist_type == "Y") {

    file_to_write << "/gps/hist/type biasy" << std::endl;

    const double norm = std::abs(
        *std::max_element(this->bin_count.begin(), this->bin_count.end()));

    const double offset = this->dimHalfY;

    const double max_element = 2 * this->dimHalfY;

    for (size_t i = 0; i < this->right_bin_edge.size(); i++) {
      file_to_write << "/gps/hist/point "
                    << (this->right_bin_edge[i] + offset) / max_element << " "
                    << this->bin_count[i] / norm << std::endl;
    }

    file_to_write << "# EOF " << std::endl;

    file_to_write.close();
    return 0;

  } else if (this->hist_type == "Theta") {

    file_to_write << "/gps/ang/type user" << std::endl;
    file_to_write << "/gps/hist/type theta" << std::endl;

    const double norm =
        std::accumulate(this->bin_count.begin(), this->bin_count.end(), 0.);
#ifdef DEBUG2
    std::cout << "Theta normalization = " << norm << std::endl;
#endif
    for (size_t i = 0; i < this->right_bin_edge.size(); i++) {
      double actual_value = (double)this->bin_count[i] / norm;
      if (actual_value < 1e-10) {
        actual_value = 0;
      }

      file_to_write << "/gps/hist/point " << this->right_bin_edge[i] << " "
                    << actual_value << std::endl;
    }

    file_to_write << "# EOF " << std::endl;

    file_to_write.close();
    return 0;

  } else if (this->hist_type == "Phi") {

    // file_to_write << "/gps/ang/type user" << std::endl;
    file_to_write << "/gps/hist/type phi" << std::endl;
    const double norm =
        std::accumulate(this->bin_count.begin(), this->bin_count.end(), 0.);
    for (size_t i = 0; i < this->right_bin_edge.size(); i++) {
      double actual_value = (double)this->bin_count[i] / norm;
      if (actual_value < 1e-10) {
        actual_value = 0;
      }
      file_to_write << "/gps/hist/point " << this->right_bin_edge[i] << " "
                    << actual_value << std::endl;
    }

    file_to_write << "# EOF " << std::endl;

    file_to_write.close();
    return 0;
  }

  return 0;
}

void CommandGenerator::generateCommandFile(const std::string &commandFileName,
                                           const double dimensions,
                                           const std::string &configuration,
                                           const std::vector<int> available_Z,
                                           // int maxZ,
                                           const std::vector<double> counts) {
  // Open the file for writing
  std::ofstream commandFile(commandFileName);

  // Check if the file is successfully opened
  if (!commandFile.is_open()) {
    std::cerr << "Error opening the command file." << std::endl;
    return;
  }
  commandFile << "###########################################" << std::endl;
  commandFile << "# Source this file to load the macro files generated "
              << std::endl;
  if (configuration == "100") {
    commandFile << "# File generated with configuration = " << configuration
                << "\n#equivalent to right after the the last element of the "
                   "beamline.\n"
                << "# Use this if your simulatsions requires a lot of space"
                << std::endl;
  } else if (configuration == "1500") {
    commandFile
        << "# File generated with configuration = " << configuration
        << "\n#equivalent to scoring the PS at the end of the beamline.\n"
        << "# Use this if your simulations requires small space" << std::endl;
  } else if (configuration == "2000") {
    commandFile
        << "# File generated with configuration = " << configuration
        << "\n#equivalent to scoring the PS at the end of the beamline.\n"
        << "# Use this if your simulations requires small space" << std::endl;
  } else {
    commandFile << "# File generated with configuration = NA!"
                << "\nThis is an error in the software configuration.\n"
                << std::endl;
  }
  commandFile << "# scoring dimensions = " << dimensions << " mm" << std::endl;
  commandFile << "# This file is automatically generated." << std::endl;
  commandFile << "###########################################" << std::endl;

  // commandFile << "/run/initialize" << std::endl;

  // Write the commands to the file
  // for (int Z = 0; Z <= maxZ; ++Z) {
  for (auto Z : available_Z) {
    // Iterate through the characteristics: Ek, X, Y, Theta, Phi
    // FIXME: make this dynamic
    for (const std::string &feature : {"Ek", "X", "Y", "Theta", "Phi"}) {
      // Generate the path for the .mac file for the current characteristic and
      // Z
      std::ostringstream filePath;
      // XXX: this should account the minium number of points. (required
      // points).
      filePath << "hist_geant_" << feature << "_" + configuration << "_Z_" << Z
               << ".mac";

      // Write the command to the file
      commandFile << "/control/execute " << filePath.str() << std::endl;
    }
  }

  commandFile << "###########################################" << std::endl;
  commandFile << "# EOF" << std::endl;

  // Close the file
  commandFile.close();

  std::cout << "Command file generated: " << commandFileName << std::endl;
}
