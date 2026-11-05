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

#include "G4_hist.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

/**
 * @brief Constructor; see the header for the full documentation.
 *
 * Stores the output location, the histogram data and the metadata, and
 * parses the particle identifier as a charge Z when possible (fZ = -1 for
 * special particle names).
 * */
G4_hist::G4_hist(std::filesystem::path path_to_file_to_write,
                 std::string hist_type, std::string name,
                 std::string description, std::string comments,
                 std::string particle, std::vector<double> right_bin_edge,
                 std::vector<double> bin_count, double relative_intensity,
                 bool is_first_source)
    : fpath_to_file_to_write(path_to_file_to_write), fhist_type(hist_type),
      fname(name), fdescription(description), fcomments(comments),
      fparticle(particle), fright_bin_edge(right_bin_edge),
      fbin_count(bin_count), frelative_intensity(relative_intensity),
      fIsFirstSource(is_first_source) {
  /**
   * @brief parse fZ if possible
   * */
  try {
    fZ = std::stoi(particle);
  } catch (const std::invalid_argument) {
    fZ = -1;
  } catch (const std::out_of_range) {
    // not a charge number: treat it like a special particle name
    fZ = -1;
  }
}

/**
 * @brief Print the histogram data (index, bin center, value) on stdout.
 *
 * @return 1 on success
 * */
int G4_hist::print_data() {
  for (size_t i = 0; i < this->fbin_count.size(); i++) {
    const double center =
        0.5 * (this->fright_bin_edge[i] + this->fright_bin_edge[i + 1]);
    std::cout << i << ", " << center << ", " << this->fbin_count[i]
              << std::endl;
  }
  return 1;
}

/**
 * @brief Debug helper: save (bin center, value) pairs to the output file as
 * CSV.
 *
 * @return 1 on success
 * */
int G4_hist::save_data_as_csv() {
  std::ofstream hist_geant(fpath_to_file_to_write);
  for (size_t i = 0; i < this->fbin_count.size(); i++) {
    const double center =
        0.5 * (this->fright_bin_edge[i] + this->fright_bin_edge[i + 1]);
    hist_geant << center << ", " << this->fbin_count[i] << std::endl;
  }
  return 1;
}

/**
 * @brief Write the Geant4 GPS macro file for this histogram.
 *
 * Emits, according to the histogram type: the source intensity and the
 * particle/source-plane definition plus the arbitrary energy distribution
 * (Ek), the biased transverse position distribution (X/Y), or the user
 * angular distribution (Theta/Phi). See the header for the per-quantity
 * details.
 *
 * @return 0 on success, 1 if the file cannot be opened
 * @throw std::runtime_error if the histogram is empty (undefined
 * normalization) or no isotope is defined for the particle
 * */
int G4_hist::write_g4_ps_file() {
  std::ofstream file_to_write(this->fpath_to_file_to_write, std::ofstream::out);
  if (!file_to_write.is_open()) {
    std::cerr << "Error: cannot open " << this->fpath_to_file_to_write
              << " for writing (G4_hist::write_g4_ps_file)" << std::endl;
    return 1;
  }

  file_to_write << "# " << this->fname << std::endl;
  file_to_write << "# " << this->fdescription << std::endl;
  file_to_write << "# " << this->fcomments << std::endl;

  file_to_write << "##########################" << std::endl;

  if (this->fhist_type == "Ek") {

    // Mass number A of the particle: needed to convert the Ekin stored in
    // MeV/u to the total kinetic energy per particle in GeV that GPS expects.
    double imass = 1.0;
    if (this->fZ >= 1) {
      const auto a_it = this->associated_A.find(this->fZ);
      if (a_it == this->associated_A.end()) {
        throw std::runtime_error(
            "No isotope defined for Z = " + std::to_string(this->fZ) +
            " (G4_hist::write_g4_ps_file)");
      }
      // imass = static_cast<double>(a_it->second);
      imass = G4NucleiProperties::GetNuclearMass(a_it->second, this->fZ) /
              CLHEP::amu_c2;
    }

    // Geant4 GPS always starts with a default source 0 (a 1 MeV geantino with
    // intensity 1.0). The first particle configures that source directly;
    // every other particle creates an additional source. Otherwise the
    // geantino source 0 would keep 50 % of all generated primaries, since
    // the relative abundances sum to 1.
    file_to_write << (this->fIsFirstSource ? "/gps/source/intensity "
                                           : "/gps/source/add ")
                  << this->frelative_intensity << std::endl;
    this->define_source(file_to_write, fparticle, fZ);

    file_to_write << "/gps/pos/type Plane" << std::endl;
    file_to_write << "/gps/pos/shape Square" << std::endl;
    file_to_write << "/gps/pos/centre 0 0 -" + this->beam_offset + " mm"
                  << std::endl; // FIXME: make this sign-depentend
    file_to_write << "/gps/direction 0 0 1" << std::endl;

    file_to_write << "/gps/pos/halfx " << this->dimHalfX << " mm" << std::endl;
    file_to_write << "/gps/pos/halfy " << this->dimHalfY << " mm" << std::endl;

    // Deterministic angular state for this source: if no Theta/Phi macro
    // follows, the source stays a perfectly collimated beam along
    // /gps/direction instead of inheriting the state of the previous source.
    file_to_write << "/gps/ang/type planar" << std::endl;

    file_to_write << "/gps/ene/type Arb" << std::endl;
    file_to_write << "/gps/ene/diffspec 1" << std::endl;
    file_to_write << "/gps/hist/type arb" << std::endl;

    const double normalization =
        std::accumulate(this->fbin_count.begin(), this->fbin_count.end(), 0.);
    if (normalization <= 0.) {
      throw std::runtime_error(
          "The Ekin histogram is empty: the normalization is undefined "
          "(G4_hist::write_g4_ps_file)");
    }

    // The Ekin histogram is in MeV/u (per nucleon); GPS expects the total
    // kinetic energy per particle in GeV:  E[GeV] = E[MeV/u] * A / 1000
    for (size_t i = 0; i < this->fbin_count.size(); i++) {
      const double center =
          0.5 * (this->fright_bin_edge[i] + this->fright_bin_edge[i + 1]);
      file_to_write << "/gps/hist/point " << (center * imass) << " "
                    << this->fbin_count[i] / std::abs(normalization)
                    << " # ek = " << center << " MeV/u" << std::endl;
    }
    file_to_write << "/gps/hist/inter Lin" << std::endl;

    file_to_write << "# EOF " << std::endl;

    file_to_write.close();
    return 0;

  } else if (this->fhist_type == "X") {

    file_to_write << "/gps/hist/type biasx" << std::endl;

    const double norm = std::abs(
        *std::max_element(this->fbin_count.begin(), this->fbin_count.end()));
    if (norm <= 0.) {
      throw std::runtime_error(
          "The X histogram is empty: the normalization is undefined "
          "(G4_hist::write_g4_ps_file)");
    }

    const double offset = this->dimHalfX;

    const double max_element = 2 * this->dimHalfX;
    // GPS maps a sampled value v to x = v*2*halfx - halfx, so a data
    // position x has to be written as v = (x + halfx) / (2*halfx)
    for (size_t i = 0; i < this->fbin_count.size(); i++) {
      const double center =
          0.5 * (this->fright_bin_edge[i] + this->fright_bin_edge[i + 1]);
      file_to_write << "/gps/hist/point " << (center + offset) / max_element
                    << " " << this->fbin_count[i] / norm << std::endl;
    }

    file_to_write << "# EOF " << std::endl;

    file_to_write.close();
    return 0;

  } else if (this->fhist_type == "Y") {

    file_to_write << "/gps/hist/type biasy" << std::endl;

    const double norm = std::abs(
        *std::max_element(this->fbin_count.begin(), this->fbin_count.end()));
    if (norm <= 0.) {
      throw std::runtime_error(
          "The Y histogram is empty: the normalization is undefined "
          "(G4_hist::write_g4_ps_file)");
    }

    const double offset = this->dimHalfY;

    const double max_element = 2 * this->dimHalfY;
    for (size_t i = 0; i < this->fbin_count.size(); i++) {
      const double center =
          0.5 * (this->fright_bin_edge[i] + this->fright_bin_edge[i + 1]);
      file_to_write << "/gps/hist/point " << (center + offset) / max_element
                    << " " << this->fbin_count[i] / norm << std::endl;
    }

    file_to_write << "# EOF " << std::endl;

    file_to_write.close();
    return 0;

  } else if (this->fhist_type == "Theta") {

    file_to_write << "/gps/ang/type user" << std::endl;
    file_to_write << "/gps/hist/type theta" << std::endl;

    const double norm =
        std::accumulate(this->fbin_count.begin(), this->fbin_count.end(), 0.);
    if (norm <= 0.) {
      throw std::runtime_error(
          "The Theta histogram is empty: the normalization is undefined "
          "(G4_hist::write_g4_ps_file)");
    }
#ifdef DEBUG2
    std::cout << "Theta normalization = " << norm << std::endl;
#endif
    // The theta values are in radians, which is the unit Geant4 uses for the
    // user angular distributions (they are fed directly into sin/cos).
    for (size_t i = 0; i < this->fbin_count.size(); i++) {
      const double center =
          0.5 * (this->fright_bin_edge[i] + this->fright_bin_edge[i + 1]);
      double actual_value = this->fbin_count[i] / norm;
      if (actual_value < 1e-10) {
        actual_value = 0;
      }

      file_to_write << "/gps/hist/point " << center << " " << actual_value
                    << std::endl;
    }

    file_to_write << "# EOF " << std::endl;

    file_to_write.close();
    return 0;

  } else if (this->fhist_type == "Phi") {

    // /gps/ang/type user is set by the Theta macro; the converter only
    // generates the Phi macro when the Theta macro exists as well (see
    // CommandGenerator::generateCommandFile).
    file_to_write << "/gps/hist/type phi" << std::endl;
    const double norm =
        std::accumulate(this->fbin_count.begin(), this->fbin_count.end(), 0.);
    if (norm <= 0.) {
      throw std::runtime_error(
          "The Phi histogram is empty: the normalization is undefined "
          "(G4_hist::write_g4_ps_file)");
    }
    // The phi values are in radians (see the Theta branch).
    for (size_t i = 0; i < this->fbin_count.size(); i++) {
      const double center =
          0.5 * (this->fright_bin_edge[i] + this->fright_bin_edge[i + 1]);
      double actual_value = this->fbin_count[i] / norm;
      if (actual_value < 1e-10) {
        actual_value = 0;
      }
      file_to_write << "/gps/hist/point " << center << " " << actual_value
                    << std::endl;
    }

    file_to_write << "# EOF " << std::endl;

    file_to_write.close();
    return 0;
  }

  return 0;
}

/**
 * @brief Write the particle definition commands of the GPS source.
 *
 * Emits /gps/particle <name> for special particles (Z = -1), neutron
 * (Z = 0) or proton (Z = 1); for Z >= 2 it emits /gps/particle ion plus
 * /gps/ion Z A, where A is the most abundant stable isotope
 * (G4_hist::associated_A) and the nuclear mass in u is written as a
 * comment.
 *
 * @param to_write open file stream to write on
 * @param particle particle name as in the ROOT file
 * @param Z charge Z, -1 for special particles
 * @throw std::runtime_error if no isotope is defined for the charge Z
 * */
void G4_hist::define_source(std::ofstream &to_write, std::string particle,
                            int Z) {
  if (Z == -1) {
    // exotic particles
    to_write << "/gps/particle " << particle << std::endl;
  } else if (Z == 0) {

    to_write << "/gps/particle neutron" << std::endl;
  } else if (Z == 1) {
    const auto a_it = this->associated_A.find(Z);
    if (a_it == this->associated_A.end()) {
      throw std::runtime_error("No isotope defined for Z = " +
                               std::to_string(Z) + " (G4_hist::define_source)");
    }
    const int A = a_it->second;
    double imass = G4NucleiProperties::GetNuclearMass(A, Z) / CLHEP::amu_c2;

    to_write << "/gps/particle proton # with u = " << imass << std::endl;
  } else {
    const auto a_it = this->associated_A.find(Z);
    if (a_it == this->associated_A.end()) {
      throw std::runtime_error("No isotope defined for Z = " +
                               std::to_string(Z) + " (G4_hist::define_source)");
    }
    const int A = a_it->second;
    double imass = G4NucleiProperties::GetNuclearMass(A, Z) / CLHEP::amu_c2;
    to_write << "/gps/particle ion" << std::endl;
    to_write << "/gps/ion " << Z << " " << A << " # with u = " << imass
             << std::endl;
  }
}

/**
 * @brief Generate the sim.mac file sourcing every generated GPS macro; see
 * the header for the full documentation.
 *
 * For each particle, the Ek macro is required (it creates the source: a
 * missing Ek macro skips the whole particle), X/Y are optional (a missing
 * position falls back to the GPS default uniform distribution) and Theta +
 * Phi must be present together (otherwise the source stays collimated).
 * */
void CommandGenerator::generateCommandFile(
    const std::string &commandFileName, const double scoring_dimensions_mm,
    const std::string &configuration, const std::vector<std::string> particles
    // int maxZ,
) {
  // Open the file for writing
  std::ofstream commandFile(commandFileName);

  // Check if the file is successfully opened
  if (!commandFile.is_open()) {
    std::cerr << "Error opening the command file." << std::endl;
    return;
  }
  commandFile << "###########################################" << std::endl;
  commandFile << "# Phase space sourcing for Geant4 GPS." << std::endl;
  commandFile << "# NOTE: run the Geant4 application from the directory "
                 "containing this file and the hist_geant_*.mac files, "
                 "using this file as macro."
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
    commandFile << "# File generated with configuration = " << configuration
                << std::endl;
  }
  commandFile << "# scoring dimensions = " << scoring_dimensions_mm << " mm"
              << std::endl;
  commandFile << "# This file is automatically generated." << std::endl;
  commandFile << "###########################################" << std::endl;

  // commandFile << "/run/initialize" << std::endl;

  // The macro files are generated in the same folder as sim.mac; the
  // existence checks below (and the /control/execute lines) therefore use
  // names relative to that folder.
  const std::filesystem::path outDir =
      std::filesystem::path(commandFileName).parent_path();
  auto fileExists = [&outDir](const std::string &fileName) {
    return std::filesystem::exists(outDir / fileName);
  };

  // Write the commands to the file
  // for (int Z = 0; Z <= maxZ; ++Z) {
  for (auto particl : particles) {
    // The Ek macro is required: it creates the source. Without it the
    // X/Y/Theta/Phi macros of this particle would silently configure the
    // PREVIOUS particle's source, so the whole particle is skipped.
    std::ostringstream ekPath;
    ekPath << "hist_geant_Ek_" << configuration << "_particle_" << particl
           << ".mac";
    if (!fileExists(ekPath.str())) {
      std::cerr << "Warning: " << ekPath.str() << " is missing: particle "
                << particl
                << " is skipped entirely (its other macros would corrupt the "
                   "previous source)."
                << std::endl;
      continue;
    }
    commandFile << "/control/execute " << ekPath.str() << std::endl;

    // Iterate through the position characteristics: X, Y
    for (const std::string &feature : {"X", "Y"}) {
      std::ostringstream filePath;
      filePath << "hist_geant_" << feature << "_" << configuration
               << "_particle_" << particl << ".mac";
      if (fileExists(filePath.str())) {
        commandFile << "/control/execute " << filePath.str() << std::endl;
      } else {
        std::cerr << "Warning: " << filePath.str() << " is missing: the "
                  << feature << " position of particle " << particl
                  << " falls back to the GPS default (uniform over the "
                     "plane)."
                  << std::endl;
      }
    }

    // Theta and Phi must be used together (they form the GPS user angular
    // distribution). If only one of them exists, both are skipped and the
    // source stays perfectly collimated (see the /gps/ang/type planar reset
    // in the Ek macro).
    std::ostringstream thetaPath, phiPath;
    thetaPath << "hist_geant_Theta_" << configuration << "_particle_" << particl
              << ".mac";
    phiPath << "hist_geant_Phi_" << configuration << "_particle_" << particl
            << ".mac";
    const bool hasTheta = fileExists(thetaPath.str());
    const bool hasPhi = fileExists(phiPath.str());
    if (hasTheta && hasPhi) {
      commandFile << "/control/execute " << thetaPath.str() << std::endl;
      commandFile << "/control/execute " << phiPath.str() << std::endl;
    } else if (hasTheta || hasPhi) {
      std::cerr << "Warning: Theta/Phi macros incomplete for particle "
                << particl << " (Theta: " << (hasTheta ? "yes" : "no")
                << ", Phi: " << (hasPhi ? "yes" : "no")
                << "): both are skipped, the source stays perfectly "
                   "collimated."
                << std::endl;
    }
  }

  commandFile << "###########################################" << std::endl;
  commandFile << "# EOF" << std::endl;

  // Close the file
  commandFile.close();

  std::cout << "Command file generated: " << commandFileName << std::endl;
}
