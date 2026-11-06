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
 * @file G4_PS_hist.cpp
 *
 * @brief Class implementation containing histogram information
 *
 * @author enrico.pierobon
 * Contact: e.pierobon@gsi.de
 *
 */

#include "G4_PS_hist.hpp"
#include <numeric>
#include <stdexcept>
#include <utility>

/**
 * @brief Constructor; see the header for the full documentation.
 *
 * Validates the input, computes the abundance histogram key
 * (fabundance_root_hist_string) and the index of the selected scoring
 * dimensions (index_of_dimensions), and creates the output directory if it
 * does not exist yet.
 * */
G4_PS_hist::G4_PS_hist(std::filesystem::path path_to_root_file,
                       std::filesystem::path outputFolder,
                       const std::string configuration,
                       const std::string dimensions)
    : fpath_to_root_file(path_to_root_file), foutputFolder(outputFolder),
      fconfiguration(configuration), fdimensions(dimensions) {

  if (this->REQUIRED_POINTS_PER_HIST <= 1) {
    std::cerr << "Warning: required points per hist cannot be <= 1, got "
              << this->REQUIRED_POINTS_PER_HIST << G4endl;
    throw std::invalid_argument("\"REQUIRED_POINTS_PER_HIST\" variable is not "
                                "valid, please adjust it in the code");
  }
#ifdef DEBUG
  std::cout << "Creating class with parameters: " << std::endl;
  std::cout << "\t path to root file = " << path_to_root_file << std::endl;
#endif
  std::string abundance_root_hist_string = "";
  if (this->fdimensions == "") {
    abundance_root_hist_string = "Ekin_" + configuration;
    this->index_of_dimensions = 0;

  } else {
    abundance_root_hist_string = "Ekin_" + dimensions + "_" + configuration; /*
                                            string from where the relative
                             abundance is calculated this histogram will be used
                             to estimate the retlative abundance for each Z
                                            */
    if (this->fdimensions == this->AVAILABLE_DIMENSIONS[1]) {
      index_of_dimensions = 1;
    } else if (this->fdimensions == this->AVAILABLE_DIMENSIONS[2]) {
      index_of_dimensions = 2;
    } else {
      // error
      std::cerr << "Error: cannot associate the dimensions of the scoring to a "
                   "valid root key"
                << std::endl;
      throw std::runtime_error(
          "index of dimensions not set G4_PS_hist::G4_PS_hist()");
    }
  }
  std::cout << "root key to calculate relative abundance = "
            << abundance_root_hist_string << std::endl;
  this->fabundance_root_hist_string = abundance_root_hist_string;

  // make sure the output directory exists
  std::error_code ec;
  std::filesystem::create_directories(this->foutputFolder, ec);
  if (ec) {
    throw std::runtime_error("Cannot create the output directory \"" +
                             this->foutputFolder.string() +
                             "\": " + ec.message());
  }
}

/**
 * @brief Compute the relative abundances, then generate all the output
 * files.
 *
 * @return the list of processed particle identifiers; a single "-1" if no
 * Ekin histogram matched the selection
 * */
std::vector<std::string> G4_PS_hist::process_file() {
  this->calculate_relative_abundance(this->fabundance_root_hist_string);
  std::vector<std::string> max_z = this->read_root_file();
  return max_z;
}

/**
 * @brief Compute the relative abundances from the ROOT histograms whose key
 * contains tree_to_use.
 *
 * Each particle abundance is its histogram total content divided by the sum
 * of the contents of all the matched histograms, so the result in
 * frelative_abundance sums to 1.
 *
 * @throw std::runtime_error if the file cannot be opened or the total
 * content is zero
 * */
void G4_PS_hist::calculate_relative_abundance(const std::string tree_to_use) {
  std::vector<std::string> particles;
  std::vector<double> counts;
  std::unique_ptr<TFile> myFile(
      TFile::Open(this->fpath_to_root_file.string().c_str()));

  if (!myFile || myFile->IsZombie()) {
    std::cerr << "Error: file " << this->fpath_to_root_file.string()
              << " not existing or corrupted" << std::endl;
    throw std::runtime_error("File not found or corrupted");
  }
  TIter next(myFile->GetListOfKeys());
  TKey *key;
  std::tuple<int, std::string, std::string> tree_parser_obj;

  while ((key = static_cast<TKey *>(next()))) {

    std::string tree_name = key->GetName();

    if (tree_name.find(tree_to_use) != std::string::npos) {
/* if the filename matches the tree to use file: */
#ifdef DEBUG_REL_ABUNDANCE
      std::cout << "using: " << tree_name << " to calculate relative abundance"
                << std::endl;
      std::cout << "Found match with " << tree_to_use << std::endl;
#endif // DEBUG2

#ifdef DEBUG2
      std::cout << "Parsing tree object: " << tree_name << std::endl;
#endif
      tree_parser_obj = this->parse_tree_key(tree_name);
#ifdef DEBUG2
      std::cout << "Parsing complete of: " << tree_name << std::endl;
#endif
      particles.push_back(
          std::get<1>(tree_parser_obj)); // push bach the information on Z

      /*
       * add the information on the counts
       * */
      counts.push_back(
          this->get_hist_total_counts((TH1D *)myFile->Get(key->GetName())));
    }

  } // iteration over the tree keys
  /*calculate the total number of elements */
  double total_counts = 0;
  /* previus definition using the total counts */
  for (auto a : counts) {
    total_counts += a;
  }

  if (total_counts <= 0.) {
    std::cerr << "Error: the total number of counts of the abundance "
                 "histograms is zero (\""
              << tree_to_use << "\"); the relative abundances are undefined."
              << std::endl;
    throw std::runtime_error(
        "All abundance histograms are empty: cannot calculate the relative "
        "abundances (G4_PS_hist::calculate_relative_abundance)");
  }

  /* new definition using the max. */
  /* NOTE: this it not working correctly on reshooting. */
  // total_counts = *std::max_element(counts.begin(), counts.end());
  /* calculate the relative abundance */
  std::vector<double> relative_counts_vector;

  for (auto a : counts) {
    relative_counts_vector.push_back(a / total_counts);
  }

  this->frelative_abundance = std::make_pair(particles, relative_counts_vector);
  for (size_t i = 0; i < this->frelative_abundance.first.size(); i++) {
    std::cout << "particle = " << this->frelative_abundance.first[i]
              << " counts = " << this->frelative_abundance.second[i]
              << std::endl;
  }

  std::cout << "End of calculating relative abundances." << std::endl;
}

/**
 * @brief Read the ROOT file and generate the CSV and GPS macro files.
 *
 * Iterates over the tree keys, keeps the histograms matching the requested
 * dimensions/configuration and not belonging to "Null", re-bins each of
 * them and writes one CSV file (trimmed data) and one Geant4 GPS macro per
 * kept histogram. See the header for the full behavior (skipping rules,
 * intensity assignment).
 *
 * @return the list of particle identifiers that have an Ekin histogram; a
 * single "-1" if none matched
 * */
std::vector<std::string> G4_PS_hist::read_root_file() {
  std::vector<std::string>
      processed_particle; // vector of the Z that have a histogram and
                          // have been processed

  // Open the file
  std::unique_ptr<TFile> myFile(
      TFile::Open(this->fpath_to_root_file.string().c_str()));

  if (!myFile || myFile->IsZombie()) {
    std::cerr << "Error: file " << this->fpath_to_root_file.string()
              << " not existing or corrupted" << std::endl;
    return std::vector<std::string>{"-1"};
  }

  // iterate over the keys
  TIter next(myFile->GetListOfKeys());
  TKey *key;
  std::tuple<int, std::string, std::string> config_Z;
  // std::vector<double> relative_abundance; // relative abundace array

  // Geant4 GPS always starts with a default source 0 (a 1 MeV geantino with
  // intensity 1.0). The first particle must therefore configure that source
  // with /gps/source/intensity; adding N sources on top of it with
  // /gps/source/add would make the geantino 50 % of all primaries, since the
  // relative abundances sum to 1.
  bool first_source_written = false;

  // read the tuple values
  while ((key = static_cast<TKey *>(next()))) {
    config_Z = this->parse_tree_key(std::string(key->GetName()));
    int setup = std::get<0>(config_Z);
    std::string particle = std::get<1>(config_Z);
    std::string name = std::get<2>(config_Z);
    std::string file_name = key->GetName();
    bool isAGoodKey = false;
    if (this->fdimensions == "") {
      // check that it does not contains wide or very_wude
      if (file_name.find(G4_PS_hist::AVAILABLE_DIMENSIONS[1]) ==
              std::string::npos &&
          file_name.find(G4_PS_hist::AVAILABLE_DIMENSIONS[2]) ==
              std::string::npos) {
        isAGoodKey = true;
      } else {
        isAGoodKey = false;
      }
    } // empty string
    else { // strig is not empy
      if (this->fdimensions == G4_PS_hist::AVAILABLE_DIMENSIONS[1]) {
        // check that the second substring is not contained
        if ((file_name.find(G4_PS_hist::AVAILABLE_DIMENSIONS[2]) ==
             std::string::npos) &&
            (file_name.find(G4_PS_hist::AVAILABLE_DIMENSIONS[1]) !=
             std::string::npos)) {
          isAGoodKey = true;
        } else {
          isAGoodKey = false;
        }
      } // user selected first element of available dimensions
      else if (this->fdimensions == G4_PS_hist::AVAILABLE_DIMENSIONS[2]) {
        // check if the substring is there
        if (file_name.find(G4_PS_hist::AVAILABLE_DIMENSIONS[2]) !=
            std::string::npos) {
          isAGoodKey = true;
        } else {
          isAGoodKey = false;
        }
      } // user selected second elment of the available dimensions
      else {
        // fail, the program do not hanlde this substring selection
        std::cerr << "Error: the program cannot decodode the selected "
                     "dimensions. Available selections are: [";
        for (auto a : G4_PS_hist::AVAILABLE_DIMENSIONS) {
          std::cerr << a << ", ";
        }
        std::cerr << "]" << std::endl;
        std::cerr << "Key under evaluation = " << file_name
                  << ". Could not find " << this->fdimensions << "."
                  << std::endl;
        throw std::invalid_argument("Invalid dimension selection for tree key "
                                    "decoding. G4_PS_hist::read_root_file. ");
      }
    }
    if (isAGoodKey && particle != "Null" &&
        setup ==
            std::stoi(
                this->fconfiguration)) // full string can be found in the file
    {
      std::cout << "Processing root tree: " << key->GetName() << std::endl;
      // opening the histogram
      // guard against non-TH1 objects (an invalid cast would be undefined
      // behavior)
      TObject *obj = myFile->Get(key->GetName());
      if (!obj || !obj->InheritsFrom(TH1::Class())) {
        std::cerr << "Warning: object '" << key->GetName()
                  << "' is not a TH1 "
                     "(got "
                  << (obj ? obj->ClassName() : "null")
                  << "); it is being skipped." << std::endl;
        continue;
      }
      TH1D *hist1 = (TH1D *)obj;
      // getting the dimension of the bins

      // std::pair<std::array<double, n_bins>, std::array<double, n_bins>>
      // root_hist_values;
      std::pair<std::vector<double>, std::vector<double>> root_hist_values;
      std::pair<std::vector<double>, std::vector<double>>
          root_hist_values_trimmed;

      // extract the root values
      root_hist_values = this->convert_root_hist_to_array(hist1);

#ifdef DEBUG2
      // DEBUG
      std::cout << "Root hist = " << std::endl;
      for (size_t a = 0; a < root_hist_values.first.size(); a++) {
        std::cout << root_hist_values.first[a] << ", "
                  << root_hist_values.second[a] << std::endl;
      }
      std::cout << "END" << std::endl;
#endif

      if (std::reduce(root_hist_values.second.begin(),
                      root_hist_values.second.end()) <=
          this->REQUIRED_POINTS_PER_HIST) {
        std::cout << "\t Warning: no valuable data points found for this "
                     "histogram: histogram has less ("
                  << std::reduce(root_hist_values.second.begin(),
                                 root_hist_values.second.end())
                  << ") then required counts ("
                  << this->REQUIRED_POINTS_PER_HIST << ")" << std::endl;
        continue;
      }
      root_hist_values_trimmed = this->trim_zeros_from_hist(root_hist_values);
#ifdef DEBUG2
      // DEBUG
      std::cout << "Root hist trimmed = " << std::endl;
      for (size_t a = 0; a < root_hist_values_trimmed.first.size(); a++) {
        std::cout << root_hist_values_trimmed.first[a] << ", "
                  << root_hist_values_trimmed.second[a] << std::endl;
      }
#endif

      if (root_hist_values_trimmed.first.size() < 2) {
        // A single (center, density) point does not define a distribution:
        // Geant4's arbitrary energy distribution is also broken by a single
        // point (division by zero in the normalization of the cumulative
        // area).
        std::cout << "\t Warning: fewer than 2 points with data left after "
                     "trimming; no distribution can be built for this "
                     "histogram, it is being skipped."
                  << std::endl;
        continue;
      }

      // create a new binning: at most N_POINTS_MAX edges, but at least 3
      // (the GPS arbitrary distribution needs at least 2 points, i.e. 3
      // edges)
      size_t N_points = std::min(root_hist_values_trimmed.first.size(),
                                 static_cast<size_t>(G4_PS_hist::N_POINTS_MAX));
      if (N_points > root_hist_values_trimmed.first.size()) {
        std::cout << "\tInfo: root points ("
                  << root_hist_values_trimmed.first.size()
                  << ") are less than the number of required points, adjusting "
                     "to not oversample."
                  << std::endl;
      }
      if (N_points < 3) {
        N_points = 3;
      }
      const double min_element =
          *std::min_element(std::begin(root_hist_values_trimmed.first),
                            std::end(root_hist_values_trimmed.first));
      const double max_element =
          *std::max_element(std::begin(root_hist_values_trimmed.first),
                            std::end(root_hist_values_trimmed.first));

      std::vector<double> new_binning =
          this->linspace(min_element, max_element, N_points);
#ifdef DEBUG2
      // DEBUG:
      std::cout << "Got a linspace:" << std::endl;
      for (size_t a = 0; a < new_binning.size(); a++) {
        std::cout << a << ", " << new_binning[a] << std::endl;
      }

#endif

      std::vector<double> rebin_y_count =
          this->rebin_histogram(root_hist_values_trimmed, new_binning);

#ifdef DEBUG2
      // DEBUG
      std::cout << "Rebinned histogram:" << std::endl;
      for (size_t a = 0; a < rebin_y_count.size(); a++) {
        std::cout << 0.5 * (new_binning[a] + new_binning[a + 1]) << ", "
                  << rebin_y_count[a] << std::endl;
      }
#endif

      std::filesystem::path file_name_root =
          foutputFolder / (std::string("hist_root_") + name + std::string("_") +
                           std::to_string(setup) + std::string("_Z_") +
                           particle + std::string(".csv"));

      this->write_hist_to_csv_file(file_name_root, root_hist_values_trimmed);

      std::filesystem::path file_name_geant =
          foutputFolder /
          (std::string("hist_geant_") + name + std::string("_") +
           std::to_string(setup) + std::string("_particle_") + particle +
           std::string(".mac"));

      double relative_abundance_at_given_z = 0;

      for (size_t a = 0; a < frelative_abundance.first.size(); a++) {
        if (frelative_abundance.first[a] == particle) {
          relative_abundance_at_given_z = frelative_abundance.second[a];
          if (name == "Ek") {
#ifdef DEBUG
            std::cout << "Pushing to processed partilce =  " << particle
                      << std::endl;
#endif
            processed_particle.push_back(particle);
          }
        }
      }

      bool is_first_source = (name == "Ek") && !first_source_written;

      G4_hist G4_h(file_name_geant, name,
                   "hist key: " + std::string(key->GetName()),
                   "Z = " + particle, "", particle, new_binning, rebin_y_count,
                   relative_abundance_at_given_z, is_first_source);
      if (is_first_source) {
        first_source_written = true;
      }

      // Set the half dimensions:
      G4_h.dimHalfX = this->SOURCE_DIMENSIONS[this->index_of_dimensions];
      G4_h.dimHalfY = this->SOURCE_DIMENSIONS[this->index_of_dimensions];
#ifdef DEBUG2
      std::cout << "Data in the class:" << std::endl;
      G4_h.print_data();
      std::cout << "END" << std::endl;
#endif

#ifdef DEBUG
      G4_h.save_data_as_csv();
#endif
      G4_h.write_g4_ps_file();
    } else {
    }
  }
  if (processed_particle.empty()) {
    std::cerr << "Error: no Ekin histogram matched configuration = \""
              << this->fconfiguration << "\" and dimensions = \""
              << this->fdimensions << "\". Nothing can be generated."
              << std::endl;
    std::cerr << "Available configurations:";
    for (size_t i = 0; i < G4_PS_hist::AVAILABLE_CONFIGURATION.size(); i++) {
      std::cerr << " " << (i + 1) << "="
                << G4_PS_hist::AVAILABLE_CONFIGURATION[i];
    }
    std::cerr << std::endl;
    return std::vector<std::string>{"-1"};
  }
  return processed_particle;
}

/**
 * @brief Write the (center, density) histogram points to a CSV file.
 *
 * The file is overwritten if it already exists; the two columns are
 * whitespace separated.
 *
 * @return 0 on success
 * */
int G4_PS_hist::write_hist_to_csv_file(
    std::filesystem::path p,
    std::pair<std::vector<double>, std::vector<double>> hist_to_write) {
  std::ofstream file_to_write(p, std::ofstream::out);

  for (size_t i = 0; i < hist_to_write.first.size(); i++) {
    file_to_write << hist_to_write.first[i] << ", " << hist_to_write.second[i]
                  << std::endl;
  }

  return 0;
}

/**
 * @brief Rebin a histogram onto a new grid of bin edges.
 *
 * The old points are treated as a piecewise-linear pdf; the exact integral
 * over each new bin (trapezoids per old segment, 0 outside the old domain)
 * divided by the new bin width conserves the total integral and correctly
 * handles split bins. The scan is O(n_old + n_new).
 *
 * @return one mean density per new bin (size = bin_edges.size() - 1)
 * */
std::vector<double> G4_PS_hist::rebin_histogram(
    std::pair<std::vector<double>, std::vector<double>> h1_data,
    std::vector<double> bin_edges) {

  const auto &x_old = h1_data.first;
  const auto &y_old = h1_data.second;

  const size_t n_new_bins = bin_edges.size() >= 2 ? bin_edges.size() - 1 : 0;

  std::vector<double> result(n_new_bins, 0.0);

  if (x_old.empty() || x_old.size() != y_old.size()) {
    return result;
  }

  // Assumes uniformly spaced old bin centers.
  const double old_width = x_old.size() > 1 ? x_old[1] - x_old[0] : 1.0;

  if (!(old_width > 0.0)) {
    return result;
  }

  size_t j = 0;

  for (size_t k = 0; k < n_new_bins; ++k) {
    const double a = bin_edges[k];
    const double b = bin_edges[k + 1];

    if (!(b > a)) {
      continue;
    }

    while (j < x_old.size() && x_old[j] + 0.5 * old_width <= a) {
      ++j;
    }

    size_t i = j;
    double count = 0.0;

    while (i < x_old.size() && x_old[i] - 0.5 * old_width < b) {
      const double old_a = x_old[i] - 0.5 * old_width;
      const double old_b = x_old[i] + 0.5 * old_width;

      const double overlap =
          std::max(0.0, std::min(b, old_b) - std::max(a, old_a));

      count += y_old[i] * overlap / old_width;
      ++i;
    }

    result[k] = count / (b - a);
  }

  return result;
}

std::
    tuple<int, std::string, std::string>
    /**
     * @brief Parse a ROOT tree key into (configuration, particle, type).
     *
     * The type is the first recognized substring (Ekin/X/Y/Theta/Phi); the key
     * is then matched against `_<configuration>_(<particle>)`, where the
     * particle is a charge Z or one of G4_PS_hist::fspecialParticles. Unmatched
     * fields are "Null" / -1.
     *
     * Example: `Ekin_wide_1500_2` -> (1500, "2", "Ek"); `Ekin_1_e+` -> (1,
     * "e+", "Ek").
     * */
    G4_PS_hist::parse_tree_key(std::string s) {
// DEBUG: print the string
#ifdef DEBUG_REGEX
  std::cout << "string to analyse: " << s << std::endl;
#endif
  std::string type = "None";

  if (s.find(this->ek_string) != std::string::npos) {
    type = "Ek";
  } else if (s.find(this->x_pos_string) != std::string::npos) {
    type = "X";
  } else if (s.find(this->y_pos_string) != std::string::npos) {
    type = "Y";
  } else if (s.find(this->theta_ang_string) != std::string::npos) {
    type = "Theta";
  } else if (s.find(this->phi_ang_string) != std::string::npos) {
    type = "Phi";
  } else {
    type = "None";
  }
#ifdef DEBUG_REGEX
  std::cout << "String type = " << type << std::endl;
#endif // DEBUG_REGEX

  if (type != "None") { // look for Ekin string

    /*
     * use a regex to match the first and the second digit separated by an
     * '_' character: Ekin_1_3 will match 1 and 3 Ekin_1_4 will mathc 1 and
     * 4 etc..
     * */
    // std::regex re("(\\d+)_(\\d+)_(\\d+)"); // old

    /*
     * Update the reject to match new format:
     * format is the type of:
     * quantity_type_distance_z
     * example:  'Ekin_very_wide_100_2'
     *
     * to update:
     * work with:
     * Ekin_very_wide_100_e-
     * Ekin_very_wide_100_e+
     *
     * */

    // std::regex re("_(\\d+)_(\\d+)"); // old

    // build the particle pttern
    std::string particlePattern;

    for (const auto &p : fspecialParticles) {
      if (!particlePattern.empty())
        particlePattern += "|";

      particlePattern += p;
    }

    std::regex re("_(\\d+)_(" + particlePattern + "|\\d+)");
    std::smatch m;

    std::regex_search(s, m, re);
#ifdef DEBUG_REGEX
    std::cout << "Got a match; full = " << m[-1] << std::endl;
    std::cout << "\t g1 = config  = " << m[1] << std::endl;
    std::cout << "\t g2 = Z = " << m[2] << std::endl;
#endif // DEBUG_REGEX

    int config = -1; // -1: not found / not parseable
    std::string particle;
    try {
      config = std::stoi(m[1]);
      particle = m[2];
    } catch (std::exception &e) {
      std::cerr << "ERROR: cannot convert tree name to extract config and Z"
                << std::endl;
      std::cerr << "\t tree name: " << s << std::endl;
    }

    return std::tuple(config, particle, type);

  } else {
    // if it is not the relative abundance histogram
    if (s != "Z_1") {
      std::cerr << "Cannot associate tree name to known class parameters;"
                << std::endl;
      std::cerr << "\t tree name '" << s << "' is being skipped." << std::endl;
    }
    return std::tuple(-1, "Null", "Null");
  }
}

/**
 * @brief N_elements equally spaced points in [min, max] (endpoints
 * included).
 *
 * Example: linspace(1, 10, 4) -> {1, 4, 7, 10}.
 * */
std::vector<double> G4_PS_hist::linspace(const double min, const double max,
                                         const uint N_elements) {
  std::vector<double> to_return;

  if (N_elements < 2) {
    return to_return;
  }

  for (size_t i = 0; i < N_elements; i++) {
    to_return.push_back((max - min) * i / (N_elements - 1) + min);
  }
  return to_return;
}

std::
    pair<std::vector<double>, std::vector<double>>
    /**
     * @brief Convert a ROOT histogram into (bin center, density) arrays.
     *
     * The density is the bin content divided by the bin width, so the
     * representation also stays correct for non-uniform binning. Underflow and
     * overflow are excluded; degenerate (zero-width) bins are skipped.
     * */
    G4_PS_hist::convert_root_hist_to_array(TH1D *hist) {
  std::vector<double> bin_center;
  std::vector<double> bin_density;

  const int n_bins = hist->GetNbinsX();
  for (int i = 1; i <= n_bins; i++) {
    const double left = hist->GetBinLowEdge(i);
    const double right = hist->GetBinLowEdge(i + 1);
    const double width = right - left;
    if (width <= 0.) {
      continue; // degenerate bin
    }
    bin_center.push_back(0.5 * (left + right));
    // density (content per unit of width): makes the representation correct
    // also for non-uniform binning
    bin_density.push_back(hist->GetBinContent(i) / width);
  }
  return std::make_pair(bin_center, bin_density);
}

/**
 * @brief Remove the leading and trailing zero-density points from a
 * (center, density) histogram.
 *
 * If every point is zero an empty pair is returned.
 * */
std::pair<std::vector<double>, std::vector<double>>
G4_PS_hist::trim_zeros_from_hist(
    std::pair<std::vector<double>, std::vector<double>> h1_data) {
  const std::vector<double> &counts = h1_data.second;
  const size_t size = h1_data.first.size();
  std::vector<double> trimmed_center;
  std::vector<double> trimmed_counts;

  if (size == 0) {
    return std::make_pair(trimmed_center, trimmed_counts);
  }

  // first point with data
  size_t lowerFirst = 0;
  while (lowerFirst < size && counts[lowerFirst] == 0.) {
    lowerFirst++;
  }
  // last point with data
  size_t upperFirst = size - 1;
  while (upperFirst > lowerFirst && counts[upperFirst] == 0.) {
    upperFirst--;
  }

  if (lowerFirst > upperFirst) {
    // all points are zero
    std::cout << "No valid points found" << std::endl;
    return std::make_pair(trimmed_center, trimmed_counts);
  }

  for (size_t j = lowerFirst; j <= upperFirst; j++) {
    trimmed_center.push_back(h1_data.first[j]);
    trimmed_counts.push_back(counts[j]);
  }

  return std::make_pair(trimmed_center, trimmed_counts);
}

bool /**
      * @brief true if ref is inside [x_prec, x_succ).
      *
      * @note legacy helper, currently unused.
      * */
G4_PS_hist::isWithinInterval(const double x_prec, const double x_succ,
                             const int ref) {

  const double ref_d = ref;

  /*
   * if ref is included in the interval
   * */
  if (ref_d >= x_prec && ref_d < x_succ) {
    return true;
  } else {
    return false;
  }
}

bool /**
      * @brief true if |ref - x_prec| < |ref - x_succ|.
      *
      * @note legacy helper, currently unused.
      * */
G4_PS_hist::isUpperEdgeCloser(const double x_prec, const double x_succ,
                              const int ref) {
  const double ref_d = ref;

  const double delta_prec = std::abs(ref_d - x_prec);
  const double delta_succ = std::abs(ref_d - x_succ);

  if (delta_prec < delta_succ) {
    return true;
  } else {
    return false;
  }
}

/**
 * @brief calcualte of of all the contents inside the hisgraam, including under
 * and overflow
 * */
double G4_PS_hist::get_hist_total_counts(TH1D *hist) {

  double hist_counts = 0;
  for (size_t i = 0; i <= hist->GetNbinsX() + 1; i++) {
    hist_counts += hist->GetBinContent(i);
  }
  hist_counts += hist->GetBinContent(0);                     // undeflow
  hist_counts += hist->GetBinContent(hist->GetNbinsX() + 1); // overlof
  return hist_counts;
}

/**
 * @brief calculate the itegral of the distribution
 * @return integral of the histogram (bin-dependend)
 * */

double G4_PS_hist::get_hist_integral(TH1D *hist) {
  return hist->Integral("width");
}
