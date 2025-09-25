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

G4_PS_hist::G4_PS_hist(std::filesystem::path path_to_root_file,
                       std::filesystem::path outputFolder,
                       const std::string configuration,
                       const std::string dimensions)
    : path_to_root_file(path_to_root_file), outputFolder(outputFolder),
      configuration(configuration), dimensions(dimensions) {
#ifdef DEBUG
  std::cout << "Creating class with parameters: " << std::endl;
  std::cout << "\t path to root file = " << path_to_root_file << std::endl;
#endif
  std::string abundance_root_hist_string = "";
  if (this->dimensions == "") {
    abundance_root_hist_string = "Ekin_" + configuration;
    this->index_of_dimensions = 0;

  } else {
    abundance_root_hist_string = "Ekin_" + dimensions + "_" + configuration; /*
                                            string from where the relative
                             abundance is calculated this histogram will be used
                             to estimate the retlative abundance for each Z
                                            */
    if (this->dimensions == this->AVAILABLE_DIMENSIONS[1]) {
      index_of_dimensions = 1;
    } else if (this->dimensions == this->AVAILABLE_DIMENSIONS[2]) {
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
  this->abundance_root_hist_string = abundance_root_hist_string;
}

std::vector<int> G4_PS_hist::process_file() {
  this->calculate_relative_abundance(this->abundance_root_hist_string);
  std::vector<int> max_z = this->read_root_file();
  return max_z;
}

void G4_PS_hist::calculate_relative_abundance(const std::string tree_to_use) {
  std::vector<unsigned int> Z;
  std::vector<double> counts;
  std::unique_ptr<TFile> myFile(
      TFile::Open(this->path_to_root_file.string().c_str()));

  if (!myFile || myFile->IsZombie()) {
    std::cerr << "Error: file " << this->path_to_root_file.string()
              << " not existing or corrupted" << std::endl;
    throw std::runtime_error("File not found or corrupted");
  }
  TIter next(myFile->GetListOfKeys());
  TKey *key;
  std::tuple<int, int, std::string> tree_parser_obj;

  while ((key = (TKey *)next())) {

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
      Z.push_back(
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

  /* new definition using the max. */
  /* NOTE: this it not working correctly on reshooting. */
  // total_counts = *std::max_element(counts.begin(), counts.end());
  /* calculate the relative abundance */
  std::vector<double> relative_counts_vector;

  for (auto a : counts) {
    relative_counts_vector.push_back(a / total_counts);
  }

  this->relative_abundance = std::make_pair(Z, relative_counts_vector);
  for (size_t i = 0; i < this->relative_abundance.first.size(); i++) {
    std::cout << "Z = " << this->relative_abundance.first[i]
              << " counts = " << this->relative_abundance.second[i]
              << std::endl;
  }

  std::cout << "End of calculating relative abundances." << std::endl;
}

std::vector<int> G4_PS_hist::read_root_file() {
  std::vector<int> processed_z; // vector of the Z that have a histogram and
                                // have been processed

  // Open the file
  std::unique_ptr<TFile> myFile(
      TFile::Open(this->path_to_root_file.string().c_str()));

  if (!myFile || myFile->IsZombie()) {
    std::cerr << "Error: file " << this->path_to_root_file.string()
              << " not existing or corrupted" << std::endl;
    return std::vector<int>{-1};
  }

  // iterate over the keys
  TIter next(myFile->GetListOfKeys());
  TKey *key;
  TCanvas c1;
  std::tuple<int, int, std::string> config_Z;
  // std::vector<double> relative_abundance; // relative abundace array

  // read the tuple values
  while ((key = (TKey *)next())) {
    config_Z = this->parse_tree_key(std::string(key->GetName()));
    int setup = std::get<0>(config_Z);
    int Z = std::get<1>(config_Z);
    std::string name = std::get<2>(config_Z);
    std::string file_name = key->GetName();
#ifdef DEBUG3
    if (Z == 21 && std::string(key->GetName()) == "X_wide_200_21") { // DEBUG
#else
    // get a good file name else skip:
    /* if no dimensions are specified, then use a key like:
     * quantity_position_z (e.g. phi_100_1)
     * */
    /* else use the normal key:
     * quantity_dimensions_position_z (e.g. phi_wide_100_1)
     * */
    bool isAGoodKey = false;
    if (this->dimensions == "") {
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
      if (this->dimensions == G4_PS_hist::AVAILABLE_DIMENSIONS[1]) {
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
      else if (this->dimensions == G4_PS_hist::AVAILABLE_DIMENSIONS[2]) {
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
                  << ". Could not find " << this->dimensions << "."
                  << std::endl;
        throw std::invalid_argument("Invalid dimension selection for tree key "
                                    "decoding. G4_PS_hist::read_root_file. ");
      }
    }
    if (isAGoodKey && Z != -1 &&
        setup ==
            std::stoi(
                this->configuration)) // full string can be found in the file
    {
#endif
      std::cout << "Processing root tree: " << key->GetName() << std::endl;
      // opening the histogram
      TH1D *hist1 = (TH1D *)myFile->Get(key->GetName());
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
      for (auto a = 0; a < root_hist_values.first.size(); a++) {
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
      for (auto a = 0; a < root_hist_values_trimmed.first.size(); a++) {
        std::cout << root_hist_values_trimmed.first[a] << ", "
                  << root_hist_values_trimmed.second[a] << std::endl;
      }
#endif

      if (root_hist_values_trimmed.first.size() == 0) {
        std::cout << "\t Warning: no valuable data points found for this "
                     "histogram: information is null."
                  << std::endl;
        continue;
      }

      // create a new binning
      uint N_points = G4_PS_hist::N_POINTS_MAX;
      if (root_hist_values_trimmed.first.size() <= N_points) {
        N_points = root_hist_values_trimmed.first.size();
        std::cout
            << "\tWarning: root points (" << N_points
            << ") are less then required points, adjusting to not oversample."
            << std::endl;
      }
      const double min_element =
          *std::min_element(std::begin(root_hist_values_trimmed.first),
                            std::end(root_hist_values_trimmed.first));
      const double max_element =
          *std::max_element(std::begin(root_hist_values_trimmed.first),
                            std::end(root_hist_values_trimmed.first));

      std::vector<double> new_binning;
      if (max_element == min_element) {
        new_binning.push_back(min_element);
        std::cout << "Single point detected, single bin used as " << min_element
                  << "." << std::endl;

      } else {
        std::cout << "Asking for linear binning with: [" << min_element << ", "
                  << max_element << "] with " << N_points << "-Points"
                  << std::endl;
        new_binning = this->linspace(min_element, max_element, N_points);
      }
#ifdef DEBUG2
      // DEBUG:
      std::cout << "Got a linspace:" << std::endl;
      for (auto a = 0; a < new_binning.size(); a++) {
        std::cout << a << ", " << new_binning[a] << std::endl;
      }

#endif

      std::vector<double> rebin_y_count =
          this->rebin_hisogram(root_hist_values_trimmed, new_binning);

#ifdef DEBUG2
      // DEBUG
      std::cout << "Rebinned histogram:" << std::endl;
      for (auto a = 0; a < rebin_y_count.size() - 1; a++) {
        std::cout << new_binning[a] << ", " << rebin_y_count[a] << std::endl;
      }
#endif

      std::string file_name_root =
          std::string(outputFolder) + std::string("hist_root_") + name +
          std::string("_") + std::to_string(setup) + std::string("_Z_") +
          std::to_string(Z) + std::string(".csv");

      this->write_hist_to_csv_file(file_name_root, root_hist_values_trimmed);

      std::string file_name_geant =
          std::string(outputFolder) + std::string("hist_geant_") + name +
          std::string("_") + std::to_string(setup) + std::string("_Z_") +
          std::to_string(Z) + std::string(".mac");

      double relative_abundance_at_given_z = 0;

      for (size_t a = 0; a < relative_abundance.first.size(); a++) {
        if (relative_abundance.first[a] == Z) {
          relative_abundance_at_given_z = relative_abundance.second[a];
          if (name == "Ek") {
#ifdef DEBUG
            std::cout << "Pushing to processed Z " << Z << std::endl;
#endif
            processed_z.push_back(Z);
          }
        }
      }

      G4_hist G4_h(file_name_geant, name,
                   "hist key: " + std::string(key->GetName()),
                   "Z = " + std::to_string(Z), "", Z, new_binning,
                   rebin_y_count, relative_abundance_at_given_z);

      // Set the half dimensions:
      // XXX: to check if this is ok sending new simulations
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
#ifdef DEBUG // FIXME: restore remove this
      std::cout << "Skipping file: " << std::string(key->GetName())
                << std::endl;
#endif
    }
  }
  return processed_z;
}

int G4_PS_hist::write_hist_to_csv_file(
    std::filesystem::path p,
    std::pair<std::vector<double>, std::vector<double>> hist_to_write) {
  std::ofstream file_to_write(p, std::ofstream::out);

  for (uint i = 0; i < hist_to_write.first.size(); i++) {
    file_to_write << hist_to_write.first[i] << ", " << hist_to_write.second[i]
                  << std::endl;
  }

  return 0;
}

std::vector<double> G4_PS_hist::rebin_hisogram(
    std::pair<std::vector<double>, std::vector<double>> h1_data,
    std::vector<double> bin_edges) {
  const uint bin_counts = bin_edges.size() + 1;

  if (std::abs(bin_edges.back() - h1_data.first.back()) >
      0.00001) { // prevent rounding error
    std::cerr << "Error in rebinning: the provided bin for the new histogram "
                 "is smaller than the original hisogram"
              << std::endl;
    std::cerr << "Last root historam bin = " << h1_data.first.back()
              << " while the last histogram bin for the new binning = "
              << bin_edges.back() << std::endl;
    std::cerr << "Absolute error: "
              << std::abs(h1_data.first.back() - bin_edges.back()) << std::endl;
    throw std::invalid_argument("Invalid binning for rebinning");
  }

  // initiazlize the vector
  std::vector<double> bin_count(bin_counts, 0.);
  const double precision = 1000;

  size_t j = 0; // index for the root binning
  size_t i = 0; // index for the G4 binning
  double new_bin_edge_max =
      std::round(bin_edges[i + 1] * precision) / precision;
  double new_bin_edge_min = std::round(bin_edges[i] * precision) / precision;
  do {
    /*
     * update the bin edge of the root histogram
     * */
    double binc = h1_data.second[j];
    double bin_edge_min = std::round(h1_data.first[j] * precision) / precision;
    double bin_edge_max =
        std::round(h1_data.first[j + 1] * precision) / precision;

    /*
     * check if it fails in the binning,
     * the bin is entirely contained in the new binning
     * */
#ifdef DEBUG2
    std::cout << "Bin edges: \t ro [" << bin_edge_min << ", " << bin_edge_max
              << "]" << std::endl;
    std::cout << "Bin edges: \t g4 [" << new_bin_edge_min << ", "
              << new_bin_edge_max << "]" << std::endl;
#endif
    if (bin_edge_min >= new_bin_edge_min && bin_edge_max <= new_bin_edge_max) {
      // no overap: sum to the bin content
      bin_count[i] += binc;
#ifdef DEBUG2
      std::cout << "Summing " << i << "-binc normally to " << bin_count[i]
                << std::endl;
#endif
      i++;
      new_bin_edge_max = std::round(bin_edges[i + 1] * precision) / precision;
      new_bin_edge_min = std::round(bin_edges[i] * precision) / precision;
    }
    /*
     * overlap the binning boundaries
     * */
    else if (bin_edge_min >= new_bin_edge_min &&
             bin_edge_max > new_bin_edge_max) {
      /*
       * do the ratio of the counts
       * */
      double const delta_x1 = new_bin_edge_max - bin_edge_min;
      double const delta_x2 = bin_edge_max - new_bin_edge_max;
      if (binc != 0) {
        bin_count[i] += binc * std::abs(delta_x1);
        bin_count[i + 1] += binc * std::abs(delta_x2);
      } else {
        bin_count[i] += 0;
      }
#ifdef DEBUG2
      std::cout << "Weighted sum " << i << "-binc = " << binc << " to "
                << std::endl;
      std::cout << "delta_x1 = " << delta_x1
                << " , bin_c =  " << binc * delta_x1 << std::endl;
      std::cout << "delta_x2 = " << delta_x2
                << " e bin_c =  " << binc * delta_x2 << std::endl;
      std::cout << "to bin[" << i << "] = " << bin_count[i] << ", bin[" << i + 1
                << "] = " << bin_count[i + 1] << std::endl;
#endif

      /*
       * update the bin edge on the geant4 histogram
       * */
      i++;
      new_bin_edge_max = std::round(bin_edges[i + 1] * precision) / precision;
      new_bin_edge_min = std::round(bin_edges[i] * precision) / precision;
    }

    // end of the cycle
    if (j == h1_data.first.size() - 3) {
#ifdef DEBUG2
      std::cout << "Summing " << i + 1 << "-binc (" << h1_data.second[j + 2]
                << ", j + 2 = " << j + 2 << ") normally to " << bin_count[i + 1]
                << " End of the cycle " << std::endl;
#endif
      bin_count[i + 1] += h1_data.second[j + 2];
    }

    j++;
  } while (j < h1_data.first.size() - 1 && i < bin_edges.size());

  return bin_count;
}

std::tuple<int, int, std::string> G4_PS_hist::parse_tree_key(std::string s) {
// DEBUG: print the string
#ifdef DEBUG2
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
#ifdef DEBUG2
  std::cout << "String type = " << type << std::endl;
#endif // DEBUG2

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
     * */

    std::regex re("_(\\d+)_(\\d+)"); // new

    std::smatch m;

    std::regex_search(s, m, re);
#ifdef DEBUG2
    std::cout << "Got a match; full = " << m[0] << std::endl;
    std::cout << "\t g1 = config  = " << m[1] << std::endl;
    std::cout << "\t g2 = Z = " << m[2] << std::endl;
#endif // DEBUG2

    int config, Z;
    try {
      config = std::stoi(m[1]);
      Z = std::stoi(m[2]);
    } catch (std::exception &e) {
      std::cerr << "ERROR: cannot convert tree name to extract config and Z"
                << std::endl;
      std::cerr << "\t tree name: " << s << std::endl;
    }

    return std::tuple(config, Z, type);

  } else {
    // if it is not the relative abundance histogram
    if (s != "Z_1") {
      std::cerr << "Cannot associate tree name to known class parameters;"
                << std::endl;
      std::cerr << "\t tree name '" << s << "' is being skipped." << std::endl;
    }
    return std::tuple(-1, -1, "Null");
  }
}

std::vector<double> G4_PS_hist::linspace(const double min, const double max,
                                         const uint N_elements) {
  std::vector<double> to_return;

  for (size_t i = 0; i < N_elements; i++) {
    to_return.push_back((max - min) * i / (N_elements - 1) + min);
  }
  return to_return;
}

std::pair<std::vector<double>, std::vector<double>>
G4_PS_hist::convert_root_hist_to_array(TH1D *hist) {
  std::vector<double> bin_right_edge;
  std::vector<double> bin_count;
  for (size_t i = 1; i <= hist->GetNbinsX(); i++) {
    double bin_right = hist->GetBinLowEdge(i);
    bin_right_edge.push_back(bin_right);
    bin_count.push_back((double)hist->GetBinContent(i - 1));
  }
  return std::make_pair(bin_right_edge, bin_count);
}

std::pair<std::vector<double>, std::vector<double>>
G4_PS_hist::trim_zeros_from_hist(
    std::pair<std::vector<double>, std::vector<double>> h1_data) {
  unsigned long long int i = h1_data.first.size(); // get the total size
  std::vector<double> counts = h1_data.second;
  std::vector<double> right_bin = h1_data.first;
  std::vector<double> counts_base = h1_data.second;
  std::vector<double> right_bin_base = h1_data.first;
  std::vector<double> trimmed_righ_bin;
  std::vector<double> trimmed_counts;

  /* upper removal */
  do {
    i--;
  } while (counts[i] == 0 || i <= 0);
  const unsigned long long int upperFirst = i;

  if (i == 0) {
    std::cout << "No valid points found" << std::endl;
    return std::make_pair(trimmed_righ_bin, trimmed_counts);
  }

  /* lower removal */
  i = 0;
  do {
    i++;
  } while (counts[i] == 0 || i == h1_data.first.size() - 1);

  /*
   * Not sure this is necessary as if there are only zeros, the first loop
   * should take care of this
   * */
  if (i == h1_data.first.size() - 1) {
    std::cout << "No valid points found" << std::endl;
    return std::make_pair(trimmed_righ_bin, trimmed_counts);
  }

  const unsigned long long int lowerFirst = i;

  for (auto j = lowerFirst; j <= upperFirst; j++) {
    trimmed_righ_bin.push_back(right_bin_base[j]);
    trimmed_counts.push_back(counts_base[j]);
  }

  return std::make_pair(trimmed_righ_bin, trimmed_counts);
}

bool G4_PS_hist::isWithinInterval(const double x_prec, const double x_succ,
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

bool G4_PS_hist::isUpperEdgeCloser(const double x_prec, const double x_succ,
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

double G4_PS_hist::get_hist_total_counts(TH1D *hist) {

  double hist_counts = 0;
  for (size_t i = 1; i <= hist->GetNbinsX(); i++) {
    hist_counts += hist->GetBinContent(i);
  }

  return hist_counts;
}
