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

/**
 * @file G4_PS_hist.hpp
 *
 * @brief Class containing histogram information
 *
 * @author enrico.pierobon
 * Contact: e.pierobon@gsi.de
 *
 */

// #include "TFile.h"
#include "G4_hist.hpp"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TTree.h" //TMVA namespace
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
/*
 *
 * Class implementing the G4PS histogram
 *
 * */

class G4_PS_hist {

public:
  G4_PS_hist(std::filesystem::path path_to_root_file,
             std::filesystem::path outputFolder,
             const std::string configuration, const std::string dimensions);
  /*
   * class that will handle root files and re-sample them
   * @params path_to_root_file: root file to open
   * @params outputFolder: output folder where the output will be generated
   * @param configuration: configuration string that defines the root hist (and
   * their configuration) 0 associated to 100 mm from the last element of the
   * beamline 1 associated to 2 m from the last element of the beamline
   * */

  std::filesystem::path path_to_root_file; // path to the file to be readed

  std::filesystem::path outputFolder;

  std::string configuration; // choose the configuration 1: 2m away from the
                             // last element of the beamline
  std::string dimensions;    // the dimensions of the input scoring
  /*
   * @brief start the reading and conversion process
   *
   * Call the functions to calculate the relative abundance and call the
   * function to process the root trees and generate an output files as intended
   *
   * @return the list of Z that has been processed
   * If -1 is returned the execution has falied
   * */
  std::vector<int> process_file();

  /*
   * @brief read root file and import the data
   *
   * Open a root file and import the data of the histogram
   *
   *
   * @param none
   * @return vector of procesed Z
   *
   * */
  std::vector<int> read_root_file();

  /*
   * @brief write the root histogram as a csv file:
   * the CSV file will contain as first column the left-bin edge
   * the second column is the actual bin counts
   *
   * return 1 if successfull
   */
  int write_hist_to_csv_file(
      std::filesystem::path,
      std::pair<std::vector<double>, std::vector<double>>);

  static const uint N_POINTS_MAX =
      1024; // Maxinum number of point allowed in the histogram
  //
  // static const uint REQUIRED_POINTS_PER_HIST = 100;
  static const uint REQUIRED_POINTS_PER_HIST =
      0; // required points to generate the histogram

  static inline const std::vector<std::string> AVAILABLE_DIMENSIONS = {
      "", "wide", "very_wide"};
  static inline const std::vector<double> SOURCE_DIMENSIONS = {

      37.5, 50., 100.}; // mm

  static inline const std::vector<std::string> AVAILABLE_CONFIGURATION = {
      "100", "1500", "2000"};

  std::pair<std::vector<unsigned int>, std::vector<double>> relative_abundance;

private:
  // define which key to select to calculate the relative abundance in the PS
  std::string abundance_root_hist_string;

  const std::string ek_string = "Ekin"; // substring on root files telling that
                                        // kinet energy i contained in the data
                                        // TODO: make this selectable
  //
  int index_of_dimensions; // This index is used in the contructor to find which
                           // dimensions of the scoring to apply

  const std::string x_pos_string = "X"; // substring on root files telling that
  const std::string y_pos_string = "Y"; // substring on root files telling that
  const std::string theta_ang_string =
      "Theta"; // substring on root files telling that
  const std::string phi_ang_string =
      "Phi"; // substring on root files telling that
  /*
   * @brief rebin the root histogram with a new binning using the whole range:
   * from 0 to 1024
   *
   * Giving the root histogram, create a new histogram with a new binning
   * the histogram new binning is provided as parameter
   *
   * @param hist: root histogram to be rebinned
   * @param new_binning: new binning as an array of double with maximum size of
   * 1024. the binning elements are the left margin of the binning
   * @return return the histogram counts
   * */
  std::array<double, N_POINTS_MAX>
      rebin_histogram(TH1D, std::array<double, N_POINTS_MAX>);

  /** @brief check if the key name coming from the tree and return its type
   *
   *  Check, depending on the sring S, which kind of particle it is.
   *  return the charge of the particle
   *  @param string s as the key name of the tree entry
   *  @return charge of the particle extracted from the string
   */
  // std::pair<int, int> parse_tree_key(std::string s);

  /** @brief check if the key name coming from the tree and return its type
   *
   *  Check, depending on the sring S, which kind of particle it is.
   *  return the charge of the particle and the type of the scoring
   *  @param string s as the key name of the tree entry
   *  @return charge of the particle extracted from the string
   */
  std::tuple<int, int, std::string> parse_tree_key(std::string s);

  /* @brief function to rebin the histogram

    This function rebin the hisogram splitting the couts on overlapping bins eg.
    | ...c1... | ...c2... | ...c3... |    binning1 new (geant4)
    | .... | .... | .... | .... |         binning2 old (root)
    |<-x1->|<-x2->|<-x3->|<-x4->|
    |......|xxx|                          delta11
    |..........|xx|                       delta12

    ci is the count of the i-th binning of the new (geant4) histogram
    di is the count of the i-th binning of the old (root) histogam (not shown)

    yi is the binwidth of the i-th bin of the new histogram (not shown)
    xi is the binwidth of the i-th bin of the old histogram

    then c1 = d1 + d2*delta11
    where delta11 = y1 - x1 // right intersecation of the bins

    c2 = d2 *(delta12)
    where delta12 = x1 + x2 - y1 // left intersecation of the bins
    and so on.

    In the case of bin #2 new, the entire bin #3 old is contained therefore the
  whole counts will be added.


  * @param std::pair<std::array<double, N>, <std::array<double, N>> old
  histogram to be rebinned
  * @param std::array<double, N> bin_edges of the new histogram
  *
  * @return std::array<double, N> new_bin_counts
  */
  std::vector<double>
      rebin_hisogram(std::pair<std::vector<double>, std::vector<double>>,
                     std::vector<double>);

  /*
   * @brief create a series of equally spaced elements:
   *
   * Creates an array of equally spaced elements [0, 5, 10, 15, 20, ...] \in [a,
   b]
   * givin the minimum, the maximum, the number of elements is fixed and defined
   * @param const double min: minium elements a
   * @param const double max: maximum element b
   * @param const uint: number of elements required
   * note: the elemets are included
   *
   * use like this:
      std::array<double, 30>  = this->linspace<30>(1, 10);
    creates an array of 30 elements in the range [1, 10]

   * @return an array with equally spaced points in [a, b]
   * */
  std::vector<double> linspace(const double min, const double max,
                               const uint N_elements);

  /* @brief convert a root histogram in a pair of arrays: bin left edge, conts
   @tparam SIZE_N size of the histogram
   @param hist root histogram to be converted
   @return pair containing as first the left bin edge, and as second the counts
   of that bin
  */
  std::pair<std::vector<double>, std::vector<double>>
  convert_root_hist_to_array(TH1D *hist);

  /*
   * @brief remove zeroes from the histrogam points to make it more compact and
   prone to rebinning

   starting from the end of the array, it removes zeroes untill the first non
   zero element is found
   @return a pair constaingin as first the left bin edge and as second the
   counts of that bin

  */
  std::pair<std::vector<double>, std::vector<double>>
      trim_zeros_from_hist(std::pair<std::vector<double>, std::vector<double>>);

  /* @brief convert the root hist of the relative abundance in pair os Z, counts
   * @param TH1D root hist containing the Z and their absolute abundance
   * @return a pair of (Z, counts) contained in the histogram
   * */
  std::pair<std::vector<int>, std::vector<int>>
  convert_root_hist_relative_abundance(TH1D *hist);

  /* @brief get the integral of the histogram to be used as relativa abundance
   * @param TH1D root hist from where the integral is calculated
   * @retrurn the integeal of that histogram
   * */
  double get_hist_total_counts(TH1D *hist);

  /* @brief calculate the relative abundange of the given root tree of the class
   * This function will iterate over z, and according to the selected variable
   * it will select these trees and calcualte the relative abundante
   * @param tree_to_use is the substring that the function will use to idenfity
   * which tree to use
   * @return return a  pair <Z, counts> of the relative abundance
   * note this will require normalization */
  void calculate_relative_abundance(const std::string tree_to_use);

  /*
   * @brief tell the reference is inside the interval between x_prec and x_succ:
   * @param x_prec (double) containing the x_prec value of binning to be
   * compared
   * @param x_succ (double) containing the next value in sequence of the binning
   * @param ref (int) containing the reference value to be used as
   * comparison. If will be casted as double for comparison.
   * @return true if ref \in [x_prec, x_succ) else return false
   * */
  bool isWithinInterval(const double x_prec, const double x_succ,
                        const int ref);

  bool isUpperEdgeCloser(const double x_prec, const double x_succ,
                         const int ref);
};
