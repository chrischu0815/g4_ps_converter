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

#pragma once

/**
 * @file G4_PS_hist.hpp
 *
 * @brief Class reading a ROOT phase space file and generating the
 * corresponding Geant4 GPS macro files.
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

/**
 * @brief Class implementing the G4 phase space histogram converter.
 *
 * The class reads a ROOT file containing the phase space histograms
 * (Ekin, X, Y, Theta and Phi) recorded by the upstream simulation for every
 * particle, re-bins each histogram on a regular grid of at most
 * N_POINTS_MAX edges, and writes, for every particle:
 *
 *  - a CSV snapshot of the trimmed ROOT histogram
 *    (`hist_root_<type>_<configuration>_Z_<Z>.csv`),
 *  - one Geant4 GPS macro per quantity
 *    (`hist_geant_<type>_<configuration>_particle_<Z>.mac`).
 *
 * The main entry point is process_file(), which first computes the relative
 * particle abundances (from the Ekin histogram of the requested
 * configuration/dimension) and then generates all the files above.
 *
 * Two ROOT file layouts are supported, selected through the `dimensions`
 * constructor argument:
 *
 *  - **dimensioned files**: one histogram per particle, named
 *    `<type>_<dimensions>_<configuration>_<Z>`
 *    (e.g. `Ekin_wide_1500_2`);
 *  - **dimension-less files**: one histogram per particle, named
 *    `<type>_<configuration>_<Z>` (e.g. `Ekin_1_2`, as in
 *    `root_files/join.root`).
 *
 * The particle identifier is a charge Z (0..26) or one of the special
 * particle names listed in fspecialParticles.
 *
 * @note The relative abundances are taken from the Ekin histogram of the
 * selected configuration/dimension (see fabundance_root_hist_string) and sum
 * to 1; they are written as the intensity of the corresponding GPS source.
 *
 * */
class G4_PS_hist {

public:
  /**
   * @brief Constructor.
   *
   * Stores the input/output locations and the user selection, derives the
   * ROOT key from which the relative abundances are computed
   * (fabundance_root_hist_string) and the index used to select the beam
   * dimensions (index_of_dimensions), and creates the output directory if it
   * does not exist yet.
   *
   * The abundance histogram key is:
   *  - `Ekin_<configuration>` if dimensions is the empty string,
   *  - `Ekin_<dimensions>_<configuration>` otherwise.
   *
   * @param path_to_root_file path to the ROOT file to read
   * @param outputFolder output folder where the generated files are written
   * (created if missing)
   * @param configuration configuration string selecting the histograms to
   * read; one of AVAILABLE_CONFIGURATION
   * @param dimensions dimension of the scoring plane; one of
   * AVAILABLE_DIMENSIONS (the empty string selects the dimension-less file
   * layout)
   * @throw std::invalid_argument if REQUIRED_POINTS_PER_HIST <= 1
   * @throw std::runtime_error if the dimensions string is not one of the
   * recognized values or the output directory cannot be created
   * */
  G4_PS_hist(std::filesystem::path path_to_root_file,
             std::filesystem::path outputFolder,
             const std::string configuration, const std::string dimensions);

  /**
   * @brief Path to the ROOT file to be read.
   * */
  std::filesystem::path fpath_to_root_file;

  /**
   * @brief Output folder where all the generated files are written.
   *
   * Created by the constructor if it does not exist.
   * */
  std::filesystem::path foutputFolder;

  /**
   * @brief Configuration string selected by the user.
   *
   * One of AVAILABLE_CONFIGURATION; selects the set of histograms to read
   * from the ROOT file.
   * */
  std::string fconfiguration;

  /**
   * @brief Dimensions of the input scoring plane selected by the user.
   *
   * One of AVAILABLE_DIMENSIONS; the empty string means that the ROOT file
   * does not contain a dimension field in the tree keys.
   * */
  std::string fdimensions;

  /**
   * @brief Start the reading and conversion process.
   *
   * Calls calculate_relative_abundance() with the abundance histogram key
   * computed in the constructor, then read_root_file() to generate all the
   * CSV and GPS macro output files.
   *
   * @return the list of particle identifiers (Z or special particle name)
   * that have an Ekin histogram and have been processed
   * @note an empty result, or a result containing "-1", means that no
   * histogram matched the selected configuration/dimensions and nothing was
   * generated
   * */
  std::vector<std::string> process_file();

  /**
   * @brief Read the ROOT file and generate the output files.
   *
   * Opens the ROOT file, iterates over its keys, keeps the histograms whose
   * name matches the requested dimensions and configuration and whose
   * particle is not "Null", then for every kept histogram:
   *
   *  - converts it to (bin center, density) points
   *    (convert_root_hist_to_array()),
   *  - skips it if the total integrated content is
   *    <= REQUIRED_POINTS_PER_HIST,
   *  - trims the leading/trailing empty points (trim_zeros_from_hist()) and
   *    skips it if fewer than 2 points with data remain (a single point does
   *    not define a distribution for the GPS arbitrary histogram),
   *  - re-bins it on a regular grid of at most N_POINTS_MAX edges
   *    (linspace() + rebin_histogram()),
   *  - writes the trimmed ROOT histogram to a CSV file
   *    (write_hist_to_csv_file()) and the re-binned histogram to a Geant4 GPS
   *    macro (via G4_hist::write_g4_ps_file()).
   *
   * Non-TH1 objects are skipped with a warning. The first Ekin histogram
   * encountered uses /gps/source/intensity (the GPS default source 0), all
   * others use /gps/source/add.
   *
   * @return the list of particle identifiers that have an Ekin histogram and
   * have been processed; a single element "-1" if no Ekin histogram matched
   * @throw std::invalid_argument if the dimensions string cannot be decoded
   * against the tree keys
   * */
  std::vector<std::string> read_root_file();

  /**
   * @brief Write a (center, density) histogram to a CSV file.
   *
   * The CSV file contains two whitespace-separated columns: the first is the
   * bin center, the second the bin density (content per unit of width).
   *
   * @param p path of the CSV file to write (overwritten if it exists)
   * @param hist_to_write pair of (bin centers, bin densities) to write
   * @return 0 on success
   * */
  int write_hist_to_csv_file(
      std::filesystem::path,
      std::pair<std::vector<double>, std::vector<double>>);

  /// Maximum number of points (bin edges) allowed in the re-binned
  /// histogram, i.e. at most N_POINTS_MAX - 1 = 1023 bins.
  static const uint N_POINTS_MAX =
      1024; // Maxinum number of point allowed in the histogram
  //
  /// Minimum total integrated content (sum of the bin densities) a histogram
  /// must have to be accepted; histograms with less content are skipped.
  static const uint REQUIRED_POINTS_PER_HIST =
      2; // required points to generate the histogram
  // WARNING: cannot be zero! will rise an excption

  /// Available dimensions of the scoring plane, as they appear in the ROOT
  /// file tree keys. The empty string means that the file has no dimension
  /// field in the tree keys (e.g. root_files/join.root). The user selects an
  /// entry by its 1-based index (see the CLI documentation).
  static inline const std::vector<std::string> AVAILABLE_DIMENSIONS = {
      "", "wide", "very_wide"};

  /// Half-size (mm) of the beam used for the phase space, one entry per
  /// element of AVAILABLE_DIMENSIONS, i.e. the halfx/halfy written by the
  /// generated GPS macros: 75 mm beam for the dimension-less files, 100 mm
  /// for "wide", 200 mm for "very_wide".
  static inline const std::vector<double> SOURCE_DIMENSIONS = {

      37.5, 50., 100.}; // mm

  /// Available configurations, as they appear in the ROOT file tree keys:
  /// "100" (2m from the last element of the beamline),
  /// "1500" (1m from the target), "2000" (2m from the target),
  /// "1" (full-GCR style files, e.g. root_files/join.root).
  /// The user selects an entry by its 1-based index (see the CLI
  /// documentation).
  static inline const std::vector<std::string> AVAILABLE_CONFIGURATION = {
      "100", "1500", "2000", "1"};

  /**
   * @brief Relative abundance of the particles, computed by
   * calculate_relative_abundance().
   *
   * First element: particle identifiers (Z as string or special particle
   * name); second element: the corresponding relative abundance, normalized
   * so that all the values sum to 1.
   * */
  std::pair<std::vector<std::string>, std::vector<double>> frelative_abundance;

  /**
   * @brief List of special particle names (in string form) that can appear
   * as the particle identifier of a ROOT tree key, in addition to the charge
   * Z. The names are used both to match the tree keys
   * (parse_tree_key()) and as the Geant4 particle name (e.g. `e+`).
   * */
  inline static const std::vector<std::string> fspecialParticles = {
      "e-", "e\\+", "gamma", "mu-", "mu\\+", "pi-", "pi\\+"};

private:
  /// Access for the unit tests (see test/test_G4_PS_hist.cpp).
  friend class G4_PS_hist_Test;

  /**
   * @brief ROOT key (substring) of the histogram used to compute the
   * relative abundances.
   *
   * `Ekin_<configuration>` for the dimension-less layout,
   * `Ekin_<dimensions>_<configuration>` otherwise.
   * */
  std::string fabundance_root_hist_string;

  /**
   * @brief Substring identifying the kinetic energy histograms in the ROOT
   * file tree keys.
   * @note TODO: make this selectable.
   * */
  const std::string ek_string = "Ekin";

  /**
   * @brief Index of the selected dimensions inside AVAILABLE_DIMENSIONS /
   * SOURCE_DIMENSIONS.
   *
   * Set in the constructor: 0 for the empty string, 1 for "wide", 2 for
   * "very_wide". Used by read_root_file() to select the beam half-size
   * written in the generated GPS macros.
   * */
  int index_of_dimensions; // This index is used in the contructor to find which
                           // dimensions of the scoring to apply

  /// Substring identifying the X position histograms in the ROOT file tree
  /// keys.
  const std::string x_pos_string = "X";

  /// Substring identifying the Y position histograms in the ROOT file tree
  /// keys.
  const std::string y_pos_string = "Y";

  /// Substring identifying the Theta angle histograms in the ROOT file tree
  /// keys.
  const std::string theta_ang_string = "Theta";

  /// Substring identifying the Phi angle histograms in the ROOT file tree
  /// keys.
  const std::string phi_ang_string = "Phi";

  /**
   * @brief Parse a ROOT tree key and return (configuration, particle, type).
   *
   * The type is determined by the first recognized substring in the key
   * (Ekin, X, Y, Theta or Phi); the key is then matched against the pattern
   * `_<configuration>_(<particle>)`, where the particle is one of
   * fspecialParticles or a charge Z.
   *
   * Example keys: `Ekin_1_2` (dimension-less), `Ekin_wide_1500_2`
   * (dimensioned), `Ekin_1_e+` (special particle).
   *
   * @param s key name of the tree entry
   * @return tuple (configuration, particle, type):
   *  - configuration: integer parsed from the key, or -1 if not parseable
   *  - particle: charge Z as string, a special particle name, or "Null"
   *  - type: "Ek", "X", "Y", "Theta" or "Phi", or "Null" if the key is not a
   *    recognized histogram
   * */
  std::tuple<int, std::string, std::string> parse_tree_key(std::string s);

  /**
   * @brief Rebin the (center, density) points of an old histogram onto a new
   * grid.
   *
   * The old data are interpreted as a piecewise-linear pdf over the bin
   * centers; the integral of that pdf over every new bin (defined by
   * bin_edges) is computed exactly (trapezoids on each old segment) and
   * divided by the new bin width, so the total integral is conserved and
   * split bins are handled with the correct overlap fractions. Outside the
   * domain of the old pdf the density is 0. The scan is O(n_old + n_new).
   *
   * @param h1_data pair of (bin centers, bin densities) of the old histogram
   * @param bin_edges edges of the new histogram (size >= 2, sorted)
   * @return one mean density per new bin (size = bin_edges.size() - 1)
   * */
  std::vector<double>
      rebin_histogram(std::pair<std::vector<double>, std::vector<double>>,
                      std::vector<double>);

  /**
   * @brief Create a vector of N_elements equally spaced values in [min, max].
   *
   * Both endpoints are included, e.g. linspace(1, 10, 4) returns
   * {1, 4, 7, 10}.
   *
   * @param min minimum value a
   * @param max maximum value b
   * @param N_elements number of elements to create (N_elements < 2 gives an
   * empty vector)
   * @return a vector with N_elements equally spaced points in [min, max]
   * */
  std::vector<double> linspace(const double min, const double max,
                               const uint N_elements);

  /**
   * @brief Convert a ROOT histogram into a pair of arrays (bin center,
   * density).
   *
   * The density is the bin content divided by the bin width, which keeps the
   * representation correct also for non-uniform binning.
   *
   * @param hist ROOT histogram to be converted
   * @return pair whose first element is the bin center and second element
   * the bin density of each real bin (underflow/overflow excluded);
   * degenerate (zero-width) bins are skipped
   * */
  std::pair<std::vector<double>, std::vector<double>>
  convert_root_hist_to_array(TH1D *hist);

  /**
   * @brief Remove the leading and trailing zero-density points from a
   * (center, density) histogram.
   *
   * Makes the histogram more compact and better conditioned for the
   * rebinning. If every point is zero, an empty pair is returned.
   *
   * @param h1_data pair of (bin centers, bin densities) to trim
   * @return pair of (bin centers, bin densities) restricted to the first and
   * last points with data
   * */
  std::pair<std::vector<double>, std::vector<double>>
      trim_zeros_from_hist(std::pair<std::vector<double>, std::vector<double>>);

  /**
   * @brief Convert the ROOT histogram of the relative abundances into a pair
   * of (Z, counts) vectors.
   *
   * @param hist ROOT histogram containing the Z and their absolute abundance
   * @return pair of (Z, counts) contained in the histogram
   * @deprecated legacy declaration, no longer implemented or used; the
   * abundances are computed by calculate_relative_abundance() from the Ekin
   * histograms instead.
   * */
  std::pair<std::vector<int>, std::vector<int>>
  convert_root_hist_relative_abundance(TH1D *hist);

  /**
   * @brief Get the total content (sum of the bin contents) of a ROOT
   * histogram.
   *
   * @param hist ROOT histogram from which the integral is calculated
   * @return the sum (total counts) of that histogram, including over-under
   * flows
   * */
  double get_hist_total_counts(TH1D *hist);

  /**
   * @brief Get the integral of the a ROOT histogram.
   *
   * @param hist ROOT histogram from which the integral is calculated
   * @return the integral (total counts) of that histogram
   * */
  double get_hist_integral(TH1D *hist);

  /**
   * @brief Calculate the relative abundance of the particles.
   *
   * Iterates over the keys of the ROOT file and keeps every histogram whose
   * name contains tree_to_use. The relative abundance of each particle is
   * its histogram total content divided by the sum of the contents of all
   * the kept histograms, so the resulting abundances sum to 1. The result is
   * stored in frelative_abundance.
   *
   * @param tree_to_use substring identifying the abundance histograms to use
   * (fabundance_root_hist_string)
   * @throw std::runtime_error if the input file cannot be opened or the
   * total content of the abundance histograms is zero
   * */
  void calculate_relative_abundance(const std::string tree_to_use);

  /**
   * @brief Tell whether a reference value is inside the half-open interval
   * [x_prec, x_succ).
   *
   * @param x_prec lower value of the binning to be compared
   * @param x_succ next value in sequence of the binning
   * @param ref reference value, cast to double for the comparison
   * @return true if ref \in [x_prec, x_succ), else false
   * @note currently unused (legacy helper).
   * */
  bool isWithinInterval(const double x_prec, const double x_succ,
                        const int ref);

  /**
   * @brief Tell whether the lower edge of the interval [x_prec, x_succ) is
   * closer to a reference value than the upper one.
   *
   * @param x_prec lower value of the binning to be compared
   * @param x_succ next value in sequence of the binning
   * @param ref reference value, cast to double for the comparison
   * @return true if |ref - x_prec| < |ref - x_succ|, else false
   * @note currently unused (legacy helper).
   * */
  bool isUpperEdgeCloser(const double x_prec, const double x_succ,
                         const int ref);
};
