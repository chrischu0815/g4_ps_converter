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

/**
 * @file G4_hist.hpp
 *
 * @brief Class hosting a single re-binned phase space histogram and writing
 * its Geant4 GPS macro, plus the CommandGenerator that assembles the main
 * sim.mac file.
 *
 */

/**
 * @brief Class containing a single re-binned phase space histogram.
 *
 * The class hosts one histogram of the phase space (one quantity, e.g. Ek or
 * X, of one particle) and knows how to write it as a Geant4 General Particle
 * Source (GPS) macro file:
 *
 *  - `Ek` -> arbitrary kinetic energy distribution (in GeV per particle),
 *  - `X`/`Y` -> biased transverse position distribution over the source
 *    plane,
 *  - `Theta`/`Phi` -> user angular distribution (in radians).
 *
 * The histogram values are the re-binned contents produced by
 * G4_PS_hist::rebin_histogram() on a grid of N+1 edges (fright_bin_edge)
 * with N points (fbin_count).
 * */
class G4_hist {
public:
  /**
   * @brief Constructor.
   *
   * Stores the output location, the histogram data and the metadata; parses
   * the particle identifier as a charge Z (fZ) when possible.
   *
   * @note TODO: replace the Z-based identification with an index mapped to
   * the particle types/names.
   *
   * @param path_to_file_to_write path of the GPS macro file to write
   * @param hist_type type of the histogram: one of "Ek", "X", "Y", "Theta",
   * "Phi" (selects the GPS commands emitted by write_g4_ps_file())
   * @param name short name of the histogram, written as a comment header in
   * the macro
   * @param description description written as a comment header in the macro
   * @param comments comments written as a comment header in the macro
   * @param particle particle identifier: the charge Z as string (e.g. "2")
   * or a special particle name (e.g. "e+"); fZ is -1 for the latter
   * @param right_bin_edge bin edges of the re-binned histogram: N+1 edges
   * for N points
   * @param bin_count histogram value per bin (density, i.e. content per unit
   * of width); the written macro is normalized from these values
   * @param relative_intensity relative intensity of the particle (Ek macro
   * only): the values over all particles sum to 1 and are written as the
   * GPS source intensity
   * @param is_first_source true for the first particle of the phase space:
   * its source is configured with /gps/source/intensity (the GPS default
   * source 0 is a 1 MeV geantino), all other particles use
   * /gps/source/add
   * */
  G4_hist(std::filesystem::path path_to_file_to_write, std::string hist_type,
          std::string name, std::string description, std::string comments,
          std::string particle, std::vector<double> right_bin_edge,
          std::vector<double> bin_count, double relative_intensity,
          bool is_first_source = false);

  /**
   * @brief Print the histogram data (index, bin center, value) on stdout.
   *
   * @return 1 on success
   * */
  int print_data();

  /**
   * @brief Save the histogram data (bin center, value per line) to the
   * output file as CSV.
   *
   * Used only for debugging (see the DEBUG flag); the file is overwritten
   * if it exists.
   *
   * @return 1 on success
   * */
  int save_data_as_csv();

  /**
   * @brief Write the class content to a Geant4 GPS macro file.
   *
   * The file starts with the name/description/comments as comment lines,
   * then the GPS commands of the histogram type:
   *
   *  - **Ek**: the source intensity (/gps/source/intensity for the first
   *    particle, /gps/source/add otherwise), the particle definition
   *    (define_source()), a square source plane of half-size
   *    dimHalfX/dimHalfY centered at (0, 0, -beam_offset) mm along
   *    (0, 0, 1), a planar (perfectly collimated) angular reset, and the
   *    arbitrary energy distribution: E[GeV] = E[MeV/u] * A / 1000 with
   *    linear interpolation, the bin values normalized to sum to 1.
   *  - **X** / **Y**: biased position distribution (biasx/biasy). GPS maps
   *    a sampled value v to x = v * 2 * halfx - halfx, so a data position x
   *    is written as v = (x + halfx) / (2 * halfx); the bin values are
   *    normalized to their maximum.
   *  - **Theta**: user angular distribution; sets /gps/ang/type user and
   *    /gps/hist/type theta. Values are in radians (fed directly to the
   *    sin/cos of the GPS), normalized to sum to 1; values < 1e-10 are
   *    written as 0.
   *  - **Phi**: /gps/hist/type phi with the same conventions as Theta;
   *    /gps/ang/type user is set by the Theta macro, so the two macros must
   *    be used together.
   *
   * If the file already exists, it is overwritten.
   *
   * @return 0 on success, 1 if the file cannot be opened
   * @throw std::runtime_error if the histogram is empty (undefined
   * normalization) or no isotope is defined for the particle Z
   * */
  int write_g4_ps_file();

  /**
   * @brief Mass number A of the most abundant stable isotope for each
   * charge Z.
   *
   * The phase space neglects the isotope contribution (see docs/readme.md),
   * so a single A per Z is used to scale the per-nucleon energy (MeV/u) to
   * the total kinetic energy per particle. Z = 43 (Tc) and Z = 61 (Pm) have
   * no stable isotope: the standard reference isotopes 99 and 147 are used.
   * */
  std::map<int, int> associated_A = {
      {1, 1},    {2, 4},    {3, 7},    {4, 9},    {5, 11},   {6, 12},
      {7, 14},   {8, 16},   {9, 19},   {10, 20},  {11, 23},  {12, 24},
      {13, 27},  {14, 28},  {15, 31},  {16, 32},  {17, 35},  {18, 40},
      {19, 39},  {20, 40},  {21, 45},  {22, 48},  {23, 51},  {24, 52},
      {25, 55},  {26, 56},  {27, 59},  {28, 58},  {29, 63},  {30, 64},
      {31, 69},  {32, 74},  {33, 75},  {34, 80},  {35, 79},  {36, 84},
      {37, 85},  {38, 88},  {39, 89},  {40, 91},  {41, 93},  {42, 96},
      {43, 99},  {44, 98},  {45, 103}, {46, 106}, {47, 107}, {48, 114},
      {49, 115}, {50, 120}, {51, 121}, {52, 130}, {53, 127}, {54, 132},
      {55, 133}, {56, 138}, {57, 139}, {58, 140}, {59, 141}, {60, 142},
      {61, 147}, {62, 152}, {63, 153}, {64, 158}, {65, 159}, {66, 164},
      {67, 165}, {68, 167}, {69, 169}, {70, 170}, {71, 175}, {72, 180},
      {73, 181}, {74, 184}, {75, 187}, {76, 192}, {77, 193}, {78, 195},
      {79, 197}, {80, 202}, {81, 204}, {82, 208}};

  /**
   * @brief Intensity per charge Z.
   *
   * @note legacy placeholder, currently unused: the source intensity is
   * taken from the relative abundances computed by
   * G4_PS_hist::calculate_relative_abundance().
   * */
  std::map<int, int> source_intensity = {
      {0, 1},  {1, 1},  {2, 1},  {3, 1},  {4, 1},  {5, 1},  {6, 1},  {7, 1},
      {8, 1},  {9, 1},  {10, 1}, {11, 1}, {12, 1}, {13, 1}, {14, 1}, {15, 1},
      {16, 1}, {17, 1}, {18, 1}, {19, 1}, {20, 1}, {21, 1}, {22, 1}, {23, 1},
      {24, 1}, {25, 1}, {26, 1}, {27, 1}, {28, 1}};

  /**
   * @brief Half-size (mm) of the source plane in the X direction, written as
   * /gps/pos/halfx.
   *
   * Must match the half-size of the scoring surface; set by the caller
   * (G4_PS_hist::read_root_file()) from
   * G4_PS_hist::SOURCE_DIMENSIONS. The default (200 mm) is kept as a
   * fallback for direct use of the class.
   * */
  double dimHalfX = 200; // mm

  /**
   * @brief Half-size (mm) of the source plane in the Y direction, written as
   * /gps/pos/halfy.
   *
   * See dimHalfX.
   * */
  double dimHalfY = 200; // mm

private:
  /// Access for the unit tests (see test/test_G4_hist.cpp).
  friend class G4_hist_Test;

  /// Path of the GPS macro file to write.
  std::filesystem::path fpath_to_file_to_write;

  /// Type of the histogram: "Ek", "X", "Y", "Theta" or "Phi".
  std::string fhist_type;

  /// Short name of the histogram (comment header of the macro).
  std::string fname;

  /// Description of the histogram (comment header of the macro).
  std::string fdescription;

  /// Comments of the histogram (comment header of the macro).
  std::string fcomments;

  /// Particle identifier as in the ROOT file: charge Z as string or special
  /// particle name.
  std::string fparticle;

  /// Charge Z parsed from fparticle; -1 for special particle names.
  int fZ;

  /// Bin edges of the re-binned histogram: N+1 edges for N points.
  std::vector<double> fright_bin_edge;

  /// Histogram value per bin (density, i.e. content per unit of width).
  std::vector<double> fbin_count;

  /// Relative intensity of the particle (written as the GPS source
  /// intensity for the Ek macro).
  double frelative_intensity;

  /// true for the first particle of the phase space: Geant4 GPS starts with
  /// a default source 0 (a 1 MeV geantino, intensity 1.0), so the first
  /// particle must configure it with /gps/source/intensity instead of adding
  /// a new source with /gps/source/add.
  bool fIsFirstSource;

  /**
   * @brief Beam offset (mm) from the origin.
   *
   * The source plane is centered at (0, 0, -beam_offset) mm, i.e. 200 mm
   * upstream of the origin, in the generated phase space file.
   * */
  const std::string beam_offset = "200";

  /// @name Beam generator dimensions
  /// @{
  /// NOTE: dimHalfX/dimHalfY (public) must match the scoring surface
  /// dimensions.
  /// @}

  /**
   * @brief Write the particle definition commands of the source.
   *
   * Emits, depending on the charge Z:
   *  - special particle (Z = -1): /gps/particle <name>,
   *  - Z = 0: /gps/particle neutron,
   *  - Z = 1: /gps/particle proton,
   *  - Z >= 2: /gps/particle ion + /gps/ion <Z> <A>, where A is the most
   *    abundant stable isotope (associated_A) and the nuclear mass in u
   *    (G4NucleiProperties) is written as a comment.
   *
   * @param to_write open file stream to write on
   * @param particle particle identifier (special particle name if Z = -1)
   * @param Z charge Z, -1 for special particles
   * @throw std::runtime_error if no isotope is defined for the charge Z
   * */
  void define_source(std::ofstream &to_write, std::string particle, int Z);
};

/**
 * @brief Generate the sim.mac file containing all the sourcing.
 *
 * This file will execute all the macro files that are generated. Executing
 * this file in Geant4 loads all the generated GPS macros.
 * */
class CommandGenerator {
public:
  /**
   * @brief Generate the sim.mac file containing all the sourcing.
   *
   * Writes a header with the configuration and the scoring dimensions, then
   * for every particle:
   *
   *  - the Ek macro is **required**: it creates the source. If it is missing
   *    the whole particle is skipped (with a warning), because its
   *    X/Y/Theta/Phi macros would silently configure the previous
   *    particle's source;
   *  - the X and Y macros are optional: a missing macro is warned about and
   *    the corresponding position falls back to the GPS default (uniform
   *    over the plane);
   *  - the Theta and Phi macros must be used **together** (they form the GPS
   *    user angular distribution): if only one of them exists, both are
   *    skipped (with a warning) and the source stays perfectly collimated
   *    (see the /gps/ang/type planar reset in the Ek macro).
   *
   * The existence checks and the /control/execute lines use names relative
   * to the directory of the command file, i.e. the generated macros must be
   * in the same folder as sim.mac.
   *
   * @param commandFileName file name to open and write on
   * @param scoring_dimensions_mm half-size (mm) of the scoring plane,
   * written in the file header as comment
   * @param configuration defines which configuration is used and
   * consequently set up the file name of the sourced files
   * @param particles list of the particles (Z or special name) generated in
   * the macros. The sourcing is emitted only for the macro files that
   * actually exist on disk.
   * */
  static void generateCommandFile(const std::string &commandFileName,
                                  const double scoring_dimensions_mm,
                                  const std::string &configuration,
                                  const std::vector<std::string> particle);
};
