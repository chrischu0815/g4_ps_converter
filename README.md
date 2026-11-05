# Geant4 Phase Space Converter

#### Program to generate a Geant4 Phase space from root output files

#### Conversion from ROOT to G4 data format ("G4_PS_converter")

The data analysis consists in the following steps:

1. Reading the ROOT file and import the data of the histograms ("Ekin", "X", "Y", "Theta", "Phi")

2. Rebinning the root histograms with a new binning using the whole range, from 0 to 1024 if possible.

3. Writing the G4 phase space using the GPS macro (.mac) format which includes information on kinetic energy, X-Y positions, and angular distributions of all particles (with charge, Z, from 0 to 26) of the radiation field to be simulated in G4.

## Requirements

- CMake >= 3.10 (3.10 is the minimum for `gtest_discover_tests`; required for a warning-free configure on CMake 4.x), with a C++17 compiler (GCC or Clang)
- ROOT >= 6.30 with TMVA (on this machine ROOT 6.40.04 is installed system-wide; the `ro` alias activates the ROOT environment)
- Geant4 >= 11.2 (only `G4NucleiProperties` is used, for isotope/nuclear data)
- GoogleTest (optional: the unit tests are skipped if it is not found)
- Doxygen (optional: to regenerate the API documentation)

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

A `Debug` build enables `-g -fsanitize=address,undefined` for local debugging.

The compiled sources in `lib/` are built as the static library `g4_ps_lib`; the resulting executable is `build/G4_PS_creator`.

## How to run

Call the program with:

```bash
G4_PS_creator <inputFile> <outputPath> [configuration] [dimensions]
```

- `<inputFile>` is the file containing the phase space in root format. A working sample is provided in `root_files/Mod1000_1_1e5.root` (`root_files/join.root` is kept as a legacy example of the older key naming). The file should contain the histograms `EkinU` (kinetic energy in MeV/u), `X` and `Y` (position in mm, plane orthogonal to the beam direction z), `Theta` and `Phi` (spherical coordinates, radians).

- `<outputPath>` is the output path where the geant4 macros will be generated (created if missing). The simulation has to be run from this directory, using `sim.mac` as the Geant4 macro.

- `[configuration]` / `[dimensions]` are optional 1-based indices selecting the user configuration and the scoring-plane dimensions. Run `G4_PS_creator` without arguments to list the available selections (defaults: configuration `100`, normal scoring plane).

Example:

```bash
./build/G4_PS_creator root_files/join.root /tmp/PS_out/
```

The output contains a `sim.mac` that is the main macro that will call all the others in succession, plus one `hist_geant_*` macro and one `hist_root_*.csv` per generated histogram.

## Unit tests

61 unit tests (GoogleTest, discovered by CTest) cover the histogram rebinning, the ROOT-to-array conversion, and the GPS macro generation:

```bash
cd build
ctest --output-on-failure
```

The tests live in `test/` and build against the same `g4_ps_lib` static library as the executable.

## Documentation

Doxygen-generated API documentation of the public classes is produced as HTML in `docs/doxygen/html/` (git-ignored; open `index.html`). Regenerate it with:

```bash
doxygen Doxyfile
```

## Continuous integration

`.github/workflows/ci.yml` configures and builds the project, runs the unit tests, and performs a smoke conversion of `root_files/Mod1000_1_1e5.root` on every push and pull request. The workflow targets a self-hosted Linux runner with ROOT and Geant4 preinstalled (see the comments in the file for the GitHub-hosted alternative).

## Helper scripts

- `bash/run.sh` — end-to-end GSI workflow: convert a batch of ROOT files, run the Geant4 simulation, and validate the re-sampled phase space with `python/check_PS.py`. The paths inside must be adapted to the local setup.
- `python/` — analysis helpers (`check_PS.py`, normalization scripts).

## Available phase space

The available phase space is in the folder `./macros/generated/`.
A documentation on how to use the generated phase space is also available in [./docs/readme.md](./docs/readme.md).

## Recent improvements

- **Fixed a bug in `rebin_histogram`** (`lib/G4_PS_hist.cpp`): old segments straddling a new bin edge were silently dropped, zeroing out part of the distribution. As a consequence, some histograms were left with fewer than 2 usable points and skipped entirely. The fix restores them: the bundled `Mod1000_1_1e5.root` sample now converts to 169 histograms instead of 133.
- **Unit test suite**: 61 GoogleTest/CTest tests in `test/`, integrated through CMake (`enable_testing()` + `gtest_discover_tests`).
- **CMake restructure**: the sources in `lib/` now compile into the static library `g4_ps_lib`, shared by the executable and the tests; C++17 is now required.
- **Doxygen documentation** of the full public API (zero warnings), HTML output in `docs/`.

## Versions and license

- ROOT `6.30/04`+, *LGPL 2.1*, [license](https://root.cern/about/license/)

- Geant4 `11.4.3`, *Geant4 license*, [license](https://geant4.org/download/license)
