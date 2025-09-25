# Geant4 Phase Space Converter
#### Program to generate a Geant4 Phase space from root output files

#### Conversion from ROOT to G4 data format ("G4_PS_converter")

The data analysis consists in the following steps:

1. Reading the ROOT file and import the data of the histograms ("Ekin", "X", "Y", "Theta", "Phi")

2. Rebinning the root histograms with a new binning using the whole range, from 0 to 1024 if possible.

3. Writing the G4 phase space using the GPS macro (.mac) format which includes information on kinetic energy, X-Y positions, and angular distributions of all particles (with charge, Z, from 0 to 26) of the radiation field to be simulated in G4.


### How to run:

Call the program with:
```bash
G4_PS_creator <inputFile> <outputPath> <parameters>
```

- `<inputFile>` is the file containing the phase space in root format. A sample of such file is provided in `root_files/join.root`.

- `<outputPath>` is the output path where the geant4 macro will be generated.

The output contains a `sim.mac` that is the main macro that will call all the others insuccession.

## Available phase space:

The available phase space is in the folder `./macros/generated/`. 
A documentation on how to use the generated phase space is also available in `./docs/` folder.

### Versions and license:

- ROOT `6.30/04`, *LGPL 2.1*, [license](https://root.cern/about/license/)

- Geant4 `geant4-v11.2.0`, *Geant4 license*, [license](https://geant4.org/download/license)
