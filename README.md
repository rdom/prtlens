## Synopsis
```
prtlens [OPTION] [ARGUMENT] ... [OPTION] [ARGUMENT]

example:
./prtlens -l 3 -theta 0 -m "2 eV" -e 100 -s 1 -z 300

```
## Options
```

-o    output file name
      focalplane.root (default)

-g    geometry configuration
                1    tank with air
                2    tank with oil (default)

-z   z position of the photodetection plane [mm]
                200 (default)

-s    Gaussian smearing of the laser beam (sigma) [mm]
                1 (default)

-l    focusing system
                0    no lens
                1    2-layer spherical lens
                2    2-layer cylindrical lens
                3    3-layer spherical lens (default)
                4    1-layer spherical lens (with air-gap)
                5    test lens
                6    3-layer cylindrical lens		
                10   ideal lens (thickness = 0, ideal focusing)

-theta    polar angle between particle beam and bar radiator [deg]
                0 (default)

-phi  azimuth angle between particle beam and bar radiator [deg]
                0 (default)

-e    number of simulated photons
                1 (default)

-m    particle momentum [GeV/c]
                "2 eV" (default)

-seed    seed number for the random generator 

-b    batch mode
               1    run silent (without GUI)

```

## Installation
```
Geant4 and root should be installed and initialized.

git clone https://github.com/rdom/prtdirc
git clone https://github.com/rdom/prtlens

cd prtlens
mkdir build
cd build
cmake ..
make -j4

```

## Example of script usage from macro folder
```
root loadlib.C drawHP.C


seq 200 5 400 | xargs -n1 -P4 -I{} ../build/prtlens -l 1 -theta 10 -m "2 eV" -e 10000 -s 2 -z {} -b 1 -o hits_{}.root

seq 200 5 400 | xargs -n1 -P4 -I{} root -b -q loadlib.C fit_sim.C"(\"hits_{}.root\")"

hdd out.root hits*out.root

```