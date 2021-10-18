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
                1    tank with air (default)
                2    tank with oil

-z   z position of the photodetection plane [mm]
                200 (default)

-beamx    x shift for the laser beam position [mm]
                0 (default)

-beamy    y shift for the laser beam position [mm]
                0 (default)

-s    Gaussian smearing of the laser beam (sigma) [mm]
                1 (default)

-d    Divergence of the laser beam [rad]
                0 (default)
		
-l    focusing system
                0    no lens
                1    2-layer spherical lens
                2    2-layer cylindrical lens
                3    3-layer spherical lens (default)
                4    1-layer spherical lens (with air-gap)
                5    test lens
                6    3-layer cylindrical lens
                9    2-layer sperical achromat lens
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

tank dimensions etc. can be configured in src/PrtDetectorConstruction.cxx

```

## Installation
```
Geant4 and root should be installed and initialized.

git clone https://github.com/rdom/prttools
git clone https://github.com/rdom/prtlens

cd prtlens
mkdir build
cd build
cmake ..
make -j4

```

## Example of event display
```
./prtlens -g 1 -l 3 -theta 10 -m "3 eV" -e 100 -s 1 -z 300 -b 0
```
![alt text](https://github.com/rdom/prtlens/raw/master/pics/example_1.png)
![alt text](https://github.com/rdom/prtlens/raw/master/pics/example_2.png)

## Example of script usage from macro folder
```
root loadlib.C drawHP.C


seq 200 5 400 | xargs -n1 -P4 -I{} ../build/prtlens -l 1 -theta 10 -m "2 eV" -e 10000 -s 2 -z {} -b 1 -o hits_{}.root

seq 200 5 400 | xargs -n1 -P4 -I{} root -b -q loadlib.C fit_sim.C"(\"hits_{}.root\")"

hdd out.root hits*out.root

```