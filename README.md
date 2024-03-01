# Apollonian-Prime
This code deals with prime components of Apollonian circle packings, and is used to produce data for the paper "Prime and thickened prime components in Apollonian circle packings", by Elena Fuchs, Holley Friedlander, Piper Harris, Catherine Hsu, James Rickards, Katherine Sanden, Damaris Schindler, and Katherine E. Stange (preprint forthcoming).

A prime component is a tangency-connected set of circles in a primitive Apollonian circle packing for which every curvature is an odd positive prime number. A thickened prime component is the union of a prime component with every circle tangent to the component.

## Main methods

* ```parthickenedmult```: parallel implementation to find the multipliticies of all curvatures in a thickened prime component in a given range. ```thickened``` is a (slower) non-parallel implementation of this.
* ```primerootred```: returns the prime root quadruple attached to v, or 0 if v contains no odd primes.
* ```primeroots_bin```: computes the number of prime roots, prime curvatures, and total curvatures in a given curvature range in a single packing. We ignore symmetry, so this method will not work on the strip packing.
* ```thickened_bin```: computes the number of (prime) curvatures in a thickened prime component, binning the data.

## Installation
System requirements:
* Linux/Mac: none that I am aware of;
* Windows: you need to be running PARI/GP through WSL. See this [tutorial](https://pari.math.u-bordeaux.fr/PDF/PARIwithWindows.pdf) for how to set this up.

In addition to this package, you need to download and compile the most up to date version of [Apollonian](https://github.com/JamesRickards-Canada/Apollonian) (follow the instructions listed there). This location can be anywhere, but you need to remember the directory ```D``` that you installed it in. In the directory for Apollonian-Prime (this project!), you need to create soft links to three of the Apollonian files as follows:
```
ln -s D/libapol-X-Y-Z.so
ln -s D/apol.gp
ln -s D/apol.h
```
If you update your copy of PARI/GP, you will need to re-build this project and re-link the new .so file (with the updated X-Y-Z version number).

Compiling this package follows the same process: you need to know where the version of PARI/GP you want to use is installed, in particular, the file ```pari.cfg```. The default location is ```/usr/local/lib/pari/pari.cfg```.
* If this is the correct location, call ```make``` to build the project.
* Otherwise, call ```make setup``` to search for the location of the file. By default the program searches in ```/usr```, but there is a chance it is not installed there (this sometimes happens on a server). If this is the case, you can supply an alternate location.
* If the program finds potential matches, it will ask you to confirm which files are correct, and saves them to ```paricfg_loc.txt```. Once this step completes, a call to ```make``` will compile the project! Modifying the program (e.g. via ```git pull```) won't require redoing this setup, unless the version of PARI/GP you use changes, and the location of this file changes with it.

At this point, you are all set! Call ```gp aprime``` to open gp and load the methods. ```?aprime``` and ```?apol``` accesses the help,
