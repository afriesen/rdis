# rdis
C++ code for the RDIS algorithm from ["Recursive Decomposition for Nonconvex Optimization." Friesen and Domingos, IJCAI 2015](http://homes.cs.washington.edu/~pedrod/papers/ijcai15.pdf).

Main platform this was tested on: Mac OS X 10.9 64 bit
Other platforms this hopefully works on: Mac OS X 32 bit, Linux 32 bit & 64 bit

Binaries are included for Mac OS X 64 bit (built on 10.9 Yosemite) in rdis/bin/macosx64, otherwise you will have to build them (see below).

To build (if all dependencies are accessible from default system paths), simply call the following:

`cmake .`

`make -j8`

RDIS builds in <BASE_DIR>/build and there should be three executables: testRDIS, optBA, and optSinusoid.
- testRDIS: run the debug tests (if no parameters are specified on the command line), or optimize a polynomial specified in a file.
  e.g., run `./build/testRDIS`
  e.g., run `./build/testRDIS -f data/testpoly.txt`
  - see the notes in data/testpoly.txt for writing your own polynomial for RDIS to optimize
- optBA: optimize a bundle adjustment function; 
  e.g., run `./build/optBA -f data/ladybug-problem-49-7776-pre.txt --ncams 5 --npts 30`
- optSinusoid: optimize a high-dimensional sinusoid
  e.g., run `./build/optSinusoid`
- calling `./build/<executable> --help` will print the command line options

## Dependencies:
- NOTE: cmake vars `EXT_INCL_DIR` and `EXT_LIB_DIR` specify paths in which the build should look for the following libraries (if they're enabled) except for boost, instead of specifying each one for each library.
- Boost 1.55.0 (required, not included) (or potentially later, but 1.55.0 is tested) from [http://www.boost.org/](http://www.boost.org/)
  - Required boost libs: system, filesystem, program_options, and chrono.
  - You must build these yourself. If the paths to the boost include folder and libs are non-standard (e.g., not in /usr/local/include or /usr/include), then you can specify the paths when calling cmake as `cmake -D BOOST_ROOT=/path/to/boost .` or you can specify `BOOST_INCLUDEDIR` and `BOOST_LIBRARYDIR` separately.
    - Building these is easy; simply follow the Getting Started instructions at [www.boost.org](http://www.boost.org/doc/libs/release/more/getting_started/index.html).
- PaToH (not required, included, used by default) hypergraph partitioning library from [http://bmi.osu.edu/umit/software.html](http://bmi.osu.edu/umit/software.html). 
  - Pre-built libs are included in the repo already, so this should (ideally) just work.
  - Can be disabled when calling cmake with `cmake -D USE_PATOH=FALSE .`
  - You can change the search path by setting the cmake variable `PATOH_DIR`.
- HMetis (not required, not included) hypergraph partitioning library from [http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview). 
  - Nothing is included, but once you download it you can specify the path in cmake variable `HMETIS_DIR`.
- levmar (not required, not included) Levenberg-Marquardt optimizer from [http://users.ics.forth.gr/~lourakis/levmar/](http://users.ics.forth.gr/~lourakis/levmar/). 
  - Nothing is included, but once you download it you can specify the path in cmake variable `LEVMAR_DIR`.
  - Note that this library requires that [LAPACK](http://www.netlib.org/lapack/) (and [BLAS](http://www.netlib.org/blas/)) are installed and work on your system.
  
## Known Issues:
- If the minimum is on the border of the domain, RDIS will do poorly because it does not handle constraints properly (a non-infinite domain constitutes box constraints on the function).
