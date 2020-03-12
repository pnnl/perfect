-*-Mode: markdown;-*-

$HeadURL$
$Id$

PERFECT Suite
=============================================================================

The PERFECT Suite is a collection of kernels representing certain
domains of interest in radar and image processing.  The Suite has been
jointly developed by Pacific Northwest National Laboratory (PNNL) and
Georgia Tech Research Institute (GTRI) with support from the DARPA
PERFECT program.

* For more information on the PERFECT Suite:
  - https://hpc.pnnl.gov/PERFECT/
  - https://gitlab.pnnl.gov/perf-lab-hub/perfect


* Contacting: _firstname_._lastname_@pnnl.gov
  - Joseph Manzano (PNNL)
  - Nathan Tallent (PNNL)


* Contributors (alphabetical):
  - Kevin Barker (PNNL)
  - Thomas Benson (GTRI)
  - Dan Campbell (GTRI)
  - David Ediger (GTRI)
  - Nitin Gawande (PNNL)
  - Roberto Gioiosa (PNNL)
  - Adolfy Hoisie (PNNL)
  - Darren Kerbyson (PNNL)
  - Joseph Manzano (PNNL)
  - Andres Marquez (PNNL)
  - Leon Song (PNNL)
  - Nathan Tallent (PNNL)
  - Antonino Tumeo (PNNL)


* Reference:

  - Kevin Barker, Thomas Benson, Dan Campbell, David Ediger, Roberto Gioiosa, Adolfy Hoisie, Darren Kerbyson, Joseph Manzano, Andres Marquez, Leon Song, Nathan R. Tallent, and Antonino Tumeo.  PERFECT (Power Efficiency Revolution For Embedded Computing Technologies) Benchmark Suite Manual. Pacific Northwest National Laboratory and Georgia Tech Research Institute, December 2013.
  - doc/PERFECT-Suite-manual.pdf
  - https://gitlab.pnnl.gov/perf-lab-hub/perfect/perfect-suite/-/blob/master/doc/PERFECT-Suite-manual.pdf


Description:
=============================================================================

This release contains the following directories:
- doc/            : PERFECT Suite manual
- suite/pa1       : "PERFECT Application 1" application and kernels
- suite/stap/     : "Space-Time Adaptive Processing" application and kernels
- suite/sar/      : "Synthetic Aperture Radar" kernels
- suite/wami/     : "Wide Area Motion Imaging" application and kernels
- suite/required/ : Required kernels


Using:
=============================================================================

The make system is based on GNU make.  Type "make help" to see a list of
supported targets.

The make system defaults to GCC, but that can be changed by overriding
`CC` and `CFLAGS`, e.g.,
  ```sh
  # environment override
  CC=clang make
  ```
or
  ```sh
  # make override
  make CC=clang
  ```

To build and test the suite's kernels, issue the command
  ```sh
  make check
  ```
To manually select input sizes, override `INPUT_SIZE` with the values
`SMALL`, `MEDIUM`, or `LARGE`, e.g.,
  ```sh
  make INPUT_SIZE=SMALL check
  ```


