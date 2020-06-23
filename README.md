-*-Mode: markdown;-*-

$Id$


PERFECT Suite
=============================================================================

**Home**:
  - https://hpc.pnnl.gov/PERFECT/
  - https://gitlab.pnnl.gov/perf-lab-hub/perfect


**About**: The PERFECT Suite is comprised of kernels and applications
representing domains of interest in radar and image processing. The
suite was formed through a survey of applications used by all teams as
well as consultations with numerous application domain experts. Each
kernel in its own right is computationally important and provides
reasonable coverage of a domain of interest. Furthermore, each kernel
is algorithmically interesting.

The Suite has been jointly developed by Pacific Northwest National
Laboratory (PNNL) and Georgia Tech Research Institute (GTRI) with
support from the DARPA PERFECT program.


**Citation**:
  > Kevin Barker, Thomas Benson, Dan Campbell, David Ediger, Roberto Gioiosa, Adolfy Hoisie, Darren Kerbyson, Joseph Manzano, Andres Marquez, Leon Song, Nathan R. Tallent, and Antonino Tumeo.  "PERFECT (Power Efficiency Revolution For Embedded Computing Technologies) Benchmark Suite Manual." Pacific Northwest National Laboratory and Georgia Tech Research Institute, December 2013. https://hpc.pnnl.gov/PERFECT/


**Contacts**: (_firstname_._lastname_@pnnl.gov)
  - Joseph Manzano
  - Nathan R. Tallent


**Contributors**: (alphabetical)
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



Details
=============================================================================

**Manual**:
  - doc/PERFECT-Suite-manual.pdf (https://gitlab.pnnl.gov/perf-lab-hub/perfect/perfect-suite/-/blob/master/doc/PERFECT-Suite-manual.pdf)
  - See above for citation.


**Contents**: The PERFECT Suite contains the following kernels and
applications.  There are serial C and OpenMP versions for all; in
select cases there are CUDA variants.

- `suite/pa1` : "PERFECT Application 1" kernels and application.
  - Discrete Wavelet Transform
  - 2D Convolution
  - Histogram Equalization

- `suite/stap`  : "Space-Time Adaptive Processing" kernels and application.
  - System Solver
  - Inner Product
  - Outer Product

- `suite/sar` : "Synthetic Aperture Radar" kernels.
  - Interpolation 1
  - Interpolation 2
  - Back Projection (Non-Fourier SAR)

- `suite/wami` : "Wide Area Motion Imaging" kernels and application.
  - Debayer
  - Image Registration
  - Change Detection

- `suite/required` : Required (common) kernels
  - Sort
  - FFT 1D
  - FFT 2D


Using
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

