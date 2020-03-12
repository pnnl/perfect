-*-Mode: text;-*-

$HeadURL$
$Id$

=============================================================================
=============================================================================

------------------------------------------------------------
Revision History
------------------------------------------------------------

----------------------------------------
Version 0.2.0, 2014.05.30
----------------------------------------

- Release serial applications and specifications for PA1 (Perfect
  Application 1), STAP (Space-Time Adaptive Processing), and WAMI
  (Wide Area Motion Imaging).  (We do not expect to release an
  application for SAR.)

- Release parallel (OpenMP) versions of all kernels.

----------------------------------------
Version 0.1.1, 2014.04.01
----------------------------------------

- Public release

----------------------------------------
Version 0.1.0 (Pre-release), 2013.11.08
----------------------------------------

- Release Stage 4 kernels and specification: Wide Area Motion Imaging (WAMI)

----------------------------------------
Version 0.0.5 (Pre-release), 2013.09.30
----------------------------------------

- Release Stage 3 kernels: Perfect Application 1 (PA1)
- Specify input ranges for all kernels (SAR and PA1)
- Specify kernel-centric work metrics for all kernels

----------------------------------------
Version 0.0.4 (Pre-release), 2013.09.06
----------------------------------------

- Release Stage 2 kernels: Synthetic Aperture Radar (SAR)
- Specify input ranges for STAP and Required kernels
- Specify kernel-centric work metrics
- Common make system with common targets
- Bug fix: sort

----------------------------------------
Version 0.0.3 (Pre-release), 2013.08.16
----------------------------------------

- Expected revision of Required kernels
  - Use (GNU)Make instead of CMake
  - FFT kernels no longer depend on FFTW package

----------------------------------------
Version 0.0.1/0.0.2 (Pre-release), 2013.07.71
----------------------------------------

- Release Stage 1 kernels:
  - Required kernels: FFT-1D, FFT-2D, Sort
  - Space-Time Adaptive Processing (STAP)
