============================================================================

SAR kernels overview

The source files within this directory implement a reference version
of the synthetic aperture radar kernels as described in the PERFECT
Suite Manual.

The three SAR kernels are split into three executables that all take as
a command-line argument the name of the directory containing the required
input files. 

The input files are available from the TAV team and named with the
following pattern:

    $(SIZE)_kernel1_input.bin
    $(SIZE)_kernel2_input.bin
    $(SIZE)_kernel3_input.bin
    $(SIZE)_golden_kernel1_output.bin
    $(SIZE)_golden_kernel2_output.bin
    $(SIZE)_golden_kernel3_output.bin

where the first three files are intputs to the respective kernels,
the last three files are golden outputs for the respective kernels,
and $(SIZE) is small, medium, or large.

If ENABLE_CORRECTNESS_CHECKING is #defined (it is by default),
then the kernel reference implementations further compute and report
the signal-to-noise correctness metric for the output generated
by each of the kernels.

If WRITE_OUTPUT_TO_DISK is #defined (it is by default), then
the kernel reference implementations write the output generated
by the kernel to $(SIZE)_kernelX_output.bin where X is 1, 2, or 3
based on the kernel being executed and $(SIZE) is small, medium,
or large.

Kernel 1 (SAR backprojection) utilizes OpenMP to reduce run-time
on multicore systems.  OpenMP can be removed from that kernel if
needed by removing the single "#pragma omp" from sar_backprojection.c.

============================================================================

Compiling and running the kernels

The kernels can be compiled as follows from the sar/kernels/ directory:

    make clean
    make INPUT_SIZE=LARGE

where INPUT_SIZE can be set to SMALL, MEDIUM, or LARGE for the
respective data sets.  The default if no INPUT_SIZE is specified is MEDIUM.

The kernels can then be run as follows:

    ./pfa-interp1/sar_kernel1_driver ../inout/
    ./pfa-interp2/sar_kernel2_driver ../inout/
    ./bp/sar_kernel3_driver ../inout/

Input and output file formats
============================================================================

In order to simplify the kernels, the necessary inputs are packed into a
single input file.  All of the fields within the input files are little
endian and the floating point values are stored in IEEE format.  Complex
data is stored in interleaved format (i.e., real and imaginary components
are adjacent in memory).  The constants references below (e.g., N_PULSES)
are defined in common/sar_params.h.

The file formats are as follows:

    input_kernel1.bin:

        Size (bytes)         | Format  | Description
        ------------------------------------------------------------------
        8*N_PULSES*N_RANGE   | float   | phase history data (complex, range-major order)
        8*N_PULSES           | double  | input start coordinates per pulse
        8*N_PULSES           | double  | input coordinate spacing per pulse
        8*PFA_NOUT_RANGE     | double  | output coordinates
        4*T_PFA              | float   | window (filter) weights

    input_kernel2.bin:

        Size (bytes)              | Format  | Description
        ------------------------------------------------------------------
        8*N_PULSES*PFA_NOUT_RANGE | float   | phase history data (complex, range-major order)
        8*N_PULSES*PFA_NOUT_RANGE | double  | input coordinates (pulse-major order)
        8*PFA_NOUT_AZIMUTH        | double  | output coordinates
        4*T_PFA                   | float   | window (filter) weights

    input_kernel3.bin :     

        Size (bytes)                   | Format  | Description
        ------------------------------------------------------------------
        8                              | double  | Carrier frequency, Hertz
        8                              | double  | R0 (distance to first range bin, meters)
        8                              | double  | dR (distance per range bin, meters)
        24*N_PULSES                    | double  | platform positions (x, y, z)
        8*N_PULSES*N_RANGE_UPSAMPLED   | float   | phase history data (complex, range-major order)

The output files include only the raw generated data and not parameters or
input coordinates.  All of the output files are complex, single-precision
floating point values consistent with the output of the associated kernel
(e.g., pixels for kernel 3 and resampled data for kernels 1 and 2).
