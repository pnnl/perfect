============================================================================

STAP kernels overview

The source files within this directory implement a reference version
of the space-time adaptive processing (STAP) kernels as described in
the PERFECT Suite Manual.

The three STAP kernels are split into three executables that all take as
a command-line argument the name of the directory containing the required
input files. 

The input files are available from the TAV team.  Each input data
set includes the following files:

    $(SIZE)_input.bin -- input to the first kernel (radar data cube)
    $(SIZE)_steering_vectors.bin -- input steering vectors, as described
        in the suite documentation
    $(SIZE)_kernel1_output.bin -- reference output data for kernel 1;
        can also be used as input to kernel 2
    $(SIZE)_kernel2_output.bin -- reference output data for kernel 2;
        can also be used as input to kernel 3
    $(SIZE)_kernel3_output.bin -- reference output data for kernel 3

where $(SIZE) is either small, medium, or large.  The input set
size used by the kernels depends upon the size specified during
compilation.  When compiling, the command line argument
INPUT_SIZE can be specified as either SMALL, MEDIUM, or
LARGE to choose between the different sample cases.  The default
sample case is MEDIUM.  For example, the following builds the
kernels for the LARGE input case and subsequently runs all of
the kernels.

    make clean
    make INPUT_SIZE=LARGE
    make check

If ENABLE_CORRECTNESS_CHECKING is #defined (it is by default),
then the reference implementation further computes and reports
the signal-to-noise correctness metric for the output generated
by each of the kernels.
