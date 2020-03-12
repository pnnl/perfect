PA1 Kernel Benchmarks README
===

The PA1 kernels can be compiled for three different input sizes (SMALL,
MEDIUM, and LARGE).  To compile:

  $ make INPUT_SIZE={SMALL,MEDIUM,LARGE}

The binary executable must be run from a directory where ../input/ is
valid and contains the input files.  To run the benchmark kernel and
validate the output automatically:

  $ make check

Input images for the PA1 kernels can be found in the input/ directory.
Images are grayscale with an 8-bit color depth.  (Note that the kernel
benchmark specifications require up to 16 bits of dynamic range.)
Three input sizes are provided:

  - SMALL: 640 x 480 pixels
  - MEDIUM: 1920 x 1080 pixels
  - LARGE: 3840 x 2160 pixels

JPEG versions are given for reference.  Kernel reference implementations
use input files in Octave matrix format.  To visualize an input image
in Octave, the following syntax may be used:

  octave:1> load input_small.mat 
  octave:2> input = uint8(input);
  octave:3> imshow(input)

Golden output files are provided in the output/ directory for
comparison. The kernels process each input frame 30 times and produce
30 output frames in Octave matrix format.  The output frames can be
visualized in a similar manner, however the kernel specifications use
up to 16 bits of dynamic range, so conversion may be necessary:

  octave:1> load histeq_output.small.mat 
  octave:2> output = output / 256;
  octave:3> output = uint8(output);
  octave:4> imshow(output)

