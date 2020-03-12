WAMI Kernel Benchmarks README
===

The WAMI kernels can be compiled for three different input sizes (SMALL,
MEDIUM, and LARGE).  To compile:

  $ make INPUT_SIZE={SMALL,MEDIUM,LARGE}

To run the benchmark kernel and validate the output automatically:

  $ make check
