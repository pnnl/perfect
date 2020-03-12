
disp (" *  ")
disp (" *  This script will evaluate the Lucas-Kanade benchmark results against the golden output.")
disp (" *  ")

arg_list = argv();

if (nargin < 1)
  disp (" * ERROR: invalid arguments")
  quit(-1)
endif

if ( strcmp(arg_list{1}, "SMALL") )
  disp (" *  Loading the golden file: inout/small_golden.mat")
  load ../../../inout/small_golden.mat;
elseif ( strcmp(arg_list{1}, "MEDIUM") )
  disp (" *  Loading the golden file: inout/medium_golden.mat")
  load ../../../inout/medium_golden.mat;
elseif ( strcmp(arg_list{1}, "LARGE") )
  disp (" *  Loading the golden file: inout/large_golden.mat")
  load ../../../inout/large_golden.mat;
endif

disp (" *  Loading the benchmark output: output.mat")
load output.mat;

snr = 20 * log10( norm(golden) / norm(golden-output) );

disp (" *")
fprintf(stdout, ' *  SNR: %.2f\n', snr);
disp (" *")

