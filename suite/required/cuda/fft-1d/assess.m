
disp (" *  ")
disp (" *  This script will evaluate the FFT benchmark results against the internal fft() function.")
disp (" *  ")

disp (" *  Loading the input file: random_input.mat")
load random_input.mat

disp (" *  Loading the benchmark output: fft_output.mat")
load fft_output.mat

snr = 20*log10(norm(fft(input)) / norm(fft(input)-output));

disp (" *")
fprintf(stdout, ' *  SNR: %.2f\n', snr);
disp (" *")
