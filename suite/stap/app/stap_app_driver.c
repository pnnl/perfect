/* -*-Mode: C;-*- */

/**BeginCopyright************************************************************
 *
 * $HeadURL$
 * $Id$
 *
 *---------------------------------------------------------------------------
 * Part of PERFECT Benchmark Suite (hpc.pnnl.gov/projects/PERFECT/)
 *---------------------------------------------------------------------------
 *
 * Copyright ((c)) 2014, Battelle Memorial Institute
 * Copyright ((c)) 2014, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * 1. Battelle Memorial Institute (hereinafter Battelle) and Georgia Tech
 *    Research Corporation (GTRC) hereby grant permission to any person
 *    or entity lawfully obtaining a copy of this software and associated
 *    documentation files (hereinafter "the Software") to redistribute
 *    and use the Software in source and binary forms, with or without
 *    modification.  Such person or entity may use, copy, modify, merge,
 *    publish, distribute, sublicense, and/or sell copies of the
 *    Software, and may permit others to do so, subject to the following
 *    conditions:
 * 
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimers.
 * 
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 * 
 *    * Other than as used herein, neither the name Battelle Memorial
 *      Institute nor Battelle may be used in any form whatsoever without
 *      the express written consent of Battelle.
 * 
 *      Other than as used herein, neither the name Georgia Tech Research
 *      Corporation nor GTRC may not be used in any form whatsoever
 *      without the express written consent of GTRC.
 * 
 *    * Redistributions of the software in any form, and publications
 *      based on work performed using the software should include the
 *      following citation as a reference:
 * 
 *      Kevin Barker, Thomas Benson, Dan Campbell, David Ediger, Roberto
 *      Gioiosa, Adolfy Hoisie, Darren Kerbyson, Joseph Manzano, Andres
 *      Marquez, Leon Song, Nathan R. Tallent, and Antonino Tumeo.
 *      PERFECT (Power Efficiency Revolution For Embedded Computing
 *      Technologies) Benchmark Suite Manual. Pacific Northwest National
 *      Laboratory and Georgia Tech Research Institute, December 2013.
 *      http://hpc.pnnl.gov/projects/PERFECT/
 *
 * 2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 *    BATTELLE, GTRC, OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *    INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **EndCopyright*************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "stap_params.h"
#include "stap_utils.h"

#include "corner_turn.h"
#include "fft-1d.h"
#include "stap_covariance_estimate.h"
#include "stap_system_solver.h"
#include "stap_apply_weighting.h"

#define ENABLE_CORRECTNESS_CHECKING

/* N_DOP is 256 for all of our scenarios */
#define N_DOP_LOG2 (8)

#if INPUT_SIZE == INPUT_SIZE_SMALL
    static const char *input_filename = "small_app_input.bin";
    static const char *output_filename = "small_app_output.bin";
    static const char *steering_vector_filename = "small_app_steering_vectors.bin";
#elif INPUT_SIZE == INPUT_SIZE_MEDIUM
    static const char *input_filename = "medium_app_input.bin";
    static const char *output_filename = "medium_app_output.bin";
    static const char *steering_vector_filename = "medium_app_steering_vectors.bin";
#elif INPUT_SIZE == INPUT_SIZE_LARGE
    static const char *input_filename = "large_app_input.bin";
    static const char *output_filename = "large_app_output.bin";
    static const char *steering_vector_filename = "large_app_steering_vectors.bin";
#else
    #error "Unhandled value for INPUT_SIZE"
#endif

int main(int argc, char **argv)
{
    complex (*datacube)[N_PULSES][N_RANGE] = NULL;
    complex (*datacube_pulse_major_padded)[N_RANGE][N_DOP] = NULL;
    complex (*doppler_datacube)[N_DOP][N_RANGE] = NULL;
    complex (*covariances)[N_BLOCKS][N_CHAN*TDOF][N_CHAN*TDOF] = NULL;
    complex (*cholesky_factors)[N_BLOCKS][N_CHAN*TDOF][N_CHAN*TDOF] = NULL;
    complex (*adaptive_weights)[N_BLOCKS][N_STEERING][N_CHAN*TDOF] = NULL;
    complex (*steering_vectors)[N_CHAN*TDOF] = NULL;
    complex (*output)[N_DOP][N_RANGE] = NULL;

    char *input_directory = NULL;

    const size_t num_datacube_elements = N_CHAN * N_PULSES * N_RANGE;
    const size_t num_datacube_pulse_major_padded_elements = N_CHAN *
        N_RANGE * N_DOP;
    const size_t num_doppler_datacube_elements = N_CHAN * N_DOP * N_RANGE;
    const size_t num_covariance_elements = (TDOF*N_CHAN) * (TDOF*N_CHAN) *
        N_DOP * N_BLOCKS;
    const size_t num_adaptive_weight_elements = N_DOP * N_BLOCKS *
        N_STEERING * (N_CHAN*TDOF);
    const size_t num_steering_vector_elements = N_STEERING *
        (N_CHAN*TDOF);
    const size_t num_output_elements = N_STEERING *
        N_DOP * N_RANGE;

    int chan, range;

#ifdef ENABLE_CORRECTNESS_CHECKING
    complex (*gold_output)[N_DOP][N_RANGE] = NULL;
#endif

    if (argc != 2)
    {
        fprintf(stderr, "%s <directory-containing-input-files>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    input_directory = argv[1];

    datacube = XMALLOC(sizeof(complex) * num_datacube_elements);
    datacube_pulse_major_padded = XMALLOC(sizeof(complex) *
        num_datacube_pulse_major_padded_elements);
    doppler_datacube = XMALLOC(sizeof(complex) * num_doppler_datacube_elements);
    covariances = XMALLOC(sizeof(complex) * num_covariance_elements);
    cholesky_factors = XMALLOC(sizeof(complex) * num_covariance_elements);
    adaptive_weights = XMALLOC(sizeof(complex) * num_adaptive_weight_elements);
    steering_vectors = XMALLOC(sizeof(complex) * num_steering_vector_elements);
    output = XMALLOC(sizeof(complex) * num_output_elements);

    memset(datacube_pulse_major_padded, 0, sizeof(complex) *
        num_datacube_pulse_major_padded_elements);
    memset(doppler_datacube, 0, sizeof(complex) * num_doppler_datacube_elements);
    memset(covariances, 0, sizeof(complex) * num_covariance_elements);
    memset(cholesky_factors, 0, sizeof(complex) * num_covariance_elements);
    memset(adaptive_weights, 0, sizeof(complex) * num_adaptive_weight_elements);
    memset(steering_vectors, 0, sizeof(complex) * num_steering_vector_elements);
    memset(output, 0, sizeof(complex) * num_output_elements);

#ifdef ENABLE_CORRECTNESS_CHECKING
    gold_output = XMALLOC(sizeof(complex) * num_output_elements);
#endif

    printf("Loading input files from %s...\n", input_directory);

    read_complex_data_file(
        (complex *) datacube,
        input_filename,
        input_directory,
        num_datacube_elements);

    read_complex_data_file(
        (complex *) steering_vectors,
        steering_vector_filename,
        input_directory,
        num_steering_vector_elements);

#ifdef ENABLE_CORRECTNESS_CHECKING
    read_complex_data_file(
        (complex *) gold_output,
        output_filename,
        input_directory,
        num_output_elements);
#endif

    /*
     * Re-order the data from [channel, pulse, range] order (i.e. range-major order)
     * to [channel, range, pulse] order (i.e. pulse-major order) for the FFT.  This
     * function also pads with zeros in the pulse dimension from N_PULSES to N_DOP.
     */
    corner_turn_to_pulse_major_order_and_pad(
        datacube_pulse_major_padded,
        datacube);

    printf("Applying FFTs to convert to Doppler space...\n");
    for (chan = 0; chan < N_CHAN; ++chan)
    {
        for (range = 0; range < N_RANGE; ++range)
        {
            /* Apply a 1D FFT to convert to Doppler space */
            int rc = fft((float *) &datacube_pulse_major_padded[chan][range][0],
                N_DOP, N_DOP_LOG2, FFT_FORWARD);
            if (rc != 0)
            {
                fprintf(stderr, "Error: %s:%u: FFT returned error code %d.\n",
                    __FILE__, __LINE__, rc);
                exit(EXIT_FAILURE);
            }
        }
    }

    /*
     * Re-order the data back to range-major order so that it has proper dimensions
     * for the following covariance estimation. Note that the datacube can be
     * re-ordered in different ways (and the kernel implementations appropriately
     * updated as an optimization), but using range-major keeps the full application
     * kernels identical to the standalone kernels in terms of interface.
     */
    corner_turn_to_range_major_order(
        doppler_datacube,
        datacube_pulse_major_padded);
    
    /* Kernel 1: Covariance Estimation (Outer Products) */
    printf("Estimating covariance (via outer products)\n");
    stap_compute_covariance_estimate(
        covariances,
        doppler_datacube);

    /* Kernel 2: Weight Generation (Linear System Solver) */
    printf("Computing adaptive weights (via a linear system solver)\n");
    stap_system_solver(
        adaptive_weights,
        covariances,
        steering_vectors,
        cholesky_factors);

    /* Kernel 3: Adaptive Weighting (Inner Products) */
    printf("Applying adaptive weighting (via inner products)\n");
    stap_apply_weighting(
        output,
        doppler_datacube,
        adaptive_weights,
        steering_vectors);

#ifdef ENABLE_CORRECTNESS_CHECKING
    {
        double snr;
        printf("Computing correctness metrics...\n");
        snr = calculate_snr(
            (complex *) gold_output,
            (complex *) output,
            num_output_elements);
        printf("\tSNR after STAP application : %.2f dB\n", snr);
    }
    FREE_AND_NULL(gold_output);
#endif

    FREE_AND_NULL(datacube);
    FREE_AND_NULL(datacube_pulse_major_padded);
    FREE_AND_NULL(doppler_datacube);
    FREE_AND_NULL(covariances);
    FREE_AND_NULL(cholesky_factors);
    FREE_AND_NULL(adaptive_weights);
    FREE_AND_NULL(steering_vectors);
    FREE_AND_NULL(output);

    return 0;
}
