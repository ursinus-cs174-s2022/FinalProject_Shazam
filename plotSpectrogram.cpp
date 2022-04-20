#include "dsp.h"
#include "dspviz.h"
#include "simplecanvas/simplecanvas.h"
#include "AudioFile.h"
#include <stdio.h>
#include <stdlib.h>

struct AudioInfo {
    const char* audiopath = "ExampleQueries/baddayClean.wav";
    const char* impath = "out.png";
    int win = 1024; // Window length of STFT
    int hop = 512; // Hop length of STFT
    int maxBin = 128; // Maximum bin of STFT to display (<= win/2+1)
    int scale = 4; // Factor by which to scale up pixels when displaying spectrogram
};

/**
 * @brief Parse the command line arguments that specify parameters
 *        for computing and saving the spectrogram
 * 
 * @param argc Number of command line arguments
 * @param argv Array of command line arguments
 * @return AudioInfo 
 */
AudioInfo parseArgs(int argc, char** argv) {
    AudioInfo info;
    argv++, argc--;
    while (argc > 0) {
        if ((*argv)[0] == '-') {
            if (strcmp(*argv, "--help") == 0) {
                printf("Usage: ./plotSpectrogram --audiopath <audio file path> --impath <path to output image> --win <Window length> --hop <Hop length> -- maxBin <Maximum number of bins to use>\n");
                exit(0);
            }
            else if (strcmp(*argv, "--audiopath") == 0) {
                argv++; argc--;
                info.audiopath = (const char*)*argv;
            }
            else if (!strcmp(*argv, "--impath")) {
                argv++; argc--;
                info.impath = (const char*)*argv;
            }
            else if (!strcmp(*argv, "--win")) {
                argv++; argc--;
                info.win = atoi(*argv);
            }
            else if (!strcmp(*argv, "--hop")) {
                argv++; argc--;
                info.hop = atoi(*argv);
            }
            else if (!strcmp(*argv, "--maxBin")) {
                argv++; argc--;
                info.maxBin = atoi(*argv);
            }
            else if (!strcmp(*argv, "--scale")) {
                argv++; argc--;
                info.scale = atoi(*argv);
            }
            else { 
                fprintf(stderr, "Invalid option: %s\n", *argv);
            }
        }
        argv++, argc--; 
    }
    return info;
}

int main(int argc, char** argv) {
    AudioInfo info = parseArgs(argc, argv);

    // Load audio samples
    AudioFile<double> file;
    file.load(info.audiopath);
    int N;
    double* samples = extractAudioSamples(file, &N);
    
    // Create spectrogram
    DSP dsp(info.win);
    int NWin;
    double** S = dsp.spectrogram(samples, N, info.win, info.hop, true, &NWin);

    // Plot the spectrogram and save to an image
    SimpleCanvas canvas = plotSpectrogram(S, info.maxBin, NWin, info.scale, info.scale);
    canvas.write(info.impath);

    // Clean up memory
    deleteSpectrogram(S, info.win);
    delete[] samples;
    return 0;
}