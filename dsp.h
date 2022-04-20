/**
 * Programmer: Chris Tralie
 * Purpose: Vanilla C++ implementation of the Fast-Fourier Transform
 * and double-Time Fourier Transform
 */
#ifndef DSP_H
#define DSP_H

#include <complex>
#include <iostream>
#include <math.h>

using namespace std;

typedef complex<double> cdouble;
#define FFT_FORWARD 0
#define FFT_INVERSE 1
#define PI 3.1415926535897932384626433832795

class DSP {
	private:
		cdouble*** W;//Cache the complex coefficients for the FFT
		double* hannwindow;
		int fftsize;

		/**
		 * Initialize the complex coefficients for the FFT
		 * @param fftsize The length of the fft
		 */
		void initCoeffs(int fftsize);

		/**
		 * Initialize the coefficients in Hann window
		 * @param fftsize The length of the fft
		 */
		void initWindow(int fftsize);

		/**
		 * Perform an in-place Cooley-Tukey FFT
		 * @param toReturn Array that holds FFT coefficients
		 * @param N Length of array (assumed to be power of 2)
		 * @param inverse Whether this is a forward or inverse FFT
		 */
		void performfft(cdouble* toReturn, int N, int inverse);
	
	public:
		cdouble fftres;
		DSP(int fftsize);
		~DSP();

		/**
		 * Implement the dft directly from the definition (used for speed comparison)
		 * @param sig Complex signal on which to compute dft
		 * @param toReturn The array that will hold the fourier coefficients
		 * @param N Length of signal
		 * @return Complex DFT coefficients
		 */
		void dft(cdouble* sig, cdouble* toReturn, int N);

		/**
		 * Perform the FFT on a complex signal
		 * @param sig The signal
		 * @param toReturn The array that will hold the fourier coefficients
		 * @param N Length of the signal (assumed to be power of 2)
		 * @return An N-length array with FFT coefficients
		 */
		void fft(cdouble* sig, cdouble* toReturn, int N);
	
		/**
		 * Perform the inverse FFT on an array of complex FFT coefficients
		 * @param sig The FFT coefficients
		 * @param toReturn The array that will hold the complex time series
		 * @param N Length of the FFT coefficients (assumed to be power of 2)
		 * @return An N-length array with FFT coefficients
		 */
		void ifft(cdouble* sig, cdouble* toReturn, int N);
		
		/**
		 * Helper function to create a complex array out of an array of 
		 * real amplitude samples
		 * @param data An array of doubles for the audio data
		 * @param res Array holding the result
		 * @param N Total number of samples in data
		 * @param start Index to start in the array
		 * @param win Length of the window
		 * @param useWindow Whether to use the window
		 */
		void toWindowedComplexArray(double* data, cdouble* res, int N, int start, int win, bool useWindow);

		/**
		 * Perform a double-time fourier transform on a bunch of samples
		 * @param sig Samples in the signal
		 * @param N Length of signal
		 * @param win Window length
		 * @param hop Hop length
		 * @param useWindow Whether to use the window
		 * @param NWin Number of windows (returned by reference)
		 * @return An NWin x win 2D array of complex doubles
		 */
		cdouble** stft(double* sig, int N, int win, int hop, bool useWindow, int* NWin);

		/**
		 * Perform a magnitude double-time fourier transform on a bunch of samples
		 * @param sig Samples in the signal
		 * @param N Length of signal
		 * @param win Window length
		 * @param hop Hop length
		 * @param useWindow Whether to use the window
		 * @param NWin Number of windows (returned by reference)
		 * @return A win x NWin 2D array of complex doubles
		 */
		double** spectrogram(double* sig, int N, int win, int hop, bool useWindow, int* NWin);
};

/**
 * Free the memory associated to an STFT
 * @param S STFT
 * @param NWin Window length
 */
void deleteSTFT(cdouble** S, int NWin);

/**
 * Free the memory associated to a spectrogram
 * @param S Spectrogram
 * @param win Window length
 */
void deleteSpectrogram(double** S, int win);

#endif