#include "dsp.h"

/**
 * Compute the closest power of 2 greater than or equal
 * to some number
 * @param a The number
 * @return The closest power of 2 >= a
 */
int getClosestPowerOf2(int a) {
	double lg = log((double)a) / log(2.0);
	int power = (int)lg;
	if ((double)((int)lg) < lg) {
		power++;
	}
	return power;
}

/**
 * Compute a version of x which is bit-reversed
 * @param x A 32-bit int to reverse
 * @param length Length of bits
 * @return Bit-reversed version of x
 */
int bitReverse(int x, int length) {
	int toReturn = 0;
	int mirror = length / 2;
	for (int mask = 1; mask <= length; mask <<= 1, mirror >>= 1) {
		if ((mask & x) > 0)
			toReturn |= mirror;
	}
	return toReturn;
}


/**
 * Rearrange the terms in-place so that they're sorted by the least
 * significant bit (this is the order in which the terms are accessed
 * in the FFT)
 * @param a An array of complex numbers
 * @param N Number of elements in the array
 */
void rearrange(cdouble* a, int N) {
	for (int i = 0; i < N; i++) {
		int j = bitReverse(i, N);
		if (j > i) { //Don't waste effort swapping two mirrored
		//elements that have already been swapped
			cdouble temp = a[j];
			a[j] = a[i];
			a[i] = temp;
		}
	}
}

/**
 * Initialize the complex coefficients for the FFT
 * @param fftsize The length of the fft
 */
void DSP::initCoeffs(int fftsize) {
	int maxlevel = getClosestPowerOf2(fftsize) + 1;
	W = new cdouble**[maxlevel+1];
	for (int level = 1; level <= maxlevel; level++) {
		int FFTSize = 1 << level;
		W[level] = new cdouble*[2];
		W[level][0] = new cdouble[FFTSize >> 1];
		W[level][1] = new cdouble[FFTSize >> 1];
		for (int i = 0; i < FFTSize >> 1; i++) {
			double iangle = (double)i * 2.0 * PI / (double)FFTSize;
			double fangle = (-1.0) * iangle;
			W[level][FFT_FORWARD][i] = cdouble(cos(fangle), sin(fangle));
			W[level][FFT_INVERSE][i] = cdouble(cos(iangle), sin(iangle)); 
		}
	}
}

/**
 * Initialize the coefficients in Hann window
 * @param fftsize Lenght of fft
 */
void DSP::initWindow(int fftsize) {
	int N = fftsize;
	hannwindow = new double[N];
	for (int n = 0; n < N; n++) {
		double angle = 2.0*PI * n / (double)(N - 1);
		//Do a hann window for now
		hannwindow[n] = 0.54 - 0.46*cos(angle);
	}
}

/**
 * Perform an in-place Cooley-Tukey FFT
 * @param toReturn Array that holds FFT coefficients
 * @param N Length of array (assumed to be power of 2)
 * @param inverse Whether this is a forward or inverse FFT
 */
void DSP::performfft(cdouble* toReturn, int N, int inverse) {
	rearrange(toReturn, N);
	//Do the trivial FFT size of 2 first
	for (int i = 0; i < N; i += 2) {
		cdouble temp = toReturn[i];
		toReturn[i] = temp + toReturn[i + 1];
		toReturn[i + 1] = temp - toReturn[i + 1];
	}
	int Mindex = 2;//Index used to access the cached complex
	//coefficients
	for (int level = 2; level < N; level <<= 1) {
		int FFTSize = level << 1;
		for (int start = 0; start < N; start += FFTSize) {
			//This is a little chunk of an FFT of size "FFTSize"
			//to do in-place with the merging algorithm
			//NOTE: "level" gives the length between mirrored terms
			for (int i = 0; i < level; i++) {
				cdouble coeff = W[Mindex][inverse][i];
				cdouble first = toReturn[start + i];
				cdouble second = coeff*toReturn[start + i + level];
				toReturn[start + i] = first + second;
				toReturn[start + i + level] = first - second;
			}
		}
		Mindex++;
	}
}
	

DSP::DSP(int fftsize) {
	this->fftsize = fftsize;
	this->initCoeffs(fftsize);
	this->initWindow(fftsize);
}
DSP::~DSP() {
	int maxlevel = getClosestPowerOf2(fftsize) + 1;
	// Clean up FFT coefficients
	for (int level = 1; level <= maxlevel; level++) {
		for (int type = 0; type < 2; type++) {
			delete[] W[level][type];
		}
		delete[] W[level];
	}
	delete[] W;
	// Clean up window coefficients
	delete[] hannwindow;
}

/**
 * Implement the dft directly from the definition (used for speed comparison)
 * @param sig Complex signal on which to compute dft
 * @param toReturn The array that will hold the fourier coefficients
 * @param N Length of signal
 * @return Complex DFT coefficients
 */
void DSP::dft(cdouble* sig, cdouble* toReturn, int N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double angle = -2.0 * PI * (double)i * (double)j / (double)N;
			cdouble coeff(cos(angle), sin(angle));
			toReturn[i] = coeff*sig[i];
		}
	}
}

/**
 * Perform the FFT on a complex signal
 * @param sig The signal
 * @param toReturn The array that will hold the fourier coefficients
 * @param N Length of the signal (assumed to be power of 2)
 * @return An N-length array with FFT coefficients
 */
void DSP::fft(cdouble* sig, cdouble* toReturn, int N) {
	for (int i = 0; i < N; i++) {
		toReturn[i] = sig[i];
	}
	performfft(toReturn, N, FFT_FORWARD);	
}

/**
 * Perform the inverse FFT on an array of complex FFT coefficients
 * @param sig The FFT coefficients
 * @param toReturn The array that will hold the complex time series
 * @param N Length of the FFT coefficients (assumed to be power of 2)
 * @return An N-length array with FFT coefficients
 */
void DSP::ifft(cdouble* sig, cdouble* toReturn, int N) {
	for (int i = 0; i < N; i++) {
		toReturn[i] = sig[i];
		//Scale by 1/N for inverse FFT
		toReturn[i] *= cdouble(1.0/(double)N, 0);
	}
	performfft(toReturn, N, FFT_INVERSE);
}

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
void DSP::toWindowedComplexArray(double* data, cdouble* res, int N, int start, int win, bool useWindow) {
	//Make a complex array out of the real array
	for (int i = 0; i < win; i++) {
		if (start+i < N) {
			res[i] = cdouble(data[start + i], 0.0);
			if (useWindow) {
				res[i] *= hannwindow[i];
			}
		}
		else {
			//Zero pad if not a power of 2
			res[i] = cdouble(0.0, 0.0);
		}
	}
}

/**
 * Perform a short-time fourier transform on a bunch of samples
 * @param sig Samples in the signal
 * @param N Length of signal
 * @param win Window length (Assumed to be a power of 2)
 * @param hop Hop length
 * @param useWindow Whether to use the window
 * @param NWin Number of windows (returned by reference)
 * @return An NWin x win 2D array of complex doubles
 */
cdouble** DSP::stft(double* sig, int N, int win, int hop, bool useWindow, int* NWin) {
	*NWin = 1 + round((N-win)/(double)hop);
	cdouble** S = new cdouble*[*NWin];
	for (int i = 0; i < *NWin; i++) {
		S[i] = new cdouble[win];
	}
	cdouble* ffti = new cdouble[win];
	for (int i = 0; i < *NWin; i++) {
		toWindowedComplexArray(sig, S[i], N, i*hop, win, useWindow);
		fft(S[i], ffti, win);
		for (int j = 0; j < win; j++) {
			S[i][j] = ffti[j];
		}
	}
	delete[] ffti;
	return S;
}

/**
 * Perform a magnitude short-time fourier transform on a bunch of samples
 * @param sig Samples in the signal
 * @param N Length of signal
 * @param win Window length of STFT (Assumed to be a power of 2)
 * @param hop Hop length
 * @param useWindow Whether to use the window
 * @param NWin Number of windows (returned by reference)
 * @return A win x NWin 2D array of complex doubles
 */
double** DSP::spectrogram(double* sig, int N, int win, int hop, bool useWindow, int* NWin) {
	cdouble** SComplex = stft(sig, N, win, hop, useWindow, NWin);
	double** S = new double*[win];
	for (int i = 0; i < win; i++) {
		S[i] = new double[*NWin];
		for (int j = 0; j < *NWin; j++) {
			S[i][j] = abs(SComplex[j][i]);
		}
	}
	deleteSTFT(SComplex, *NWin);
	return S;
}

/**
 * Free the memory associated to an STFT
 * @param S Spectrogram
 * @param win Window length
 */
void deleteSTFT(cdouble** S, int NWin) {
	for (int i = 0; i < NWin; i++) {
		delete[] S[i];
	}
	delete[] S;
}

/**
 * Free the memory associated to a spectrogram
 * @param S Spectrogram
 * @param win Window length
 */
void deleteSpectrogram(double** S, int win) {
	for (int i = 0; i < win; i++) {
		delete[] S[i];
	}
	delete[] S;
}