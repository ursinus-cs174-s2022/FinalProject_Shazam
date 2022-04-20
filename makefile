CC=g++
CFLAGS=-std=c++11 -g -Wall 

all: plotSpectrogram

simplecanvas.o: simplecanvas/simplecanvas.h simplecanvas/simplecanvas.cpp
	$(CC) $(CFLAGS) -c simplecanvas/simplecanvas.cpp

dsp.o: dsp.cpp dsp.h
	$(CC) $(CFLAGS) -c dsp.cpp

dspviz.o: dsp.o dspviz.cpp dspviz.h simplecanvas.o
	$(CC) $(CFLAGS) -c dspviz.cpp

plotSpectrogram: dsp.o dspviz.o simplecanvas.o AudioFile.h plotSpectrogram.cpp
	$(CC) $(CFLAGS) -o plotSpectrogram dsp.o dspviz.o simplecanvas.o plotSpectrogram.cpp

clean:
	rm *.o *.exe *.stackdump plotSpectrogram