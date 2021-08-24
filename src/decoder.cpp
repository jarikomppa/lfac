#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define DR_WAV_IMPLEMENTATION
#include "../external/dr_wav.h"

char outfilename[1024];
unsigned char* rawdata;
int rawlen;
unsigned char* outdata;

void load_data(const char* fn)
{
	FILE* f = fopen(fn, "rb");
	if (!f)
	{
		printf("Can't open source file\n");
		exit(0);
	}
	fseek(f, 0, SEEK_END);
	rawlen = ftell(f);
	fseek(f, 0, SEEK_SET);
	rawdata = new unsigned char[rawlen];
	fread(rawdata, rawlen, 1, f);
	fclose(f);
}

int samplerate, channels, dimensions, window, datasize;

void decode()
{
	int* hd = (int*)rawdata;
	if (hd[0] != 'CAFL' || hd[1] != 0)
	{
		printf("Source data header not recognized. Aborting.\n");
		exit(0);
	}
	
	samplerate = hd[2];
	channels = hd[3];
	dimensions = hd[4];
	window = hd[5];
	datasize = hd[6];

	printf(
		"Samplerate: %d\n"
		"Channels:   %d\n"
		"Dimensions: %d\n"
		"Window:     %d\n"
		"Data size:  %d\n",
		samplerate, channels, dimensions, window, datasize);

	unsigned char* dictionary = rawdata + 8 * sizeof(int);
	unsigned char* data = dictionary + dimensions * 256;

	outdata = new unsigned char[datasize];

	unsigned char* out = outdata;
	int i = 0;
	int ws = window;
	int dataleft = datasize;
	while (dataleft > 0)
	{ 
		for (int j = 0; j < dimensions; j++)
		{
			*out = dictionary[data[i] * dimensions + j];
			out++;
			dataleft--;
		}
		i++;
		ws--;
		if (ws == 0)
		{
			ws = window;
			dictionary = data + window;
			data = dictionary + dimensions * 256;
			i = 0;
		}
	}
}

void save_data()
{
	drwav_data_format format;
	format.container = drwav_container_riff;     // <-- drwav_container_riff = normal WAV files, drwav_container_w64 = Sony Wave64.
	format.format = DR_WAVE_FORMAT_PCM;          // <-- Any of the DR_WAVE_FORMAT_* codes.
	format.channels = channels;
	format.sampleRate = samplerate;
	format.bitsPerSample = 8;
	drwav wav;
	if (!drwav_init_file_write(&wav, outfilename, &format, NULL))
	{
		printf("Unable to open output file\n");
		exit(0);
	}
	drwav_write_pcm_frames(&wav, datasize / channels, outdata);
	drwav_uninit(&wav);
}

int main(int parc, char** pars)
{
	printf("LFAC decoder by Jari Komppa 2021 http://iki.fi/sol\n");
	if (parc < 2)
	{		
		printf("Usage: %s infilename [outfilename]\n", pars[0]);
		return 0;
	}
	if (parc < 3)
	{
		sprintf(outfilename, "%s.wav", pars[1]);
	}
	else
		sprintf(outfilename, "%s", pars[2]);

	printf("Decoding %s to %s\n", pars[1], outfilename);

	load_data(pars[1]);

	decode();
	save_data();
	printf("Done.\n");
	return 0;
}