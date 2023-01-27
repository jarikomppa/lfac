#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream> // needed by optionparser
#define DR_WAV_IMPLEMENTATION
#include "../external/dr_wav.h"
#define DR_MP3_IMPLEMENTATION
#include "../external/dr_mp3.h"
#define DR_FLAC_IMPLEMENTATION
#include "../external/dr_flac.h"
#include "../external/stb_vorbis.c"
#include "../external/optionparser.h"
#include "../external/libsamplerate/samplerate.h"

#include "vq.hpp"

#define MAX_DIMENSIONS 64
#define MAX_GROUPS 256
unsigned char* gChunkData;
unsigned char* gOutData;
unsigned char* gUnpackedData;
unsigned int* gIndex;
unsigned int gGroup[MAX_GROUPS];
unsigned int gGroupofs[MAX_GROUPS];
unsigned int gGroupCount = 1;
unsigned int gAnalysis[MAX_GROUPS * MAX_DIMENSIONS];
unsigned int gChunks = 0;
unsigned int gDimensions = 4;
unsigned int gDatalen = 0;
unsigned char gDictionary[MAX_GROUPS * MAX_DIMENSIONS];
int gWindowOffset = 0;
int gWindowSize = 0;
FILE* gOutFile;
unsigned int gChannels;
unsigned int gSampleRate;
unsigned int gWindow;
float* gSampleData = NULL;


// It may seem a bit weird to store the data this way (as a linear
// pass would work just as well), but this way we can do weird things
// with chunks if we want to.
void insert_chunk(unsigned char* p)
{
	memcpy(gChunkData + gChunks * gDimensions, p, gDimensions);
	gChunks++;
}

void analyze_group(int g)
{
	unsigned char* chunkdata = gChunkData + gWindowOffset;
	int minval[MAX_DIMENSIONS];
	int maxval[MAX_DIMENSIONS];
	for (unsigned int i = 0; i < gDimensions; i++)
	{
		int v = chunkdata[gIndex[gGroupofs[g]] * gDimensions + i];
		minval[i] = v;
		maxval[i] = v;
	}
	for (unsigned int i = 0; i < gGroup[g]; i++)
	{
		for (unsigned int j = 0; j < gDimensions; j++)
		{
			int v = chunkdata[gIndex[gGroupofs[g] + i] * gDimensions + j];
			if (minval[j] > v) minval[j] = v;
			if (maxval[j] < v) maxval[j] = v;
		}
	}
	for (unsigned int i = 0; i < gDimensions; i++)
		gAnalysis[g * gDimensions + i] = maxval[i] - minval[i];
}

int find_largest()
{
	unsigned int v = 0, d = 0, gs = gGroup[0];
	for (unsigned int i = 0; i < gDimensions * gGroupCount; i++)
	{
		if (gAnalysis[i] > v)
		{
			d = i;
			v = gAnalysis[i];
			gs = gGroup[i / gDimensions];
		} else
		if (gAnalysis[i] == v && gs < gGroup[i/gDimensions])
		{
			d = i;
			v = gAnalysis[i];
			gs = gGroup[i / gDimensions];
		}

	}
	return d;
}

int gSortDimension = 0;

int cmpfunc(const void* a, const void* b)
{
	unsigned char* chunkdata = gChunkData + gWindowOffset;
	int av = *(int*)a;
	int bv = *(int*)b;
	return (int)chunkdata[av * gDimensions + gSortDimension] - 
		   (int)chunkdata[bv * gDimensions + gSortDimension];
}

void sort_group(int g, int d)
{
	gSortDimension = d;
	qsort(gIndex + gGroupofs[g], gGroup[g], sizeof(int), cmpfunc);
}

void load_data(const char* fn)
{
	// Let's assume wav.
	drwav_uint64 totalPCMFrameCount;
	gSampleData = drwav_open_file_and_read_pcm_frames_f32(fn, &gChannels, &gSampleRate, &totalPCMFrameCount, NULL);

	if (!gSampleData)
	{
		// not a wav, what about mp3?
		drmp3_config conf;
		gSampleData = drmp3_open_file_and_read_pcm_frames_f32(fn, &conf, &totalPCMFrameCount, NULL);
		if (gSampleData)
		{
			gSampleRate = conf.sampleRate;
			gChannels = conf.channels;
		}
	}

	if (!gSampleData)
	{
		// flac, maybe?
		gSampleData = drflac_open_file_and_read_pcm_frames_f32(fn, &gChannels, &gSampleRate, &totalPCMFrameCount, NULL);
	}

	if (!gSampleData)
	{
		// it must be an ogg then.
		short* output;
		int frames = stb_vorbis_decode_filename(fn, (int*)&gChannels, (int*)&gSampleRate, &output);	
		if (frames > 0)
		{
			gDatalen = frames * gChannels;
			gSampleData = new float[gDatalen];
			for (unsigned int i = 0; i < gDatalen; i++)
				gSampleData[i] = output[i] * (1.0f / 0x7fff);
			free(output);
			totalPCMFrameCount = frames;
		}
	}

	if (!gSampleData)
	{
		// okay, keep your secrets.
		printf("Failed to load data\n");
		exit(0);
	}

	gDatalen = ((int)(totalPCMFrameCount * gChannels) / gDimensions) * gDimensions; // rounded to "gDimensions"
}

void resample(int targetsamplerate, int mono, int fast)
{
	if (mono && gChannels != 1)
	{
		printf("Mixing %d channels to mono.\n", gChannels);
		// Mix down to mono by adding all gChannels together.
		float* t = new float[gDatalen / gChannels];
		float mag = 0;
		float srcmag = 0;
		for (unsigned int i = 0; i < (gDatalen / gChannels); i++)
		{
			t[i] = 0;
			for (unsigned int j = 0; j < gChannels; j++)
			{
				float src = gSampleData[i * gChannels + j];
				t[i] += src;
				if (abs(src) > srcmag) srcmag = abs(src);
			}
			if (abs(t[i]) > mag) mag = abs(t[i]);
		}
		delete[] gSampleData;
		gSampleData = t;
		gDatalen /= gChannels;
		gChannels = 1;
		
		// re-scale magnitude to original levels
		if (srcmag != 0)
		for (unsigned int i = 0; i < gDatalen; i++)
		{
			gSampleData[i] *= srcmag / mag;
		}
	}

	if (targetsamplerate != gSampleRate)
	{
		printf("Resampling from %d to %d (libsamplerate %s).\n", gSampleRate, targetsamplerate, fast ? "SRC_SINC_FASTEST" : "SINC_BEST_QUALITY");
		SRC_DATA d;
		d.data_in = gSampleData;
		d.input_frames = gDatalen / gChannels;
		d.src_ratio = targetsamplerate / (float)gSampleRate;
		d.output_frames = (int)((gDatalen / gChannels) * d.src_ratio);
		d.data_out = new float[d.output_frames * gChannels];
		
		if (src_simple(&d, fast ? SRC_SINC_FASTEST : SRC_SINC_BEST_QUALITY, gChannels))
		{
			printf("Resampling failed.\n");
			exit(0);
		}

		delete[] gSampleData;
		gSampleData = d.data_out;
		gDatalen = d.output_frames_gen * gChannels;
		gSampleRate = targetsamplerate;
	}

	gChunkData = new unsigned char[gDatalen];
	gIndex = new unsigned int[gDatalen / gDimensions];
	unsigned char tp[MAX_DIMENSIONS];
	for (unsigned int i = 0; i < gDatalen / gDimensions; i++)
	{
		for (unsigned int j = 0; j < gDimensions; j++)
			tp[j] = (unsigned char)((gSampleData[gDimensions * i + j] + 1.0) * 127);
		insert_chunk(tp);
	}
}

void prep_output(const char* fn)
{
	gOutFile = fopen(fn, "wb");
	int tag = 'CAFL'; // gets reversed
	fwrite(&tag, 1, 4, gOutFile);
	int version = 0;
	fwrite(&version, 1, 4, gOutFile);
	fwrite(&gSampleRate, 1, 4, gOutFile);
	fwrite(&gChannels, 1, 4, gOutFile);
	fwrite(&gDimensions, 1, 4, gOutFile);
	fwrite(&gWindow, 1, 4, gOutFile);
	fwrite(&gDatalen, 1, 4, gOutFile);
	int reserved = 0;
	fwrite(&reserved, 1, 4, gOutFile);
}

void init_encode()
{
	for (int i = 0; i < MAX_GROUPS; i++)
	{
		gGroup[i] = 0;
		gGroupofs[i] = 0;
	}

	gGroupCount = 1;

	gChunks = gWindowSize / gDimensions;

	for (unsigned int i = 0; i < gChunks; i++)
		gIndex[i] = i;

	for (unsigned int i = 0; i < MAX_GROUPS * MAX_DIMENSIONS; i++)
		gAnalysis[i] = -1;

	gGroup[0] = gChunks; // first group contains all chunks
}

enum {
	CUT_MEDIAN,
	CUT_AVERAGE,
	CUT_MEAN,
	CUT_EVEN
};

void split_group(int g, int d, int cut_type)
{
	unsigned char* chunkdata = gChunkData + gWindowOffset;
	//printf("Splitting group %3d, dimension %d, delta %3d (%d items)\n", g, d, gAnalysis[g*gDimensions+d], gGroup[g]);
	sort_group(g, d);
	int total = 0;
	int minval = 255, maxval = 0;
	for (unsigned int i = 0; i < gGroup[g]; i++)
	{
		total += chunkdata[gIndex[gGroupofs[g] + i] * gDimensions + d];		
		if (chunkdata[gIndex[gGroupofs[g] + i] * gDimensions + d] > maxval) maxval = chunkdata[gIndex[gGroupofs[g] + i] * gDimensions + d];
		if (chunkdata[gIndex[gGroupofs[g] + i] * gDimensions + d] < minval) minval = chunkdata[gIndex[gGroupofs[g] + i] * gDimensions + d];
	}
	if (total < 0)
	{
		printf("borked\n");
		exit(0);
	}
	unsigned int i = 0;
	if (cut_type == CUT_MEAN)
	{
		// average cut
		int midval = minval + (maxval - minval) / 2;
		while (chunkdata[gIndex[gGroupofs[g] + i] * gDimensions + d] < midval) i++;
	}
	if (cut_type == CUT_AVERAGE)
	{
		int split = total / 2;
		total = 0;
		// median cut
		while (total < split)
		{
			total += chunkdata[gIndex[gGroupofs[g] + i] * gDimensions + d];
			i++;
		}
	}
	if (cut_type == CUT_EVEN)
	{
		i = gGroup[g] / 2; // halve space
	}
	if (cut_type == CUT_MEDIAN)
	{
		i = gGroup[g] / 2;
		int v = chunkdata[gIndex[gGroupofs[g] + i] * gDimensions + d];
		while (i > 0 && v == chunkdata[gIndex[gGroupofs[g] + i] * gDimensions + d]) i--;
		if (i == 0)
		{
			i = gGroup[g] / 2;
			v = chunkdata[gIndex[gGroupofs[g] + i] * gDimensions + d];
			while (i < gGroup[g] && v == chunkdata[gIndex[gGroupofs[g] + i] * gDimensions + d]) i++;
			if (i == gGroup[g])
				i /= 2;
		}
	}

	if (i == gGroup[g])
		i--; // avoid splitting to empty groups
	if (i < 1) return;
	gGroup[gGroupCount] = gGroup[g] - i;
	gGroupofs[gGroupCount] = gGroupofs[g] + i;
	gGroup[g] = i;
	analyze_group(g);
	analyze_group(gGroupCount);
	gGroupCount++;
}

double dist(const unsigned char* a, const unsigned char* b)
{
	double s = 0;
	for (unsigned int i = 0; i < gDimensions; i++)
		s += ((double)a[i] - (double)b[i]) * ((double)a[i] - (double)b[i]);
	return sqrt(s);
}

void average_groups()
{
	unsigned char* chunkdata = gChunkData + gWindowOffset;
	for (unsigned int i = 0; i < gGroupCount * gDimensions; i++)
		gDictionary[i] = 0;

	// Calculate average chunk
	for (unsigned int g = 0; g < gGroupCount; g++)
	{
		for (unsigned int d = 0; d < gDimensions; d++)
		{
			if (gGroup[g])
			{
				unsigned int total = 0;
				for (unsigned int i = 0; i < gGroup[g]; i++)
					total += chunkdata[gIndex[gGroupofs[g] + i] * gDimensions + d];
				gDictionary[g * gDimensions + d] = total / gGroup[g];
			}
		}
	}
}

double verify()
{
	unsigned char* chunkdata = gChunkData + gWindowOffset;
	// decompress the data and calculate error compared to source.
	for (unsigned int i = 0; i < gChunks; i++)
	{
		for (unsigned int j = 0; j < gDimensions; j++)
		{
			gUnpackedData[i * gDimensions + j + gWindowOffset] = gDictionary[gOutData[i] * gDimensions + j];
		}
	}
	long long errsum = 0;
	for (unsigned int i = 0; i < gChunks * gDimensions; i++)
	{
		int d = abs(gUnpackedData[i + gWindowOffset] - chunkdata[i]);
		errsum += d;
	}

	return errsum / (double)(gChunks * gDimensions);
}


void map_indices()
{
	unsigned char* chunkdata = gChunkData + gWindowOffset;
	gOutData = new unsigned char[gChunks];
	for (unsigned int i = 0; i < gChunks; i++)
	{
		int idx = 0;
		double distance = dist(gDictionary, chunkdata + i * gDimensions);

		for (unsigned int g = 1; g < gGroupCount; g++)
		{
			double d = dist(gDictionary + g * gDimensions, chunkdata + i * gDimensions);
			if (d < distance)
			{
				distance = d;
				idx = g;
			}
		}
		gOutData[i] = idx;
	}
	double err = verify();
	printf("avg error %3.3f\n", err);
	fwrite(gDictionary, gDimensions * 256, 1, gOutFile);
	fwrite(gOutData, gChunks, 1, gOutFile);
}

void map_indices(int maxiter)
{
	unsigned char* chunkdata = gChunkData + gWindowOffset;
	gOutData = new unsigned char[gChunks];
	int changed = 1;
	int timeout = 0;
	double err = 1000000;
	double preverr = err;
	while (changed && timeout < maxiter && err <= preverr)
	{
		preverr = err;
		printf("%c\r", "\\-/|"[timeout % 4]);
		timeout++;
		changed = 0;
		for (unsigned int i = 0; i < gChunks; i++)
		{
			int idx = 0;
			double distance = dist(gDictionary, chunkdata + i * gDimensions);

			for (unsigned int g = 1; g < gGroupCount; g++)
			{
				double d = dist(gDictionary + g * gDimensions, chunkdata + i * gDimensions);
				if (d < distance)
				{
					distance = d;
					idx = g;
				}
			}
			if (gOutData[i] != idx)
			{
				gOutData[i] = idx;
				changed = 1;
			}
		}

		// Recalculate averages

		for (unsigned int g = 0; g < 256; g++)
		{
			for (unsigned int d = 0; d < gDimensions; d++)
			{
				int total = 0;
				int count = 0;
				for (unsigned int i = 0; i < gChunks; i++)
				{
					if (gOutData[i] == g)
					{
						total += *(chunkdata + i * gDimensions + d);
						count++;
					}
				}
				if (count != 0)
					gDictionary[g * gDimensions + d] = total / count;
				else
					gDictionary[g * gDimensions + d] = 0;
			}
		}

		//average_groups();
		err = verify();
	}
	printf("%d iterations, avg error %3.3f\n", timeout, err);
	fwrite(gDictionary, gDimensions * 256, 1, gOutFile);
	fwrite(gOutData, gChunks, 1, gOutFile);
}

void finish()
{
	fclose(gOutFile);
}

void reduce(int cut_type, int part, int parts)
{
	printf("Compressing part %d/%d with %d dimensions at %dHz, %d byte / %d grain window.\n", part, parts, gDimensions, gSampleRate, gWindowSize, gWindowSize / gDimensions);
	analyze_group(0);
	for (int i = 1; i < MAX_GROUPS; i++)
	{
		int t = find_largest();
		int g = t / gDimensions;
		int d = t % gDimensions;
		split_group(g, d, cut_type);
	}
}


void save_data(const char* fn, unsigned char*data)
{
	drwav_data_format format;
	format.container = drwav_container_riff;
	format.format = DR_WAVE_FORMAT_PCM;
	format.channels = gChannels;
	format.sampleRate = gSampleRate;
	format.bitsPerSample = 8;
	drwav wav;
	if (!drwav_init_file_write(&wav, fn, &format, NULL))
	{
		printf("Failed to open \"%s\" for writing\n", fn);
		return;
	}
	drwav_write_pcm_frames(&wav, gDatalen / gChannels, data);
	drwav_uninit(&wav);
}

void save_compare_data(const char* fn, int fast)
{
	int windowing_extra_data = gWindow ? (gDatalen / (gDimensions * gWindow)) * 256 * gDimensions : 0;
	int compress_size = gDatalen / gDimensions + 256 * gDimensions + windowing_extra_data;
	float compressionratio = compress_size / (float)gDatalen;
	int targetsamplerate = (int)(gSampleRate * compressionratio);
	
	printf("Saving compare file \"%s\" at sample rate %dHz\n", fn, targetsamplerate);

	SRC_DATA d;
	d.data_in = gSampleData;
	d.input_frames = gDatalen / gChannels;
	d.src_ratio = targetsamplerate / (float)gSampleRate;
	d.output_frames = (int)((gDatalen / gChannels) * d.src_ratio);
	d.data_out = new float[d.output_frames * gChannels];

	if (src_simple(&d, fast ? SRC_SINC_FASTEST : SRC_SINC_BEST_QUALITY, gChannels))
	{
		printf("Resampling failed.\n");
		exit(0);
	}

	unsigned char* data = new unsigned char[d.output_frames_gen * gChannels];
	for (unsigned int i = 0; i < d.output_frames_gen * gChannels; i++)
		data[i] = (unsigned char)((d.data_out[i] * 0.5 + 0.5) * 255);

	drwav_data_format format;
	format.container = drwav_container_riff;
	format.format = DR_WAVE_FORMAT_PCM;
	format.channels = gChannels;
	format.sampleRate = targetsamplerate;
	format.bitsPerSample = 8;
	drwav wav;
	if (!drwav_init_file_write(&wav, fn, &format, NULL))
	{
		printf("Failed to open \"%s\" for writing\n", fn);
		return;
	}
	drwav_write_pcm_frames(&wav, d.output_frames_gen, data);
	drwav_uninit(&wav);

	delete[] data;
	delete[] d.data_out;
}

enum optionIndex { UNKNOWN, HELP, SAMPLERATE, DIMENSIONS, MONO, WINDOW, SAVE, SAVESRC, SAVECMP, FASTRS, CUTTYPE, MAXITER, CUTS, OLDENC };
const option::Descriptor usage[] =
{
	{ UNKNOWN,		0, "", "",	option::Arg::None,				 "USAGE: encoder inputfilename outputfilename [options]\n\nOptions:"},
	{ HELP,			0, "h", "help", option::Arg::None,			 " -h --help\t Print usage and exit"},
	{ SAMPLERATE,	0, "r", "samplerate", option::Arg::Optional, " -r --samplerate=sr\t Set target samplerate (default: use source)"},
	{ DIMENSIONS,	0, "d", "dimensions", option::Arg::Optional, " -d --dimensions=dim\t Set number of dimensions (default 4)"},
	{ MONO,			0, "m", "mono", option::Arg::None,			 " -m --mono\t Mix to mono (default: use source)"},
	{ WINDOW,		0, "w", "window", option::Arg::Optional,	 " -w --window=winsize\t Set window size in grains (default: infinite)"},
	{ CUTS,         0, "z", "cuts", option::Arg::Optional,       " -z --cuts=n\t Cut source in N windows (default: one)"},
	{ SAVE,         0, "o", "saveout", option::Arg::Optional,    " -o --saveout=filename\t Save re-decompressed file (default: don't)"},
	{ SAVESRC,      0, "s", "savesrc", option::Arg::Optional,    " -s --savesrc=filename\t Save raw soure data (after resampling) (default: don't)"},
	{ SAVECMP,      0, "c", "savecmp", option::Arg::Optional,    " -c --savecmp=filename\t Save size-comparable low-freq wav (default: don't)"},
	{ FASTRS,       0, "f", "fastresample", option::Arg::None,   " -f --fastresample\t Use fast resampler (default: SINC_BEST)"},
	{ CUTTYPE,		0, "x", "cuttype", option::Arg::Optional,    " -x --cuttype=type\t Subspace cut type: even, mean, average, median. (default: median)"},
	{ MAXITER,      0, "i", "maxiter", option::Arg::Optional,    " -i --maxiter=iters\t Maximum iterations for re-centering grains (default:10)"},
	{ OLDENC,       0, "l", "oldenc", option::Arg::None,         " -l --oldenc\t Use old encoder (default: use new)"},
	{ UNKNOWN,      0, "", "", option::Arg::None,				 "Example:\n  encoder dasboot.mp3 theshoe.sad -m --gWindow=65536 -d16"},
	{ 0,0,0,0,0,0 }
};

int main(int parc, char** pars)
{
	printf("LFAC encoder by Jari Komppa 2021-2023 http://iki.fi/sol\n");

	option::Stats stats(usage, parc - 1, pars + 1);
	assert(stats.buffer_max < 16 && stats.options_max < 16);
	option::Option options[16], buffer[16];
	option::Parser parse(true, usage, parc - 1, pars + 1, options, buffer);

	if (options[UNKNOWN])
	{
		for (option::Option* opt = options[UNKNOWN]; opt; opt = opt->next())
			printf("Unknown option: %s\n", opt->name);
		printf("Run without parameters for help.\n");
		exit(0);
	}

	if (parse.error() || parc < 3 || options[HELP] || parse.nonOptionsCount() != 2)
	{
		option::printUsage(std::cout, usage);
		return 0;
	}

	if (options[DIMENSIONS] && options[DIMENSIONS].arg)
	{
		gDimensions = atoi(options[DIMENSIONS].arg);
		if (gDimensions < 2 || gDimensions > 64)
		{
			printf("Bad value for dimensions. Try something like 2 or 16.\n");
			return 0;
		}
		if (gDimensions > 16)
		{
			printf("Note: dimensions set quite high. Quality will likely be horrible.\n");
		}
	}

	if (options[WINDOW] && options[WINDOW].arg)
	{
		gWindow = atoi(options[WINDOW].arg);
		if (gWindow == 0)
		{
			// fine
		}
		else
		if (gWindow < gDimensions * gDimensions * 256 / (gDimensions - 1))
		{
			printf("Given window size (%d) would cause resulting file to be bigger than original. (break even at %d)\n", gWindow, gDimensions * gDimensions * 256 / (gDimensions - 1));
			return 0;
		}
	}

	load_data(parse.nonOption(0));

	int sr = gSampleRate;
	if (options[SAMPLERATE] && options[SAMPLERATE].arg)
		sr = atoi(options[SAMPLERATE].arg);
	if (sr < 1 || sr > 256000)
	{
		printf("Bad value for sample rate.\n");
		return 0;
	}	

	int cut_type = CUT_MEDIAN;
	if (options[CUTTYPE] && options[CUTTYPE].arg)
	{
		if (_stricmp(options[CUTTYPE].arg, "even")) cut_type = CUT_EVEN; else
		if (_stricmp(options[CUTTYPE].arg, "average")) cut_type = CUT_AVERAGE; else
		if (_stricmp(options[CUTTYPE].arg, "median")) cut_type = CUT_MEDIAN; else
		if (_stricmp(options[CUTTYPE].arg, "mean")) cut_type = CUT_MEAN; else
		{
			printf("Unknown cut type: %s\n", options[CUTTYPE].arg);
			exit(0);
		}
	}

	int maxiters = 10;
	if (options[MAXITER] && options[MAXITER].arg)
		maxiters = atoi(options[MAXITER].arg);
	if (maxiters < 1 || maxiters > 10000)
	{
		printf("Invalid number of max iterations specified\n");
		exit(0);
	}

	resample(sr, !!options[MONO], !!options[FASTRS]);

	if (options[CUTS] && options[CUTS].arg)
	{
		int cuts = atoi(options[CUTS].arg);
		if (cuts <= 0)
		{
			printf("Invalid number of cuts (%d) specified\n", cuts);
			return 0;
		}
		int w = (gDatalen / gDimensions) / cuts;
		// make sure last window doesn't end up being really tiny..
		if (w * cuts * gDimensions < gDatalen) w++;
		gWindow = w;
		printf("Cutting source to %d pieces, yielding %d grain window.\n", cuts, w);
	}


	if (gWindow > gDatalen / gDimensions)
		gWindow = gDatalen / gDimensions;

	gUnpackedData = new unsigned char[gDatalen];


	if (options[SAVESRC] && options[SAVESRC].arg && strlen(options[SAVESRC].arg) > 0)
	{
		printf("Saving source data as \"%s\"\n", options[SAVESRC].arg);
		save_data(options[SAVESRC].arg, gChunkData);
	}

	if (options[SAVECMP] && options[SAVECMP].arg && strlen(options[SAVECMP].arg) > 0)
		save_compare_data(options[SAVECMP].arg, !!options[FASTRS]);

	bool oldenc = options[OLDENC];

	unsigned int total_left = gDatalen;
	prep_output(parse.nonOption(1));
	int part = 0;
	int parts = 1;
	if (gWindow)
		parts += gDatalen / (gWindow * gDimensions);
	while (total_left > 0)
	{
		part++;
		printf("Total bytes left: %d\n", total_left);
		if (total_left > gWindow * gDimensions)
			gWindowSize = gWindow * gDimensions;
		else
			gWindowSize = total_left;
		if (gWindow == 0) gWindowSize = gDatalen;
		if (oldenc)
		{
			// old encoder
			init_encode();
			reduce(cut_type, part, parts);
			average_groups();
			map_indices(maxiters);
			verify();
		}
		else
		{
			init_encode();
			reduce(cut_type, part, parts);
			gGroupCount = 256;
			average_groups();
			// new encoder
			vq::reduce<unsigned char, 256, 0, 256>(
				gChunkData + gWindowOffset,
				gDimensions, 
				gChunks,
				gDictionary);
			map_indices();
			verify();
		}
		gWindowOffset += gWindowSize;
		if (gWindow != 0 && gWindow * gDimensions < total_left)
			total_left -= gWindow * gDimensions;
		else
			total_left = 0;
	}
	finish();
	if (options[SAVE] && options[SAVE].arg && strlen(options[SAVE].arg) > 0)
	{
		printf("Saving uncompressed data as \"%s\"\n", options[SAVE].arg);
		save_data(options[SAVE].arg, gUnpackedData);
	}
	return 0;
}