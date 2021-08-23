#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define DR_WAV_IMPLEMENTATION
#include "../external/dr_wav.h"
#define DR_MP3_IMPLEMENTATION
#include "../external/dr_mp3.h"
#define DR_FLAC_IMPLEMENTATION
#include "../external/dr_flac.h"
#include "../external/stb_vorbis.c"
#include <iostream>
#include "../external/optionparser.h"
#include "../external/libsamplerate/samplerate.h"

#define MAX_DIMENSIONS 64
#define MAX_GROUPS 256
unsigned char* gChunkdata;
unsigned char* outdata;
unsigned int* index;
int group[MAX_GROUPS];
int groupofs[MAX_GROUPS];
int groups = 1;
int analysis[MAX_GROUPS * MAX_DIMENSIONS];
int chunks = 0;
int dimensions = 4;
int datalen = 0;
unsigned char dictionary[MAX_GROUPS * MAX_DIMENSIONS];
unsigned char* unpackeddata;

int window_offset = 0;
int window_size = 0;

// It may seem a bit weird to store the data this way (as a linear
// pass would work just as well), but this way we can do weird things
// with chunks if we want to.
void insert_chunk(unsigned char* p)
{
	memcpy(gChunkdata + chunks * dimensions, p, dimensions);
	chunks++;
}

void analyze_group(int g)
{
	unsigned char* chunkdata = gChunkdata + window_offset;
	int minval[MAX_DIMENSIONS];
	int maxval[MAX_DIMENSIONS];
	for (int i = 0; i < dimensions; i++)
	{
		int v = chunkdata[index[groupofs[g]] * dimensions + i];
		minval[i] = v;
		maxval[i] = v;
	}
	for (int i = 0; i < group[g]; i++)
	{
		for (int j = 0; j < dimensions; j++)
		{
			int v = chunkdata[index[groupofs[g] + i] * dimensions + j];
			if (minval[j] > v) minval[j] = v;
			if (maxval[j] < v) maxval[j] = v;
		}
	}
	for (int i = 0; i < dimensions; i++)
		analysis[g * dimensions + i] = maxval[i] - minval[i];
}

int find_largest()
{
	int v = 0, d = 0, gs = group[0];
	for (int i = 0; i < dimensions * groups; i++)
	{
		if (analysis[i] > v)
		{
			d = i;
			v = analysis[i];
			gs = group[i / dimensions];
		} else
		if (analysis[i] == v && gs < group[i/dimensions])
		{
			d = i;
			v = analysis[i];
			gs = group[i / dimensions];
		}

	}
	return d;
}

int sortdimension = 0;

int cmpfunc(const void* a, const void* b)
{
	unsigned char* chunkdata = gChunkdata + window_offset;
	int av = *(int*)a;
	int bv = *(int*)b;
	return (int)chunkdata[av * dimensions + sortdimension] - 
		   (int)chunkdata[bv * dimensions + sortdimension];
}

void sort_group(int g, int d)
{
	sortdimension = d;
	qsort(index + groupofs[g], group[g], sizeof(int), cmpfunc);
}

FILE* outfile;
unsigned int channels;
unsigned int samplerate;
unsigned int window;
float* sampledata = NULL;

void load_data(const char* fn)
{
	// Let's assume wav.
	drwav_uint64 totalPCMFrameCount;
	sampledata = drwav_open_file_and_read_pcm_frames_f32(fn, &channels, &samplerate, &totalPCMFrameCount, NULL);
	datalen = ((int)totalPCMFrameCount / dimensions) * dimensions; // rounded to "dimensions"

	if (!sampledata)
	{
		// not a wav, what about mp3?
		drmp3_config conf;
		sampledata = drmp3_open_file_and_read_pcm_frames_f32(fn, &conf, &totalPCMFrameCount, NULL);
		if (sampledata)
		{
			samplerate = conf.sampleRate;
			channels = conf.channels;
		}
	}

	if (!sampledata)
	{
		// flac, maybe?
		sampledata = drflac_open_file_and_read_pcm_frames_f32(fn, &channels, &samplerate, &totalPCMFrameCount, NULL);
	}

	if (!sampledata)
	{
		// it must be an ogg then.
		short* output;
		int frames = stb_vorbis_decode_filename(fn, (int*)&channels, (int*)&samplerate, &output);	
		if (frames > 0)
		{
			datalen = frames * channels;
			sampledata = new float[datalen];
			for (int i = 0; i < datalen; i++)
				sampledata[i] = output[i] * (1.0f / 0x7fff);
			free(output);
		}
	}

	if (!sampledata)
	{
		// okay, keep your secrets.
		printf("Failed to load data\n");
		exit(0);
	}
}

void resample(int targetsamplerate, int mono, int fast)
{
	if (mono && channels != 1)
	{
		printf("Mixing %d channels to mono..\n", channels);
		// Mix down to mono by adding all channels together.
		float* t = new float[datalen / channels];
		for (unsigned int i = 0; i < (datalen / channels); i++)
		{
			t[i] = 0;
			for (unsigned int j = 0; j < channels; j++)
				t[i] += sampledata[i * channels + j];
		}
		delete[] sampledata;
		sampledata = t;
		datalen /= channels;
		channels = 1;
	}

	if (targetsamplerate != samplerate)
	{
		printf("Resampling from %d to %d (libsamplerate %s)..\n", samplerate, targetsamplerate, fast?"SRC_SINC_FASTEST":"SINC_BEST_QUALITY");
		SRC_DATA d;
		d.data_in = sampledata;
		d.input_frames = datalen / channels;
		d.src_ratio = targetsamplerate / (float)samplerate;
		d.output_frames = (int)((datalen / channels) * d.src_ratio);
		d.data_out = new float[d.output_frames];
		
		if (src_simple(&d, fast?SRC_SINC_FASTEST:SRC_SINC_BEST_QUALITY, channels))
		{
			printf("Resampling failed.\n");
			exit(0);
		}

		delete[] sampledata;
		sampledata = d.data_out;
		datalen = d.output_frames_gen * channels;
		samplerate = targetsamplerate;
	}

	gChunkdata = new unsigned char[datalen];
	index = new unsigned int[datalen / dimensions];
	unsigned char tp[MAX_DIMENSIONS];
	for (int i = 0; i < datalen / dimensions; i++)
	{
		for (int j = 0; j < dimensions; j++)
			tp[j] = (unsigned char)((sampledata[dimensions * i * channels + j] + 1.0) * 127);
		insert_chunk(tp);
	}
}

void prep_output(const char* fn)
{
	outfile = fopen(fn, "wb");
	int tag = 'CAFL'; // gets reversed
	fwrite(&tag, 1, 4, outfile);
	int version = 0;
	fwrite(&version, 1, 4, outfile);
	fwrite(&samplerate, 1, 4, outfile);
	fwrite(&channels, 1, 4, outfile);
	fwrite(&dimensions, 1, 4, outfile);
	fwrite(&window, 1, 4, outfile);
	fwrite(&datalen, 1, 4, outfile);
	int reserved = 0;
	fwrite(&reserved, 1, 4, outfile);
}

void init_encode()
{
	for (int i = 0; i < MAX_GROUPS; i++)
	{
		group[i] = 0;
		groupofs[i] = 0;
	}

	groups = 1;

	chunks = window_size / dimensions;

	for (int i = 0; i < chunks; i++)
		index[i] = i;

	for (int i = 0; i < MAX_GROUPS * MAX_DIMENSIONS; i++)
		analysis[i] = -1;

	group[0] = chunks; // first group contains all chunks
}

enum {
	CUT_MEDIAN,
	CUT_AVERAGE,
	CUT_MEAN,
	CUT_EVEN
};

void split_group(int g, int d, int cut_type)
{
	unsigned char* chunkdata = gChunkdata + window_offset;
	//printf("Splitting group %3d, dimension %d, delta %3d (%d items)\n", g, d, analysis[g*dimensions+d], group[g]);
	sort_group(g, d);
	int total = 0;
	int minval = 255, maxval = 0;
	for (int i = 0; i < group[g]; i++)
	{
		total += chunkdata[index[groupofs[g] + i] * dimensions + d];		
		if (chunkdata[index[groupofs[g] + i] * dimensions + d] > maxval) maxval = chunkdata[index[groupofs[g] + i] * dimensions + d];
		if (chunkdata[index[groupofs[g] + i] * dimensions + d] < minval) minval = chunkdata[index[groupofs[g] + i] * dimensions + d];
	}
	if (total < 0)
	{
		printf("borked\n");
		exit(0);
	}
	int i = 0;
	if (cut_type == CUT_MEAN)
	{
		// average cut
		int midval = minval + (maxval - minval) / 2;
		while (chunkdata[index[groupofs[g] + i] * dimensions + d] < midval) i++;
	}
	if (cut_type == CUT_AVERAGE)
	{
		int split = total / 2;
		total = 0;
		// median cut
		while (total < split)
		{
			total += chunkdata[index[groupofs[g] + i] * dimensions + d];
			i++;
		}
	}
	if (cut_type == CUT_EVEN)
	{
		i = group[g] / 2; // halve space
	}
	if (cut_type == CUT_MEDIAN)
	{
		i = group[g] / 2;
		int v = chunkdata[index[groupofs[g] + i] * dimensions + d];
		while (i > 0 && v == chunkdata[index[groupofs[g] + i] * dimensions + d]) i--;
		if (i == 0)
		{
			i = group[g] / 2;
			v = chunkdata[index[groupofs[g] + i] * dimensions + d];
			while (i < group[g] && v == chunkdata[index[groupofs[g] + i] * dimensions + d]) i++;
			if (i == group[g])
				i /= 2;
		}
	}

	if (i == group[g])
		i--; // avoid splitting to empty groups
	if (i < 1) return;
	group[groups] = group[g] - i;
	groupofs[groups] = groupofs[g] + i;
	group[g] = i;
	analyze_group(g);
	analyze_group(groups);
	groups++;
}

float dist(const unsigned char* a, const unsigned char* b)
{
	float s = 0;
	for (int i = 0; i < dimensions; i++)
		s += ((float)a[i] - (float)b[i]) * ((float)a[i] - (float)b[i]);
	return (float)sqrt(s);
}

void average_groups()
{
	unsigned char* chunkdata = gChunkdata + window_offset;
	for (int i = 0; i < groups * dimensions; i++)
		dictionary[i] = 0;

	// Calculate average chunk
	for (int g = 0; g < groups; g++)
	{
		for (int d = 0; d < dimensions; d++)
		{
			int total = 0;
			for (int i = 0; i < group[g]; i++)
				total += chunkdata[index[groupofs[g] + i] * dimensions + d];
			dictionary[g * dimensions + d] = total / group[g];
		}
	}
	
}

float verify()
{
	unsigned char* chunkdata = gChunkdata + window_offset;
	// decompress the data and calculate error compared to source.
	for (int i = 0; i < chunks; i++)
	{
		for (int j = 0; j < dimensions; j++)
		{
			unpackeddata[i * dimensions + j + window_offset] = dictionary[outdata[i] * dimensions + j];
		}
	}
	long long errsum = 0;
	for (int i = 0; i < chunks * dimensions; i++)
	{
		int d = abs(unpackeddata[i + window_offset] - chunkdata[i]);
		errsum += d;
	}

	return errsum / (double)(chunks * dimensions);
}


void map_indices(int maxiter)
{
	unsigned char* chunkdata = gChunkdata + window_offset;
	outdata = new unsigned char[chunks];
	int changed = 1;
	int timeout = 0;
	float err = 1000000;
	float preverr = err;
	while (changed && timeout < maxiter && err <= preverr)
	{
		preverr = err;
		printf("%c\r", "\\-/|"[timeout % 4]);
		timeout++;
		changed = 0;
		for (int i = 0; i < chunks; i++)
		{
			int idx = 0;
			float distance = dist(dictionary, chunkdata + i * dimensions);

			for (int g = 1; g < groups; g++)
			{
				float d = dist(dictionary + g * dimensions, chunkdata + i * dimensions);
				if (d < distance)
				{
					distance = d;
					idx = g;
				}
			}
			if (outdata[i] != idx)
			{
				outdata[i] = idx;
				changed = 1;
			}
		}

		// Recalculate averages
		for (int d = 0; d < dimensions; d++)
		{
			for (int g = 0; g < 256; g++)
			{
				int total = 0;
				int count = 0;
				for (int i = 0; i < chunks; i++)
				{
					if (outdata[i] == g)
					{
						total += *(chunkdata + i * dimensions + d);
						count++;
					}
				}
				if (count != 0)
					dictionary[g * dimensions + d] = total / count;
			}
		}
		err = verify();
	}
	printf("%d iterations, avg error %3.3f\n", timeout, err);
	fwrite(dictionary, dimensions * 256, 1, outfile);
	fwrite(outdata, chunks, 1, outfile);
}

void finish()
{
	fclose(outfile);
}

void reduce(int cut_type, int part, int parts)
{
	printf("Compressing part %d/%d with %d dimensions at %dHz, window size %d..\n", part, parts, dimensions, samplerate, window_size);
	analyze_group(0);
//	while (groups < MAX_GROUPS)
	for (int i = 1; i < MAX_GROUPS; i++)
	{
		int t = find_largest();
		int g = t / dimensions;
		int d = t % dimensions;
		split_group(g, d, cut_type);
	}
}


void save_data(const char* fn, unsigned char*data)
{
	drwav_data_format format;
	format.container = drwav_container_riff;
	format.format = DR_WAVE_FORMAT_PCM;
	format.channels = channels;
	format.sampleRate = samplerate;
	format.bitsPerSample = 8;
	drwav wav;
	if (!drwav_init_file_write(&wav, fn, &format, NULL))
	{
		printf("Failed to open \"%s\" for writing\n", fn);
		return;
	}
	drwav_write_pcm_frames(&wav, chunks * dimensions / channels, data);
	drwav_uninit(&wav);
}

void save_compare_data(const char* fn, int fast)
{
	int windowing_extra_data = window ? (datalen / (dimensions * window)) * 256 * dimensions : 0;
	int compress_size = datalen / dimensions + 256 * dimensions + windowing_extra_data;
	float compressionratio = compress_size / (float)datalen;
	int targetsamplerate = (int)(samplerate * compressionratio);
	
	printf("Saving compare file \"%s\" at sample rate %dHz\n", fn, targetsamplerate);

	SRC_DATA d;
	d.data_in = sampledata;
	d.input_frames = datalen / channels;
	d.src_ratio = targetsamplerate / (float)samplerate;
	d.output_frames = (int)((datalen / channels) * d.src_ratio);
	d.data_out = new float[d.output_frames];

	if (src_simple(&d, fast ? SRC_SINC_FASTEST : SRC_SINC_BEST_QUALITY, channels))
	{
		printf("Resampling failed.\n");
		exit(0);
	}

	unsigned char* data = new unsigned char[datalen];
	for (int i = 0; i < datalen; i++)
		data[i] = (unsigned char)((d.data_out[i] * 0.5 + 0.5) * 255);

	drwav_data_format format;
	format.container = drwav_container_riff;
	format.format = DR_WAVE_FORMAT_PCM;
	format.channels = channels;
	format.sampleRate = targetsamplerate;
	format.bitsPerSample = 8;
	drwav wav;
	if (!drwav_init_file_write(&wav, fn, &format, NULL))
	{
		printf("Failed to open \"%s\" for writing\n", fn);
		return;
	}
	drwav_write_pcm_frames(&wav, d.output_frames_gen * channels, data);
	drwav_uninit(&wav);

	delete[] data;
	delete[] d.data_out;
}

enum optionIndex { UNKNOWN, HELP, SAMPLERATE, DIMENSIONS, MONO, WINDOW, SAVE, SAVESRC, SAVECMP, FASTRS, CUTTYPE, MAXITER };
const option::Descriptor usage[] =
{
	{ UNKNOWN,		0, "", "",	option::Arg::None,				 "USAGE: encoder inputfilename outputfilename [options]\n\nOptions:"},
	{ HELP,			0, "h", "help", option::Arg::None,			 " -h --help\t Print usage and exit"},
	{ SAMPLERATE,	0, "r", "samplerate", option::Arg::Optional, " -r --samplerate=sr\t Set target samplerate (default: use source)"},
	{ DIMENSIONS,	0, "d", "dimensions", option::Arg::Optional, " -d --dimensions=dim\t Set number of dimensions (default 4)"},
	{ MONO,			0, "m", "mono", option::Arg::None,			 " -m --mono\t Mix to mono (default: use source)"},
	{ WINDOW,		0, "w", "window", option::Arg::Optional,	 " -w --window=winsize\t Set window size in grains (default: infinite)"},
	{ SAVE,         0, "o", "saveout", option::Arg::Optional,    " -o --saveout=filename\t Save re-decompressed file (default: don't)"},
	{ SAVESRC,      0, "s", "savesrc", option::Arg::Optional,    " -s --savesrc=filename\t Save raw soure data (after resampling) (default: don't)"},
	{ SAVECMP,      0, "c", "savecmp", option::Arg::Optional,    " -c --savecmp=filename\t Save size-comparable low-freq wav (default: don't)"},
	{ FASTRS,       0, "f", "fastresample", option::Arg::None,   " -f --fastresample\t Use fast resampler (default: SINC_BEST)"},
	{ CUTTYPE,		0, "x", "cuttype", option::Arg::Optional,    " -x --cuttype=type\t Subspace cut type: even, mean, average, median. (default: median)"},
	{ MAXITER,      0, "i", "maxiter", option::Arg::Optional,    " -i --maxiter=iters\t Maximum iterations for re-centering grains (default:10)"},
	{ UNKNOWN,      0, "", "", option::Arg::None,				 "Example:\n  encoder dasboot.mp3 theshoe.sad -m --window=65536 -d 16"},
	{ 0,0,0,0,0,0 }
};

int main(int parc, char** pars)
{
	printf("LFAC encoder by Jari Komppa 2021 http://iki.fi/sol\n");

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
		dimensions = atoi(options[DIMENSIONS].arg);
		if (dimensions < 2 || dimensions > 64)
		{
			printf("Bad value for dimensions. Try something like 2 or 16.\n");
			return 0;
		}
		if (dimensions > 16)
		{
			printf("Note: dimensions set quite high. Quality will likely be horrible.\n");
		}
	}

	if (options[WINDOW] && options[WINDOW].arg)
	{
		window = atoi(options[WINDOW].arg);
		if (window == 0)
		{
			// fine
		}
		else
		if (window < dimensions * dimensions * 256 / (dimensions - 1))
		{
			printf("Given window size (%d) would cause resulting file to be bigger than original. (break even at %d)\n", window, dimensions * dimensions * 256 / (dimensions - 1));
			return 0;
		}
	}

	load_data(parse.nonOption(0));
	if (window > datalen / dimensions)
		window = datalen / dimensions;
	unpackeddata = new unsigned char[datalen];
	int sr = samplerate;
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

	if (options[SAVESRC] && options[SAVESRC].arg && strlen(options[SAVESRC].arg) > 0)
	{
		printf("Saving source data as \"%s\"\n", options[SAVESRC].arg);
		save_data(options[SAVESRC].arg, gChunkdata);
	}

	if (options[SAVECMP] && options[SAVECMP].arg && strlen(options[SAVECMP].arg) > 0)
		save_compare_data(options[SAVECMP].arg, !!options[FASTRS]);

	unsigned int total_left = datalen;
	prep_output(parse.nonOption(1));
	int part = 0;
	int parts = 1;
	if (window)
		parts += datalen / (window * dimensions);
	while (total_left > 0)
	{
		part++;
		printf("total left: %d\n", total_left);
		if (total_left > window * dimensions)
			window_size = window * dimensions;
		else
			window_size = total_left;
		if (window == 0) window_size = datalen;

		init_encode();
		reduce(cut_type, part, parts);
		average_groups();
		map_indices(maxiters);
		verify();
		
		window_offset += window_size;
		if (window != 0 && window * dimensions < total_left)
			total_left -= window * dimensions;
		else
			total_left = 0;
	}
	finish();
	if (options[SAVE] && options[SAVE].arg && strlen(options[SAVE].arg) > 0)
	{
		printf("Saving uncompressed data as \"%s\"\n", options[SAVE].arg);
		save_data(options[SAVE].arg, unpackeddata);
	}
	return 0;
}