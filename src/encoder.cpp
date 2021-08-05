#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define DR_WAV_IMPLEMENTATION
#include "../external/dr_wav.h"

#define MAX_DIMENSIONS 16
#define MAX_GROUPS 256
unsigned char* chunkdata;
unsigned char* outdata;
unsigned int* index;
int group[MAX_GROUPS];
int groupofs[MAX_GROUPS];
int groups = 1;
int analysis[MAX_GROUPS * MAX_DIMENSIONS];
int chunks = 0;
int dimensions = 16;
int datalen = 1024 * 8;
unsigned char dictionary[MAX_GROUPS * MAX_DIMENSIONS];
unsigned char* unpackeddata;
char outfilename[1024];

//#define ROTATE_CHUNKS // "it's a great idea. But it doesn't work."
//#define USE_EXISTING_CHUNK // hard to tell

void insert_chunk(unsigned char* p)
{
#ifndef ROTATE_CHUNKS	
	memcpy(chunkdata + chunks * dimensions, p, dimensions);
#else
	int minval = *p;
	int minidx = 0;
	for (int i = 1; i < dimensions; i++)
		if (minval > p[i])
		{
			minval = p[i];
			minidx = i;
		}

	for (int i = 0; i < dimensions; i++)
		chunkdata[chunks * dimensions + i] = p[(minidx + i) % dimensions];
#endif

	index[chunks] = chunks;
	chunks++;
}

void analyze_group(int g)
{
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
	/*
	printf("g %d (%d): ", g, group[g]);
	for (int i = 0; i < dimensions; i++)
		printf("%d ", analysis[g * dimensions + i]);
	printf("\n");
	*/
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

void load_data(const char * fn)
{
	drwav_uint64 totalPCMFrameCount;
	float* pSampleData = drwav_open_file_and_read_pcm_frames_f32(fn, &channels, &samplerate, &totalPCMFrameCount, NULL);
	datalen = ((int)totalPCMFrameCount / dimensions) * dimensions; // rounded to "dimensions"

	if (!pSampleData)
	{
		printf("Failed to load data\n");
		exit(0);
	}
	outfile = fopen(outfilename, "wb");
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

	for (int i = 0; i < MAX_GROUPS; i++)
	{
		group[i] = 0;
		groupofs[i] = 0;
	}
	for (int i = 0; i < MAX_GROUPS * MAX_DIMENSIONS; i++)
		analysis[i] = -1;
	group[0] = datalen / dimensions; // first group contains all chunks
	chunkdata = new unsigned char[datalen];
	index = new unsigned int[datalen / dimensions];
	unsigned char tp[MAX_DIMENSIONS];
	for (int i = 0; i < datalen / dimensions; i++)
	{
		for (int j = 0; j < dimensions; j++)
			tp[j] = (unsigned char)((pSampleData[dimensions * i * channels + j] + 1.0) * 127); //rand();
		insert_chunk(tp);
	}
	//printf("%d chunks (%d max)\n", chunks, datalen / dimensions);
}

void split_group(int g, int d)
{
	//printf("Splitting group %3d, dimension %d, delta %3d (%d items)\n", g, d, analysis[g*dimensions+d], group[g]);
	sort_group(g, d);
	int total = 0;
	for (int i = 0; i < group[g]; i++)
	{
		total += chunkdata[index[groupofs[g] + i] * dimensions + d];		
	}
	if (total < 0)
	{
		printf("borked\n");
		exit(0);
	}
	int split = total / 2;
	total = 0;
	int i = 0;
	while (total < split)
	{
		total += chunkdata[index[groupofs[g] + i] * dimensions + d];
		i++;
	}
	if (i == group[g])
		i--; // avoid splitting to empty groups
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
	for (int i = 0; i < groups * dimensions; i++)
		dictionary[i] = 0;

	// Calculate average chunk
	for (int g = 0; g < groups; g++)
	{
		//printf("(");
		for (int d = 0; d < dimensions; d++)
		{
			int total = 0;
			for (int i = 0; i < group[g]; i++)
				total += chunkdata[index[groupofs[g] + i] * dimensions + d];
			dictionary[g * dimensions + d] = total / group[g];
			//printf("%d ", dictionary[g * dimensions + d]);
		}
		//printf(") ");
	}

#ifdef USE_EXISTING_CHUNK
	// Find closest chunk from original data and use that instead of average
	for (int g = 0; g < groups; g++)
	{
		int idx = 0;
		float distance = dist(dictionary + g * dimensions, chunkdata);
		for (int i = 0; i < chunks; i++)
		{
			float d = dist(dictionary + g * dimensions, chunkdata + i * dimensions);
			if (d < distance)
			{
				distance = d;
				idx = i;
			}
		}
		memcpy(dictionary + g * dimensions, chunkdata + idx * dimensions, dimensions);
	}
#endif

	/*
	for (int g = 0; g < groups; g++)
	{
		int p = g % group[g]; // pick random item
		for (int d = 0; d < dimensions; d++)
		{
			dictionary[g * dimensions + d] = chunkdata[index[groupofs[g] + p] * dimensions + d];
			//printf("%d ", dictionary[g * dimensions + d]);
		}
	}
	*/
	

	fwrite(dictionary, dimensions * 256, 1, outfile);
}

void map_indices()
{
	outdata = new unsigned char[chunks];
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
		outdata[i] = idx;
	}
/*
	for (int g = 0; g < groups; g++)
	{
		for (int i = 0; i < group[g]; i++)
		{
			outdata[index[groupofs[g] + i]] = g;
		}
	}
*/
/*
	printf("\n\n---\n\n");
	for (int i = 0; i < chunks; i++)
		printf("%d ", outdata[i]);
*/		
	fwrite(outdata, chunks, 1, outfile);
	fclose(outfile);
}

void reduce()
{
	analyze_group(0);
	while (groups < MAX_GROUPS)
	{
		int t = find_largest();
		int g = t / dimensions;
		int d = t % dimensions;
		split_group(g, d);
	}
}

void verify()
{
	unpackeddata = new unsigned char[datalen];
	for (int i = 0; i < chunks; i++)
	{
		for (int j = 0; j < dimensions; j++)
		{
			unpackeddata[i * dimensions + j] = dictionary[outdata[i] * dimensions + j];
		}
	}
	long long errsum = 0;	
	for (int i = 0; i < chunks * dimensions; i++)
	{
		errsum += abs(unpackeddata[i] - chunkdata[i]);
	}
	printf("Absolute error: %d, average error: %3.3f\n", (int)errsum, errsum / (double)(chunks * dimensions));
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
	drwav_init_file_write(&wav, outfilename, &format, NULL);
	drwav_write_pcm_frames(&wav, chunks * dimensions / channels, unpackeddata);
	drwav_uninit(&wav);
}

int main(int parc, char** pars)
{
	printf("LFAC encoder by Jari Komppa 2021 http://iki.fi/sol\n");
	if (parc < 2)
	{
		printf("Usage: %s infilename [outfilename]\n", pars[0]);
		return 0;
	}
	if (parc < 3)
	{
		sprintf(outfilename, "%s.sad", pars[1]);
	}
	else
		sprintf(outfilename, "%s", pars[2]);

	load_data(pars[1]); // 6.120 -> 5.423 (closest match)
	reduce();
	average_groups();
	map_indices();
	verify();
	save_data();
	return 0;
}