#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <functional>
#include "types.h"
#include "decompress.h"
#include "compress.h"

using namespace std;

#define worst_ratio 1.02
// Do not increase this limit, will lead to int overflow
const int MAX_INPUT_SIZE = int(INT_MAX / 8 / worst_ratio);

void PrintAbout() {
	printf("\n");
	printf("RIP compressor v2022.03.17");
	if (sizeof(void*) > 4) {
		printf(" (x64)\n");
	} else {
		printf(" (x86)\n");
	}
	printf("by Eugene Larchenko (https://gitlab.com/eugene77)\n");
	printf("\n");
}

void PrintUsage() {
	printf("Usage:\n");
	printf("  rip.exe [-d] <inputfilename> [<outputfilename>]\n");
	printf("    -d = decompress (default is compress)\n");
	printf("    -r = reverse bit order\n");
	printf("\n");
}

bool file_exists(const char* path)
{
	FILE* f = fopen(path, "r");
	if (f) {
		fclose(f);
	}
	return (f != NULL);
}

int main(int argc, char* argv[])
{
	int retCode = 10;
	string outPath = "\0";
	bool deleteResult = false;
	try
	{
		PrintAbout();

		char mode = 'c'; // compress
		int a = 1;
		char f_r = 0;

		while (a < argc) if (strcmp(argv[a], "-d") == 0) {
			mode = 'd'; a++;
		} else if (strcmp(argv[a], "-r") == 0) {
			f_r = 1; a++;
		} else break;
		if (mode == 'd' && f_r) {
			printf ("Error: can't decompress with reverse bit order\n");
			throw 2;
		}

		if (a >= argc) {
			PrintUsage();
			printf("Invalid arguments\n");
			throw 1;
		}
		string inPath = argv[a++];
		inPath += "\0";
		outPath = inPath;
		if (a < argc) {
			outPath = argv[a++];
			outPath += "\0";
		} else {
			outPath += (mode == 'c' ? ".rip\0" : ".unrip\0");
		}

		if (file_exists(outPath.c_str())) {
			// we don't want overwriting file. Don't waste time, abort now.
			printf("Error: output file already exists: %s\n", outPath.c_str());
			throw 2;
		}

		FILE* fIn = fopen(inPath.c_str(), "rb");
		if (!fIn) {
			printf("Error opening file %s\n", inPath.c_str());
			throw 2;
		}

		printf("Reading input file %s\n", inPath.c_str());
		vector<byte> data;
		while(true) {
			byte buf[4096];
			size_t q = fread(buf, 1, 4096, fIn);
			if (q <= 0) {
				break;
			}
			data.reserve(data.size() + q);
			for(size_t i=0; i<q; i++) {
				data.push_back(buf[i]);
			}
			if (data.size() > MAX_INPUT_SIZE) {
				throw exception("Input file is too large");
			}
		}
		
		fclose(fIn); fIn = NULL;

		std::function<void(vector<byte>&)> saveResult = [&](vector<byte>& result) {
			remove(outPath.c_str());
			deleteResult = true;
			FILE* fOut = fopen(outPath.c_str(), "wb");
			if (!fOut) {
				printf("Error creating output file %s\n", outPath.c_str());
				throw 3;
			}
			if (f_r) {
				byte rtable [256];
				size_t i;
				for (i = 0; i < 256; i++) rtable [i] = i << 7 | i << 5 & 0x40 | i << 3 & 0x20 | i << 1 & 0x10 |
					i >> 1 & 0x08 | i >> 3 & 0x04 | i >> 5 & 0x02 | i >> 7 & 0x01;
				for (i = 0; i < result.size(); i++) result.at (i) = rtable [result.at (i)];
			}
			size_t written = fwrite(result.data(), 1, result.size(), fOut);
			if (written != result.size()) {
				printf("Error writing output file\n");
				throw 3;
			}
			fclose(fOut);
		};

		printf("Processing\n");

		size_t datasize = data.size();
		if (mode == 'c')
		{
			// repeat the data; required for Z-function
			data.reserve(datasize * 2);
			for(size_t i=0; i < datasize; i++) {
				byte t = data[i]; data.push_back(t);
			}
		}

		vector<byte> result = mode=='c' 
			? compress(data, (int)datasize, saveResult)
			: decompress(data);

		//printf("time=%f\n", timeSpent);

		printf("Writing result to %s\n", outPath.c_str());
		saveResult(result);
		deleteResult = false;

		printf("All OK\n");
		retCode = 0;
	}
	catch (int code)
	{
		retCode = code;
	}
	catch (std::bad_alloc& e)
	{
		printf("ERROR: Couldn't allocate memory (%s)\n", e.what());
		retCode = 8;
	}
	catch (std::exception& e)
	{
		printf("ERROR: %s\n", e.what());
		retCode = 9;
	}

	if (deleteResult) {
		// delete incomplete compressed file
		remove(outPath.c_str());
	}

	return retCode;
}
