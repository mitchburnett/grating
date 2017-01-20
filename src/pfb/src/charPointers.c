#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define DEFAULT_PFB {32, 8, 320, 0, "hanning\0", "float\0", 1};

typedef struct {
	int nfft;
	int taps;
	int subbands;
	int select;
	char* window;
	char* dataType;
	int plot;
} params;

params pfbParams = DEFAULT_PFB;

int main() {

	(void) fprintf(stdout, "Hello!\n");

	char* progName = "./runPFB";
	int argCount = 11;
	if(pfbParams.plot) {
		argCount++;
	}

	char* arguments[32] = {};
	for(int i = 0; i < 32; i++){
		arguments[i] = (char*) malloc(256*sizeof(char*));
	}

	char temp[256] = {};

	arguments[0] = "-n\0";
	sprintf(temp, "%d", pfbParams.nfft);
	strncpy(arguments[1], temp, 256);

	arguments[2] = "-t\0";
	sprintf(temp, "%d", pfbParams.taps);
	strncpy(arguments[3], temp, 256);

	if(arguments[0] == NULL) {
		fprintf(stderr, "NULL!!!\n");
		return 0;
	}

	fprintf(stderr, "%s %s %s %s\n", arguments[0], arguments[1], arguments[2], arguments[3]);


	return 0;
}