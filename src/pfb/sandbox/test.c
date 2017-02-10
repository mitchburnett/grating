#include <stdio.h>
#include <string.h>

#define DEFAULT {3,2}
#define DEFAULT_PFB {32, 8, 320, 0, "hanning\0", "float\0", 1};

typedef struct {
	int a;
	int b;
} mystruct;

typedef struct {
	int nfft;
	int taps;
	int subbands;
	int select;
	char* window;
	char* dataType;
	int plot;
} params;

mystruct A = DEFAULT;
params pfbParams = DEFAULT_PFB;

int main() {
	(void) fprintf(stdout, "Hello!\n");
	(void) fprintf(stdout, "struct contains: a=%d, b=%d\n", A.a, A.b);

	char args[256] = {"-a %i -b %i"};
	sprintf(args, args, A.a, A.b);
	fprintf(stdout, "%s\n", args);

	char pfbArguments[256];// -b %i -d %s -s %i"};

	int i = 0;
	
	i = sprintf(pfbArguments, "-n %i -t %i -w %s -b %i -d %s -s %i", pfbParams.nfft,
										pfbParams.taps,
										pfbParams.window,
										pfbParams.subbands,
										pfbParams.dataType,
										pfbParams.select);

	(void) fprintf(stdout, "%i\n", i);

	fprintf(stdout, "%s\n", pfbArguments);

	if(pfbParams.plot) {
		strcat(pfbArguments, " -p");
	}

	fprintf(stdout, "%s\n", pfbArguments);

	const int count = 6;
	char* newArgs[count] = {};
	char temp[256] = {};

	newArgs[0] = "-n";
	sprintf(temp, "-%d", pfbParams.nfft);
	newArgs[1] = temp;
	fprintf(stdout, "%s %s\n", newArgs[0], newArgs[1]);


	return 0;
}