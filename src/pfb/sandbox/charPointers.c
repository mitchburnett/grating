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

	// char stackChars[2][7] = {"hello", "hi"}; // Could have also done that this way. As I understand it, this creates variables on the stack insted of malloc on the heap.

	// char* ptrPtr = stackChars[0];
	// ptrPtr++;

	// fprintf(stderr, "%s %s\n", stackChars[0], ptrPtr);

	char* progName = "./runPFB";
	int argCount = 11;
	if(pfbParams.plot) {
		argCount++;
	}

	char* arguments[32] = {};								// Initialize the two dimensional array. This is an array of pointers to char, i.e a char**.
	for(int i = 0; i < 32; i++){
		arguments[i] = (char*) malloc(256*sizeof(char*));	// The initialization of the variable arguments creates the pointers to char's, however, what they would point to has not been initialized i.e they are null, here they are initialized.
	}

	char temp[256] = {}; 									// a temp variable to store the various componetns to copy.

	arguments[0] = progName;
	arguments[1] = "-n\0";									// Had I not malloc'ed the array earlier this assignment would have worked since it is a static assignment. C allocatest the memoery for me to use it.
	sprintf(temp, "%d", pfbParams.nfft);					// Use the temp variable as a buffer for the string to be stored in arguments.
	strncpy(arguments[2], temp, 256);						// make the ptr copy from temp to arguments.

	arguments[3] = "-t\0";
	sprintf(temp, "%d", pfbParams.taps);
	strncpy(arguments[4], temp, 256);

	if(arguments[0] == NULL) {								// Since arguments had been initalized above, I had used this to prove to myself that the char* arguments[] array was initialzed but the arrayys it pointed to were null.
		fprintf(stderr, "NULL!!!\n");
		return 0;
	}

	fprintf(stderr, "%s %s %s\n", arguments[0], arguments[1], arguments[2]);


	// Now, after I have convinced myself that this works, and why, I have no verified that the temp variable
	// can done away with completely. I had tried using sprintf directly to arguments[], but I was getting an error.
	// this again, is due to the fact that I had not initalized the arrays that char* arguments[] could point to.

	arguments[5] = "-b";
	sprintf(arguments[6], "%d", pfbParams.subbands);

	arguments[7] = "-w";
	sprintf(arguments[8], "%s", pfbParams.window);

	arguments[9] = "-d";
	sprintf(arguments[10], "%s", pfbParams.dataType);

	if(pfbParams.plot){
		arguments[11] = "-p";
	}

	for(int i = 0; i <= argCount; i++) {
		fprintf(stderr, " %s", arguments[i]);
	}
	fprintf(stderr, "\n");

	return 0;
}