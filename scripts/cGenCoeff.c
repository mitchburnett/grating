/*
    -h  --help                 Display this usage information
    -n  --nfft <value>         Number of points in FFT
    -t  --taps <value>         Number of taps in PFB
    -w  --window <value>       Window to apply i.e "cheb-win", default: rect.
    -b  --sub-bands <value>    Number of sub-bands in data
    -d  --data-type <value>    Data type - "float" or "signedchar"
    -p  --no-plot              Do not plot coefficients
*/

#include <stdio.h>
#include <python2.7/Python.h>

int main()
{
	FILE* file;
	int argc;
	char * argv[7];
	char * fname = "grating_gencoeff.py";

	argc = 7;
	argv[0] = fname;
	argv[1] = "-n32";
	argv[2] = "-t8";
	argv[3] = "-whanning";
	argv[4] = "-b320";
	argv[5] = "-dfloat";
	argv[6] = "-p";

	Py_SetProgramName(argv[0]);
	Py_Initialize();
	PySys_SetArgv(argc, argv);
	file = fopen("grating_gencoeff.py", "r");
	PyRun_SimpleFile(file, "grating_gencoeff.py");
	Py_Finalize();

	return 0;

}