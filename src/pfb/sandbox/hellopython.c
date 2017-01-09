//My code file
#include <stdio.h>
#include <python2.7/Python.h>

int main()
{
    FILE* file;
    int argc;
    char * argv[3];

    argc = 3;
    argv[0] = "helloworld.py";
    argv[1] = "-m";
    argv[2] = "/tmp/targets.list";

    Py_SetProgramName(argv[0]);
    Py_Initialize();
    PySys_SetArgv(argc, argv);
    file = fopen("helloworld.py","r");
    PyRun_SimpleFile(file, "helloworld.py");
    Py_Finalize();

    return 0;
}