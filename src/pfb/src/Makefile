CC = nvcc

LIBRARIES = -lPython -lcufft

# Also known as a gencode, whenever there is a invalid device function its usually because the gencode doeesnt match the device architecture
# from what I understand by default if not specified nvcc picks the compute capability for the version of CUDA installed and so
# in general what you see in makefiles is an expansion of all pairs for multiple compatibiltiy and so 
ARCH_FLAG = -arch=sm_11


all: runPFB

runPFB: main.o pfb.o helper.o
	$(CC) -o ../bin/runPFB main.o pfb.o helper.o $(LIBRARIES)

main.o: main.cu helper.h pfb.h
	$(CC) -c main.cu $(ARCH_FLAG)

pfb.o: pfb.cu pfb.h
	$(CC) -c pfb.cu $(ARCH_FLAG)

helper.o: helper.c helper.h pfb.h
	$(CC) -c helper.c

clean:
	rm *.o ../bin/runPFB