clearvars;
fname = '../bin/output/fftIn_no_pfb.dat';
f = fopen(fname, 'r');
noPfb = fread(f, 'float32');  
fclose(f);

fname = '../bin/output/fftIn_pfb.dat';
f = fopen(fname, 'r');
pfb = fread(f, 'float32');
fclose(f);

fname = '../bin/coeff_float_8_32_320.dat';
f = fopen(fname, 'r');
coeff = fread(f, 'float32');
fclose(f);

figure(1);
plot(coeff);
figure(2);
plot(pfb);
figure(3);
plot(noPfb);
