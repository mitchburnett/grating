clearvars;
% Open spec data - pfb off.
f = fopen('fft_spec.dat', 'r');
d = fread(f, 'float');
fclose(f);


nfft = 32;
series = nfft*4;
freqs = linspace(0,128e6,nfft);
rows = length(d)/series;

XX = zeros(rows, nfft);

for i=1:rows
    
    chunk = d((i-1)*series+1:i*series);
    XX(i,:) = chunk(1:4:4*nfft);
end
keyboard;
XX = mean(XX,1);

figure;
plot(freqs, XX);