clearvars;
% Specturm data contains interleaved power spectrum data, XX, YY, XY, YX

% Open and read spec data
% pfb off.
f = fopen('spec_fft.dat', 'r');
fft_data = fread(f, 'float');
fclose(f);

% Open and read spec data -p
% pfb on
f = fopen('spec_pfb.dat', 'r');
pfb_data = fread(f, 'float');
fclose(f);

%open and read pfb spec data using dolph-cheb coeff
f = fopen('spec_cheb.dat', 'r');
cheb_data = fread(f, 'float');
fclose(f);

% filter specs
fs = 256e6;
nfft = 32;
series = nfft*4;
freqs = linspace(0,fs/2,nfft);
rows = length(fft_data)/series;

XX_fft = zeros(rows, nfft);
XX_pfb = zeros(rows, nfft);
XX_cheb = zeros(rows,nfft);

for i=1:rows    
    % Accumulate matrix of fft data
    chunk = fft_data((i-1)*series+1:i*series);
    XX_fft(i,:) = chunk(1:4:4*nfft);
    %Accumulate matrix of pfb data
    chunk = pfb_data((i-1)*series+1:i*series);
    XX_pfb(i,:) = chunk(1:4:4*nfft);
    %Accumulate matrix of cheb data
    chunk = cheb_data((i-1)*series+1:i*series);
    XX_cheb(i,:) = chunk(1:4:4*nfft);
end

%filter response comparison plot - plots 8 responses
figure; hold on;
for i=1:8:nfft
    plot(linspace(0,fs/2,length(XX_pfb)), 10*log10(XX_pfb(:,i)), '-b');
    plot(linspace(0,fs/2,length(XX_fft)), 10*log10(XX_fft(:,i)), '-r');
end

%configure axis
grid on;
xlim([0, fs/2]); ylim([0, 80]);
xtick = linspace(0, fs/2, nfft+1);
set(gca, 'XTick', xtick(1:8:end), 'FontSize', 12);
title('PFB Response', 'FontSize', 20);
xlabel('Frequency (MHz)', 'FontSize', 16);
legend('PFB', 'FFT');

% show accumulated power over all freq.
XX_pfb_mean = mean(XX_pfb,1);
figure;
plot(freqs, 10*log10(XX_pfb_mean));