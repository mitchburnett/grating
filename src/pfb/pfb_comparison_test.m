close all; clearvars;
% constants
N = 4000;
fs = 303e3; % 303 Khz
totalChannels = 25;
numEl = 64;
subbands = totalChannels*numEl;
nfft = 256;
ntaps = 8;

% load data
loadFromFile = 1;
if loadFromFile
    fname = 'data/tone_5kHz.dat';
    f = fopen(fname, 'r');
    d = fread(f, 'schar');
    fclose(f);
else
    n = 0:N-1;
    E = zeros(subbands,N) ;
    for l = 0:totalChannels-1
        tmp1 = zeros(numEl*2,N);
        f = 10e3*l + 5e3; % base f of 5Khz and jump by 10
        x_re = 127 * (0.1 * cos(2*pi*f*n/fs));
        x_im = 127 * (0.1 * sin(2*pi*f*n/fs));

        x_re = repmat(x_re, [numEl, 1]);
        x_im = repmat(x_im, [numEl, 1]);
        tmp1(1:2:end,:) = x_re;
        tmp1(2:2:end,:) = x_im;
        E(l*(numEl*2)+1:(l+1)*(numEl*2),:) = tmp1; 
    end
    d = reshape(E, [1 numEl*2*totalChannels*N]).';
end

% start
windows = ceil((N - nfft*ntaps)/nfft);
coarseChannels = 5;
fineChannels = nfft*coarseChannels;
subbands = coarseChannels * numEl;

select = 0;
c_min = select*coarseChannels;
c_max = c_min+4;

CHANNELS = zeros(subbands,nfft);

for c = c_min:c_max
    
    totalSpectra = zeros(numEl,nfft);
    
    for e = 1:numEl
        element_re = zeros(1,N);
        element_im = zeros(1,N);
        % Effectively loop on time.
        %   idx needs to becomputed in reference to the entire data set.
        %   This is why totalChannels is also used for the jump factor
        %   because instead of taking the data out just jump over it.
        idx = c*numEl*2 + 2*e - 1;
        element_re(1,:) = d(idx: totalChannels*numEl*2 : end );
        element_im(1,:) = d(idx+1 : totalChannels*numEl*2 : end);
            
        % Have acquired the time series for a single element in one of the
        % coarse channels of interst.

        % perform fine channelization and average over time for spectrometer
        % output to verify results
        spectra = zeros(windows,nfft);
        for n = 1:windows
            S_re = fft(element_re((n-1)*nfft + 1: (nfft*ntaps) + nfft*(n-1)),nfft); % need to add the imag elements
            S_im = fft(element_im((n-1)*nfft + 1: (nfft*ntaps) + nfft*(n-1)),nfft); % need to add the imag elements
            
            spectra(n,:) = fftshift(abs(S_re + j*S_im).^2); % fftshift is important since it is plotting 0 - 2pi not -pi to pi
        end
        %average the spectra
        spectra = mean(spectra,1);
        
        totalSpectra(e,:) = spectra;
    end
    
    %Channels will be a (coarseChannels*numEl)*nfft matrix where each jump
    %of numEl represents a channel.
    c_idx = (c-c_min); % maps the c idx down to place in a matrix.
    CHANNELS(c_idx*numEl + 1: numEl*(c_idx+1), : ) = totalSpectra;
    
    
end
%%

fs_n = fs/1e3; %normalized fs
coarseChannelBinWidth = fs_n/totalChannels;
fineChannelBinWidth = coarseChannelBinWidth/nfft;

ch = 0;
absCh = ch+c_min;
%faxis = linspace(absCh*coarseChannelBinWidth, (absCh+1)*coarseChannelBinWidth,nfft);
%faxis = absCh*10:fineChannelBinWidth:(absCh+1)*10-fineChannelBinWidth;
faxis = linspace(-fs_n/2, fs_n/2, nfft);
subplot(321);
plot(faxis, 10*log10(CHANNELS(1+numEl*ch,:)+.0001)); grid on;
title('Coarse Channel 1');
xlabel('Freq (Khz)');
ylabel('Magnitude dB');
xlim([-fs_n/2, fs_n/2]);
set(gca, 'XTick', [-fs_n/2:25:fs_n/2]);

ch = 1;
absCh = ch+c_min;
%faxis = linspace(absCh*coarseChannelBinWidth, (absCh+1)*coarseChannelBinWidth,nfft);
%faxis = absCh*10:fineChannelBinWidth:(absCh+1)*10-fineChannelBinWidth;
subplot(322);
plot(faxis, 10*log10(CHANNELS(1+numEl*ch,:)+.0001)); grid on;
title('Coarse Channel 2');
xlabel('Freq (Khz)');
ylabel('Magnitude dB');
xlim([-fs_n/2, fs_n/2]);
set(gca, 'XTick', [-fs_n/2:25:fs_n/2]);

ch = 2;
absCh = ch+c_min;
%faxis = linspace(absCh*coarseChannelBinWidth, (absCh+1)*coarseChannelBinWidth,nfft);
%faxis = absCh*10:fineChannelBinWidth:(absCh+1)*10-fineChannelBinWidth;
subplot(323);
plot(faxis, 10*log10(CHANNELS(1+numEl*ch,:)+.0001)); grid on;
title('Coarse Channel 3');
xlabel('Freq (Khz)');
ylabel('Magnitude dB');
xlim([-fs_n/2, fs_n/2]);
set(gca, 'XTick', [-fs_n/2:25:fs_n/2]);

ch = 3;
absCh = ch+c_min;
%faxis = linspace(absCh*coarseChannelBinWidth, (absCh+1)*coarseChannelBinWidth,nfft);
%faxis = absCh*10:fineChannelBinWidth:(absCh+1)*10-fineChannelBinWidth;
subplot(324);
plot(faxis, 10*log10(CHANNELS(1+numEl*ch,:)+.0001)); grid on;
title('Coarse Channel 4');
xlabel('Freq (Khz)');
ylabel('Magnitude dB');
xlim([-fs_n/2, fs_n/2]);
set(gca, 'XTick', [-fs_n/2:25:fs_n/2]);

ch = 4;
absCh = ch+c_min;
%faxis = linspace(absCh*coarseChannelBinWidth, (absCh+1)*coarseChannelBinWidth,nfft);
%faxis = absCh*10:fineChannelBinWidth:(absCh+1)*10-fineChannelBinWidth;
subplot(325);
plot(faxis, 10*log10(CHANNELS(1+numEl*ch,:)+.0001)); grid on;
title('Coarse Channel 5');
xlabel('Freq (Khz)');
ylabel('Magnitude dB');
xlim([-fs_n/2, fs_n/2]);
set(gca, 'XTick', [-fs_n/2:25:fs_n/2]);