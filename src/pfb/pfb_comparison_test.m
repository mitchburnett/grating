close all; clearvars;
% open data
f = fopen('full_tone_set.dat', 'r');
d = fread(f, 'schar');
fclose(f);

N = 4000;
fs = 256;
totalChannels = 25;
numEl = 64;
nfft = 32;
ntaps = 8;

windows = ceil(N - nfft*ntaps)/nfft;
coarseChannels = 5;
fineChannels = nfft*coarseChannels;
subbands = coarseChannels * numEl;

select = 1;
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
            S_re = fft(hanning(nfft*ntaps)'.*element_re((n-1)*nfft + 1: (nfft*ntaps) + nfft*(n-1)),nfft); % need to add the imag elements
            S_im = fft(hanning(nfft*ntaps)'.*element_im((n-1)*nfft + 1: (nfft*ntaps) + nfft*(n-1)),nfft); % need to add the imag elements
            
            spectra(n,:) = abs(S_re + j*S_im).^2;
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
% faxis = linspace(0,128/25,nfft);
% plot(faxis, 10*log10(CHANNELS(1,:)+.0001));
% xtick = linspace(0,5.12/nfft,nfft+1);
% xlabel('Freq (Mhz)');
% ylabel('Magnitude dB');

ch = 0;
absCh = ch+c_min;
faxis = linspace(absCh*10.24,(absCh+1)*10.24,nfft);
plot(faxis, 10*log10(CHANNELS(1+64*ch,:)+.0001));

xlabel('Freq (Mhz)');
ylabel('Magnitude dB');
