% open data
f = fopen('single_tone_5Mz.dat', 'r');
d = fread(f, 'schar');
fclose(f);

N = 4000;
totalChannels = 25;
numEl = 64;
nfft = 32;
ntaps = 8;

windows = ceil(N - nfft*ntaps)/nfft;

coarseChannels = 5;
fineChannels = nfft*coarseChannels;
subbands = coarseChannels * numEl;

select = 0;
c_min = select*coarseChannels;
c_max = select+4;

for c = c_min:c_max
    
    % TODO : Compute index for total data
    
    totalSpectra = zeros(numEl,nfft);
    
    for e = 1:numEl
        element = zeros(1,N);
        % Effectively loop on time.
        %   idx needs to becomputed in reference to the entire data set.
        %   This is why totalChannels is also used for the jump factor
        %   because instead of taking the data out just jump over it.
        element(1,:) = d(idx: coarseChannels*numEl*2 : end );
            
        % Have acquired the time series for a single element in one of the
        % coarse channels of interst.

        % perform fine channelization and average over time for spectrometer
        % output to verify results
        spectra = zeros(windows,nfft);
        for n = 1:windows
            S = fft(element((n-1)*nfft: (nfft*ntaps-1) + nfft*(n-1)),nfft); % need to add the imag elements
            spectra(n,:) = abs(S).^2;
        end
        %average the spectra
        spectra = mean(spectra,1);
        
        totalSpectra(e,:);
    end
    
    %Channels will be a (coarseChannels*numEl)*nfft matrix where each jump
    %of numEl represents a channel.
    c_idx = (c-c_min); % maps the c idx down to place in a matrix.
    CHANNELS( c_idx*numEl + 1: numEl*(c_idx+1), : ) = totalSpectra;
    
    
end




