%clearvars;
%%

N = 512*512;
totalCh = 1;
coarseCh = 1;
numEl = 1;

fs = 512;
nfft = 64;
ntaps = 8;
subbands = coarseCh*numEl; % 5 channels processed * 64 el = 320 subbands

%% Load data
% measured pfb response
f = fopen('../bin/output/outfile.dat', 'r');
output = fread(f, 'float32');  
fclose(f);

% load designed prototype filter response
f = fopen('../bin/coeff_float_8_32_1.dat', 'r');
w = fread(f, 'float32');
fclose(f);

%% Plot designed pfb response
figure(1);
freqz(w);

radAxis = (-1:1/512:1-1/512);
W = fft(w, 1024);
W = fftshift(W);
figure(2); grid on;
idx = (0:length(W)-1);
plot(radAxis, 20*log10(abs(W)/abs(max(W)))); grid on;




%% Process pfb output samples
windows = floor((N - nfft*ntaps)/nfft);
size_slice = coarseCh*(numEl*2)*nfft;

rows = length(output)/size_slice; % This is not the same size as windows and need to figure out why.

CH = zeros(numEl*coarseCh, nfft);
XX = zeros(windows, nfft);
for i = 1:windows
    slice = output((i-1)*size_slice+1:i*size_slice);
    tmp = reshape(slice, [coarseCh*(2*numEl) nfft]);
    
    
    for k = 1:coarseCh
        ch_slice = tmp((k-1)*numEl*2+1:k*numEl*2,:);
        el_re = ch_slice(1:2:end,:);
        el_im = ch_slice(2:2:end,:);
        el_spectra = abs(el_re + 1j*el_im).^2;
        XX(i,:) = el_spectra; % save the spectra to plot filter response.
        
        CH((k-1)*numEl+1:k*numEl,:) = CH((k-1)*numEl+1:k*numEl,:) + el_spectra;
    end
end
% time average CH
CH = CH/windows;


%% Plot element and channel
el_idx = 1; % [1 64]
ch_idx = 1; % [1 5]

faxis = 0:fs/nfft:fs-1/nfft;
el_data = CH((ch_idx-1)*numEl+1:ch_idx*numEl,:);
el = el_data(el_idx, :);

figure(3); clf;
subplot(321);
plot(faxis, 10*log10(el+.001)); grid on;
xlim([min(faxis), max(faxis)]);
ylim([-5, 80]);
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Coarse Channel 1 - Processed output');
set(gca, 'xtick', [0:14]*20 + 5);

for i = 2:coarseCh
    ch_idx = i;
    el_data = CH((ch_idx-1)*numEl+1:ch_idx*numEl,:);
    el = el_data(el_idx, :);
    subplot(3,2,i);
    plot(faxis, 10*log10(el+.001)); grid on;
    xlim([min(faxis), max(faxis)]);
    xlabel('Frequency (kHz)');
    ylabel('Magnitude (dB)');
    ylim([-5, 80]);
    title(['coarse channel ' num2str(i)]);
    set(gca, 'xtick', [0:14]*20 + 5);
end

%% Plot the filter response.
figure(4); hold on; grid on;
%numberOfResponses = 4;

faxis_response = 0:fs/windows:fs - 1/windows;

plot(faxis_response, 20*log10(abs(XX(:,17))));
xlim([min(faxis_response) max(faxis_response)]);
title('PFB Filter Response');
xlabel('Frequency (MHz)');
ylabel('Gain (dB)');
box on;

