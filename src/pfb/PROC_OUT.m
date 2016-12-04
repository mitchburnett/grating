clearvars;
%%

N = 4000;
totalCh = 25;
coarseCh = 5;
numEl = 64;

fs = 303;
nfft = 32;
ntaps = 8;
subbands = coarseCh*numEl; % 5 channels processed * 64 el = 320 subbands

%% Load data
f = fopen('outfile.dat', 'r');
output = fread(f, 'float32');
fclose(f);

windows = floor((N - nfft*ntaps)/nfft);
size_slice = coarseCh*(numEl*2)*nfft;

rows = length(output)/size_slice; % This is not the same size as windows and need to figure out why.

CH = zeros(numEl*coarseCh, nfft);
for i = 1:rows
    slice = output((i-1)*size_slice+1:i*size_slice);
    tmp = reshape(slice, [coarseCh*(2*numEl) nfft]);
    
    
    for k = 1:coarseCh
        ch_slice = tmp((k-1)*numEl*2+1:k*numEl*2,:);
        el_re = ch_slice(1:2:end,:);
        el_im = ch_slice(2:2:end,:);
        el_spectra = abs(el_re + j*el_im).^2;
        
        CH((k-1)*numEl+1:k*numEl,:) = CH((k-1)*numEl+1:k*numEl,:) + el_spectra;
    end
end
% time average CH
CH = CH/rows;


%% Plot element and channel
el_idx = 1; % [1 64]
ch_idx = 1; % [1 5]

faxis = linspace(0, fs, nfft);
el_data = CH((ch_idx-1)*numEl+1:ch_idx*numEl,:);
el = el_data(el_idx, :);

figure(1);
subplot(321);
plot(faxis, 10*log10(el+.001)); grid on;
xlim([0, fs]);
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');

set(gca, 'XTick', [0:20:300]);
for i = 2:5
    ch_idx = i;
    el_data = CH((ch_idx-1)*numEl+1:ch_idx*numEl,:);
    el = el_data(el_idx, :);
    subplot(3,2,i);
    plot(faxis, 10*log10(el+.001)); grid on;
    xlim([0, fs]);
    xlabel('Frequency (kHz)');
    ylabel('Magnitude (dB)');

    set(gca, 'XTick', [0:20:300]);
end

