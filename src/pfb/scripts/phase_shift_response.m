clearvars;
%%

N = 256*256;
totalCh = 1;
coarseCh = 1;
numEl = 1;

fs = 256;
nfft = 32;
ntaps = 8;
subbands = coarseCh*numEl; % 5 channels processed * 64 el = 320 subbands

windows = floor((N - nfft*ntaps)/nfft);
size_slice = coarseCh*(numEl*2)*nfft;

%% process multiple filter response outputs and average
% measured pfb response
files = 6;
response = zeros(files, windows/ntaps);
for p = 0:files-1
    % load the data
    filename = ['../bin/output/phase_out/phase_', int2str(p), '.dat'];
    f = fopen(filename, 'r');
    output = fread(f, 'float32');  
    fclose(f);

    %process outputs
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

            %CH((k-1)*numEl+1:k*numEl,:) = CH((k-1)*numEl+1:k*numEl,:) + el_spectra;
        end
    end
    % time average CH
    %CH = CH/windows;

    %% Plot the different response on top of eachother.
    figure(7); hold on; grid on;
    %numberOfResponses = 4;

    faxis_response = 0:fs/windows:fs - 1/windows;

    plot(faxis_response(1:ntaps:end), 10*log10(abs(XX(1:ntaps:end,17)).^2));
    xlim([min(faxis_response) max(faxis_response)]);
    title('PFB Filter Response');
    xlabel('Frequency (MHz)');
    ylabel('Gain (dB)');
    box on;
    
    response(p+1,:) = (XX(1:ntaps:end, 17));
end

R = mean(response);

figure(8);
plot(faxis_response(1:ntaps:end), 20*log10(abs(R))); grid on;
xlim([min(faxis_response) max(faxis_response)]);
