%%
clearvars;
close all;

N = 4000;
totalCh = 25;
coarseCh = 5;
numEl = 64;

fs = 303;
nfft = 32;
ntaps = 8;
%%
f = fopen('../bin/output/outfile.dat', 'r');
output = fread(f, 'float32');
fclose(f);

windows = floor((N - nfft*ntaps)/nfft);
len = coarseCh*numEl*2;


COARSE_CH = zeros(coarseCh*numEl*2, windows);





% %% plot the output from block mode pfb
% f = fopen('outfile.dat', 'r');
% output = fread(f, 'float32');
% fclose(f);
% 
% windows = floor((N - nfft*ntaps)/nfft);
% len = coarseCh*numEl*2;
% 
% element_spectrum = zeros(numEl, coarseCh*nfft);
% 
% for k = 1:numEl
%     fft_data = zeros(windows, nfft*coarseCh);
%     
%     for i = 1:windows
%         spectra = zeros(1,coarseCh*nfft);
%         for j = 1:nfft
%             slice = output((j-1)*len+1:j*len);
%             re = slice(1:2:end);
%             im = slice(2:2:end);
%             e0 = re(1:numEl:end)';
%             spectra(1, (j-1)*coarseCh+1:j*coarseCh) = e0;
%         end
%         fft_data(i,:) = spectra;
%     end
%     fft_data = mean(fft_data,1);
%     
%     element_spectrum(k,:) = fft_data;
% end

%% plot element spectrum
element =1;

f_axis = linspace(-fs/2, fs/2,nfft);

y = 10*log10(abs(element_spectrum(element,:)+.0001).^2);
figure(1);
plot(f_axis, y); grid on;
xlim([0,50/2]);
xlabel('Frequency MHz');
ylabel('Magnitude (dB)');

% %% use matlab to look at and compute spectrum
% f = fopen('data/tone_5kHz.dat', 'r');
% input = fread(f, 'schar');
% fclose(f);
% 
% freqChannels = zeros(channels, N);
% len = channels*numEl*2;
% 
% element_series_real = zeros(numEl,N);
% element_series_imag = zeros(numEl,N);
% 
% for i = 1:numEl
%     for j = 1:N
%         timeSlice = input((j-1)*len+1:j*len);
%         s_re = timeSlice(1:2:end);
%         s_im = timeSlice(2:2:end);
%         e_re = s_re(i:numEl:end);
%         e_re = sum(e_re);
%         e_im = s_im(i:numEl:end);
%         e_im = sum(e_im);
%         element_series_real(i,j) = e_re;
%         element_series_imag(i,j) = e_im;
%     end
% end
% %%
% d = element_series_real(1,:);
% len = ntaps*nfft;
% windows = floor((N - nfft*ntaps)/nfft);
% 
% F = zeros(windows,pfbChannels*nfft);
% 
% for i = 1:windows
%     slice = d( (i-1)*nfft+1 : (nfft*ntaps-1) +(i-1)*nfft);
%     f = fft(d,nfft);
%     F(i,:) = f;
% end
% 
% F = mean(F,1);
% f_axis = linspace(0,50/2,nfft*pfbChannels);
% y = 10*log10(abs(F).^2);
% figure(2);
% plot(f_axis,F); grid on;
% xlim([0,50/2]);
% xlabel('Frequency MHz');
% ylabel('Magnitude (dB)');

