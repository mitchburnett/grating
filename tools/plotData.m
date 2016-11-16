%%
clearvars;
close all;

% f = fopen('../src/pfb/data/tone_15Mhz.dat', 'r');
% d = fread(f, 'schar');
% fclose(f);
% 
% f = fopen('../src/pfb/outfile.dat', 'r');
% output = fread(f, 'schar');
% fclose(f);
%% Constants
N = 4000;
fs = 250;
channels = 25;
numEl = 64;
subbands = numEl*2*channels;
nfft = 32;
ntaps = 8;

%% matlab gendata
n = 0:N-1;
E = zeros(subbands,N);
for l = 0:channels-1
    tmp1 = zeros(numEl*2,N);
    f = 10*l + 5;
    x_re = 127 * (0.1 * cos(2*pi*f*n/fs));
    x_im = 127 * (0.1 * sin(2*pi*f*n/fs));
    
    x_re = repmat(x_re, [numEl, 1]);
    x_im = repmat(x_im, [numEl, 1]);
    tmp1(1:2:end,:) = x_re;
    tmp1(2:2:end,:) = x_im;
    E(l*(numEl*2)+1:(l+1)*(numEl*2),:) = tmp1; 
end

d = reshape(E, [1 numEl*2*channels*N]).';


%%

freqChannels = zeros(channels, N);
E = zeros(numEl, N);

len = channels*numEl*2;

for i = 1:N
    timeSlice = d((i-1)*len+1:i*len);
    s_real = timeSlice(1:2:end);
    s_im = timeSlice(2:2:end);
    e0 = s_real(1:numEl:end);
    freqChannels(:,i) = e0;
end

%% plot spectrum

% combine signal into sum of sinusoids
s = sum(freqChannels);

% evaulte dft to plot
S = fft(s);

P2 = abs(S/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure;
faxis = fs*(0:N/2)/N;
plot(faxis,P1);

%%
PFB_CHANNELS = 5;
len = PFB_CHANNELS*numEl*2;

outputChannels = zeros(PFB_CHANNELS, N);
for i = 1:N
    timeSlice = output((i-1)*len+1:i*len);
    s_real = timeSlice(1:2:end);
    s_im = timeSlice(2:2:end);
    e0 = s_real(1:numEl:end);
    outputChannels(:,i) = e0;
end


% combine signal into sum of sinusoids
s = sum(outputChannels);

% evaulte dft to plot
S = fft(s);

P2 = abs(S/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure;
faxis = fs*(0:N/2)/N;
plot(faxis,P1);


