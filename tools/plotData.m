%%
clearvars;

f = fopen('../src/pfb/data.dat', 'r');
d = fread(f, 'schar');
fclose(f);

f = fopen('../src/pfb/outfile.dat', 'r');
output = fread(f, 'schar');
fclose(f);
%%
N = 4000;       % time samples
channels = 25;  % frequency channels
numEl = 64;     % number of elements
fs = 256;       % sample rate MHz

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


