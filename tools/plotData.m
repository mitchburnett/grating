%%
clearvars;

f = fopen('data.dat', 'r');
d = fread(f, 'schar');
fclose(f);
%%
N = 4000;
channels = 25;
numEl = 64;

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


