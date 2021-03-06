clearvars;

N = 4000;
NUM_EL = 64;
FREQ = 25;
PFB_CHANNELS = 5;

%% load data
f = fopen('../data/full_tone_set_kHz.dat', 'r');
dataIn = fread(f, 'schar');
fclose(f);

%% 
blockDim.x = 1;
blockDim.y = NUM_EL*2;

gridDim.x = 4000;
gridDim.y = FREQ;

blockIdx.x = 0:gridDim.x-1;
blockIdx.y = 0:gridDim.y-1;

threadsPerBlock = blockDim.x*blockDim.y;

threadIdx.y = 0:threadsPerBlock-1;
threadIdx.x = 1;

dataOut = zeros(N*PFB_CHANNELS*NUM_EL*2,1);

IDX = zeros(gridDim.y, gridDim.x);
MAP_IDX = zeros(gridDim.y/5, gridDim.x);
select = 0;
fmin = select*PFB_CHANNELS;
fmax = fmin + PFB_CHANNELS-1;

for i = 0:gridDim.x-1
    for j = 0:gridDim.y-1
        absIdx = blockDim.y * (blockIdx.x(i+1) * gridDim.y + blockIdx.y(j+1)) + threadIdx.y;
        IDX(j+1,i+1) = min(absIdx);
        if(j >= fmin && j <= fmax)
            idx = mod(j,PFB_CHANNELS);
            midx = blockDim.y * (blockIdx.x(i+1) * gridDim.y/PFB_CHANNELS + blockIdx.y(idx+1)) + threadIdx.y;
            MAP_IDX(idx+1,i+1) = min(midx);
            dataOut(midx + 1) = dataIn(absIdx + 1);
        end
    end
end


%% 
len = PFB_CHANNELS*NUM_EL*2;
fs = 303;

outputChannels = zeros(PFB_CHANNELS, N);
for i = 1:N
    timeSlice = dataOut((i-1)*len+1:i*len);
    s_real = timeSlice(1:2:end);
    s_im = timeSlice(2:2:end);
    e0 = s_real(1:NUM_EL:end);
    outputChannels(:,i) = e0;
end

s = sum(outputChannels);

% evaulte dft to plot
S = fft(s);

P2 = abs(S/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure;
faxis = fs*(0:N/2)/N;
plot(faxis,P1);

%% Launch only the required number of threads

gridDim.x = N;
gridDim.y = PFB_CHANNELS;

blockDim.x = 1;
blockDim.y = NUM_EL*2;

blockIdx.x = 0:gridDim.x-1;
blockIdx.y = 0:gridDim.y-1;

threadIdx.x = 0:blockDim.x-1;
threadIdx.y = 0:blockDim.y-1;

dataOut = zeros(1,N*PFB_CHANNELS*NUM_EL*2);

IDX = zeros(N, FREQ);
MAP_IDX = zeros(N,PFB_CHANNELS);

select = 1;
fmin = select*PFB_CHANNELS;
fmax = fmin + PFB_CHANNELS-1;

for i = 0:gridDim.x-1
    for j = 0:gridDim.y-1
        absIdx = blockDim.y*(blockIdx.x(i+1)*FREQ + (fmin+blockIdx.y(j+1))) + threadIdx.y;
        IDX(i+1,j+1) = min(absIdx);
        mapIdx = blockDim.y*(blockIdx.x(i+1)*PFB_CHANNELS + blockIdx.y(j+1)) + threadIdx.y;
        MAP_IDX(i+1, j+1) = min(mapIdx);
        dataOut(mapIdx+1) = dataIn(absIdx+1);
    end
end

%% 
len = PFB_CHANNELS*NUM_EL*2;
fs = 256;

outputChannels = zeros(PFB_CHANNELS, N);
for i = 1:N
    timeSlice = dataOut((i-1)*len+1:i*len);
    s_real = timeSlice(1:2:end);
    s_im = timeSlice(2:2:end);
    e0 = s_real(1:NUM_EL:end);
    outputChannels(:,i) = e0;
end

s = sum(outputChannels);

% evaulte dft to plot
S = fft(s);

P2 = abs(S/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure;
faxis = fs*(0:N/2)/N;
plot(faxis,P1);


