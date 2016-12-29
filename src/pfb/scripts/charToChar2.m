%% load data and define constants
clearvars;
N = 4000;
NUM_EL = 64;
FREQ = 25;
PFB_CHANNELS = 5;

f = fopen('../data/data.dat', 'r');
dataIn = fread(f, 'schar');
fclose(f);

%% Launch only the required number of threads

gridDim.x = N;
gridDim.y = PFB_CHANNELS;

blockDim.x = 1;
blockDim.y = NUM_EL;

blockIdx.x = 0:gridDim.x-1;
blockIdx.y = 0:gridDim.y-1;

threadIdx.x = 0:blockDim.x-1;
threadIdx.y = 0:blockDim.y-1;

dataOut.x = zeros(1,N*PFB_CHANNELS*NUM_EL);
dataOut.y = zeros(1,N*PFB_CHANNELS*NUM_EL);

IDX = zeros(N, FREQ);
MAP_IDX = zeros(N,PFB_CHANNELS);

select = 4;
fmin = select*PFB_CHANNELS;
fmax = fmin + PFB_CHANNELS-1;

for i = 0:gridDim.x-1
    for j = 0:gridDim.y-1
        absIdx = 2*blockDim.y*(blockIdx.x(i+1)*FREQ + (fmin+blockIdx.y(j+1))) + 2*threadIdx.y;
        IDX(i+1,j+1) = min(absIdx);
        mapIdx = blockDim.y*(blockIdx.x(i+1)*PFB_CHANNELS + blockIdx.y(j+1)) + threadIdx.y;
        MAP_IDX(i+1, j+1) = min(mapIdx);
        dataOut.x(mapIdx+1) = dataIn(absIdx+1); %real channel
        dataOut.y(mapIdx+1) = dataIn(absIdx+2); %imag channel
    end
end

%% 
len = PFB_CHANNELS*NUM_EL;
fs = 256;

outputChannels = zeros(PFB_CHANNELS, N);
for i = 1:N
    timeSlice_real = dataOut.x((i-1)*len+1:i*len);
    timeSlice_img = dataOut.y((i-1)*len+1:i*len);
    %s_real = timeSlice(1:2:end);
    s_real = timeSlice_real;
    %s_im = timeSlice(2:2:end);
    s_im = timeSlice_img;
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


