clearvars;

% %%
% tap = 8;
% M = 10240;
% nfft = 32;
% 
% iter = 1;
% for k = 0:nfft-1
%     for j = 0:tap-1
%         idx(k+1,j+1) = j*M+k;
%         lin_idx(iter) = j*M+k;
%         iter = iter+1;
%     end
% end

%%
N = 4000;
fs = 303;
coarseCh = 5;
totalChannels = 25;
numEl = 64;

subbands = coarseCh*numEl;
nfft = 32;

e = 0;
c = 0;
t = 0:nfft:(N-1);

time_idx = e + (numEl)*c + (numEl*coarseCh)*t;

%% Generate some data
n = 0:N-1;
E = zeros(subbands,N) ;
for l = 0:totalChannels-1
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

d = reshape(E, [1 numEl*2*totalChannels*N]).';


%% GPU PFB implementation.
taps = 8;

blockDim.x = nfft; blockDim.y = 1; blockDim.z = 1;
gridDim.x = subbands*nfft/blockDim.x; gridDim.y =1; gridDim.z = 1;

blockIdx.x = 0:gridDim.x-1;

threadsPerBlock = blockDim.x*blockDim.y*blockDim.z;
threadIdx.x = 0:threadsPerBlock-1;

IDX = zeros(gridDim.x*taps, nfft);
for k = 0:gridDim.x-1;
    i = blockIdx.x(k+1)*gridDim.x + threadIdx.x;
    iNFFT = gridDim.x*blockDim.x;
    
    for j = 0:taps-1
        IDX(taps*k+1 + j,:) = j*iNFFT + i;
    end
end













