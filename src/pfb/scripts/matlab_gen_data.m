clearvars;
N = 4000;
fs = 303;
totalChannels = 1;
numEl = 64;
subbands = numEl*2*totalChannels;
nfft = 32;
ntaps = 8;

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