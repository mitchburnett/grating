%Simulated data and weights for beamformer

M = 38; % # of elements
M2 = 64;
F = 25; % # of frequency bins
n_beams = 14;
n_beams1 = n_beams/2;
N = 4000; % # of time samples
% theta = pi/6;
beam_angles_d = -45:(45 + 45)/N:45; % Steering vector angles
% beam_angles_d = beam_angles_d(1:end-1)+90; % Steering vector angles
beam_angles_d = ones(size(beam_angles_d))*90;
%beam_angles_d(1:2000) = 90;
%beam_angles_d(2000:end) = 85;
beam_angles = 0:(10 - 0)/4:10;
beam_angles = [-beam_angles(end:-1:2), beam_angles];
beam_angles = beam_angles(2:end-1)+90;
% beam_angles = ones(size(beam_angles))*90;
%theta_d = beam_angles(3);
c = 3*(10^8);
freq = 1550*(10^6);
lambda = c/freq;
d = lambda/2;
rf_freqs = freq - 100e6 + 500e3*(0:F-1);
%tau = d*(cosd(theta_d)/c); % Time delay corresponding to direction of arrival

weights_H = zeros(2*F*M2,n_beams);
weights = zeros(n_beams,F,M2);
weights_x = zeros(n_beams1,F,M2);
weights_y = zeros(n_beams1,F,M2);

for nb = 1:n_beams1
%     tau_beam = d*cosd(beam_angles(nb))/c;
    tau_beam_x = d*cosd(beam_angles(nb))/c;
    tau_beam_y = d*cosd(beam_angles(nb))/c;
%     phi = zeros(F,M);
    phi_x = zeros(F,M);
    phi_y = zeros(F,M);
    for f = 1:F
        for m = 0:M-1
%             phi(f,m+1) = m*2*pi*rf_freqs(f)*tau_beam; % Phase shift
            phi_x(f,m+1) = m*2*pi*rf_freqs(f)*tau_beam_x; % Phase shift
            phi_y(f,m+1) = m*2*pi*rf_freqs(f)*tau_beam_y; % Phase shift
        end
%         weights_x(nb,f,:) = exp(1j*phi_x(f,:));
%         weights_x(nb,f,:) = squeeze(weights_x(nb,f,:))./norm(squeeze(weights_x(nb,f,:)));
%         weights_y(nb,f,:) = exp(1j*phi_y(f,:));
%         weights_y(nb,f,:) = squeeze(weights_y(nb,f,:))./norm(squeeze(weights_y(nb,f,:)));
        weights_x(nb,f,1:(M/2)) = exp(1j*phi_x(f,1:(M/2)));
        weights_x(nb,f,((M/2)+1):M) = zeros(1,1,M/2);
        weights_x(nb,f,:) = squeeze(weights_x(nb,f,:))./norm(squeeze(weights_x(nb,f,:)));
        
        
        weights_y(nb,f,1:(M/2)) = zeros(1,1,M/2);
        weights_y(nb,f,((M/2)+1):M) = exp(1j*phi_y(f,1:(M/2)));
        weights_y(nb,f,:) = squeeze(weights_y(nb,f,:))./norm(squeeze(weights_y(nb,f,:)));
        
    end
    
    weights(1:n_beams1,:,:) = weights_x(:,:,:);
    weights((n_beams1+1):n_beams,:,:) = weights_y(:,:,:);
end

for nb = 1:n_beams
    weights_vec_C = reshape(squeeze(weights(nb,:,:)).', F*M2, 1); % Conjugate of the weights vector.
    
    % Interleave real and imaginary components (real then imaginary after)
    weights_C_real = real(weights_vec_C);
    weights_C_imag = imag(weights_vec_C);
    
    interleaved_w = zeros(2*F*M2,1);
    interleaved_w(1:2:end) = weights_C_real;
    interleaved_w(2:2:end) = weights_C_imag;
    weights_H(:,nb) = interleaved_w;
end

% Create metadata for weight file
offsets_el =  single([0,  0,    1,   1, 0,  -1,   -1]*2);
offsets_xel = single([0, -1, -0.5, 0.5, 1, 0.5, -0.5]*2);
offsets = [offsets_el; offsets_xel];
offsets = offsets(:);
cal_filename = '2016_06_13_16:58:04A.fits';
to_skip1 = 64 - length(cal_filename);
algorithm_name = 'Max Signal-to-Noise Ratio';
to_skip2 = 64 - length(algorithm_name);
xid = 3;

% Write to binary file
FID = fopen('weights_vec_C.bin','w');
% Write payload
fwrite(FID,single(weights_H(:)),'float');
% Write metadata
fwrite(FID,single(offsets),'float');
fwrite(FID,cal_filename, 'char*1');
if to_skip1 > 0
    fwrite(FID, char(zeros(1,to_skip1)));
end
fwrite(FID,algorithm_name, 'char*1');
if to_skip2 > 0
    fwrite(FID, char(zeros(1,to_skip2)));
end
fwrite(FID, uint64(xid), 'uint64');
fclose(FID);

% fileID2 = fopen('weights_vec_C.bin');
% weights_2 = fread(fileID2,[2*F*M*n_beams 1],'float');

s = ones(1,N); % Signal of interest amplitude. Simple test with no noise.

a = zeros(F*M2,N);
tmp_a = zeros(F,M2,N);

for h = 1:N
    tau_beam_d = d*cosd(beam_angles_d(h))/c;
    phid = zeros(F,M);
    for f = 1:F
        for m = 0:M-1
            phid(f,m+1) = m*2*pi*rf_freqs(f)*tau_beam_d; % Phase shift
        end        
    end
    tmp_a(:,1:(M/2),h) = exp(1j*phid(:,(M/2)+1:M));
    tmp_a(:,(M/2)+1:M,h) = exp(-1j*phid(:,(M/2)+1:M));
%     tmp_a(:,:,h) = exp(1j*phid(:,:));
    for f = 1:F
        tmp_a(f,:,h) = tmp_a(f,:,h)./norm(tmp_a(f,:,h));
    end
    a(:,h) = reshape(tmp_a(:,:,h).', F*M2, 1);
end

data_mat = zeros(F*M2,N);
for n = 1:N
    data_mat(:,n) = a(:,n)*s(n);
end
data_vec = reshape(data_mat,N*F*M2,1);

% Interleave real and imaginary components (real then imaginary after)
data_real = real(data_vec);
data_imag = imag(data_vec);

interleaved_d = zeros(2*F*M2*N,1);
interleaved_d(1:2:end) = data_real;
interleaved_d(2:2:end) = data_imag;

d_min = -1;
d_max = 1;
interleaved_chars = int8(ones(size(interleaved_d)));%int8(((interleaved_d - d_min)/(d_max - d_min) - 0.5)*256); 

% Write to binary file
fileID = fopen('data_vec.bin','w');
fwrite(fileID,single(interleaved_chars),'int8');
fclose(fileID);

% fileID3 = fopen('data_vec.bin');
% data_2 = fread(fileID3,[2*F*M*N 1],'float');
