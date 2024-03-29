function [phase,pow] = multiphasevec3(f,S,Fs,width, silent)
% FUNCTION [phase,pow] = multiphasevec3(f,S,Fs,width, silent)
%
% Returns the phase and power as a function of time for a range of
% frequencies (f).
%
% INPUT ARGS:
%   f = [2 4 8];   % Frequencies of interest
%   S = dat;       % Signal to process (trials X samples)
%   Fs = 256;      % Sampling frequency
%   width = 6;     % Width of Morlet wavelet (>= 5 suggested).
%   silent = 0     % Silence the warning when using this function with vectors
%
% OUTPUT ARGS:
%   phase- Phase data [trials,freqs,time]
%   power- Power data [trials,freqs,time]
%   Sfft - descrete Fourier transform of S
%
% Edit 07/13/2012: changed for loop in line 42 to parfor in order to implement parallel program. Make sure matlabpool has opened some workers
% For S- time is assumed to be the second dimensions (rows) (see line 54)
%
% REVISION HISTORY
%   08/2018 MST - Add warning and silent parameter, and vector reshape (to avoid unintentional missuse)

if isvector(S)
    S = S(:)';
end
if nargin < 5
    silent = 0;
end

nF = length(f);
nS = size(S);

dt = 1/Fs;
st = 1./(2*pi*(f/width));

% get the Morlet's wavelet for each frequency
curWaves = arrayfun( @(i) morlet( f(i), -3.5*st(i):dt:3.5*st(i), width ), 1:nF, 'UniformOutput', false );
nCurWaves = cellfun( @(w) length(w), curWaves );

Lys = nS(2) + nCurWaves - 1;    % length of convolution of S and curWaves{i}
% Ly2s = pow2(nextpow2(Lys));     % next power of two (for fft)
% Changing this to make Ly2s a vector of length nF

for i=1:length(Lys)
  Ly2s(i)=pow2(nextpow2(Lys(i))); % next power of two (for fft)
end
ind1 = ceil(nCurWaves/2);       % start index of signal after convolution

pow = zeros(nS(1), nF, nS(2));
phase = zeros(nS(1), nF, nS(2));

for i = 1:nF
    Ly2 = Ly2s(i);
    
    %%% Perform the convolution of curWaves{i} with every row of S
    % take the fft of S and curWaves{i}, multiply them, and take the ifft
    
    % Sfft = fft(S,Ly2,2);
    % curWaveFFT = fft(curWaves{i},Ly2);
    % Y = bsxfun(@times, Sfft, curWaveFFT);
    % y1 = ifft(Y,Ly2,2);
    
    % (EH - it's much quicker to do it in one line)
    y1 = ifft( bsxfun( @times, fft(double(S),Ly2,2), fft(curWaves{i},Ly2) ) ,Ly2,2);

    y1 = y1( :, ind1(i):(ind1(i)+nS(2)-1) );
    
    % find phase and power (do it inside this loop to save memory)
    pow(:,i,:) = abs(y1).^2;
    phase(:,i,:) = angle(y1);
end

if ~silent
    sz = size(pow);
    if size(pow,1) == 1
        fprintf(['\nWarning::multiphasevec3: size(pow)=[%d %d %d].\n', ...
            'This is a reminder that this function expects [trials x samples].\n', ...
            'You may want to squeeze(pow) now, since it contans [trials=1 x freq x samples].\n\n'], ...
            sz(1), sz(2), sz(3));
    end
end


