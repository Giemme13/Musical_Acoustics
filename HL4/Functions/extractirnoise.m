function [ir] = extractirnoise(x, y, nfft)
% EXTRACTIRNOISE 
% This function compute the transfer function of between a white noise
% source and the measurement point. The transfer function H is compute in
% the frequency domain as a simple deconvolution. The impulse response is
% obtained as the inverse fourier transform of H.
%
% Musical Acoustics Course
% Riccardo R. De Lucia
% 2018
% Mirco Pezzoli
% 2019-21

X = fft(x, nfft);
Y = fft(y, nfft);

%ifft computes a signal with a same number of samples as its input
ir = real(ifft(Y./X));
