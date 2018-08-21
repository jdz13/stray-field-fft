% A 1D play around to make sure that we understand the Hermitic properties
% of the FFT. THis will allow us to see when the inputs are put in properly
% and hopefully shed some light on how to do it for 2D and higher
% dimensions. 


clc
clear


Fs = 1000; % [Hz]
finaltime = 1; % [s]
t = 0:1/Fs:finaltime-(1/Fs); % Time vector in increments of 1/Fs

f = 5; % [Hz]
x = 1 + sin(2*pi*f*t) + sin(2*pi*100*t); % Generate the input signal. Make slightly more complicated - more than one frequency and DC part do these can be tested.
nfft = 2.^nextpow2(length(x)); % Calculate the next power of two, to pad the FFT to.
 
X = fft(x,nfft); % Taking the FFT.
Xgraph = X(1:round(nfft/2)); % Taking only the positive frequencies so that it's easier to look at. 

magX = abs(Xgraph)/length(x); % Normalising (energy dist) by the number of samples. 

fv = (0:nfft/2-1)*Fs/nfft; % Generating the frequency axis values. 

inv = ifftn(X); % Looking at the IFFT. Already a power 2 so nfft not needed.

% look at the input signal
figure(19)
subplot(3,1,1)
plot(t,x)
title 'Input signal'
xlabel 'Time [s]'
ylabel 'Amplitude'

% Show the FFT - we're only really interested in the positive side.
subplot (3,1,2)
plot(fv, (magX))
title 'FFT results'
xlabel 'Frequency [Hz]'
ylabel 'FFT amplitude'

% Show the iFFT.
% If the FFT is Hermitic, and the appropriate BC's have been applied, then
% this should be the same as the input (when the zeros are taken out). 
subplot(3,1,3)
plot(t, inv(1:length(t)))
title 'IFFT results'
xlabel 'Time [s]'
ylabel 'IFFT Amplitude'