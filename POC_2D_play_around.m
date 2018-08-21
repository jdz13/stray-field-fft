% A 2D play around to make sure that we understand the Hermitic properties
% of the FFT. THis will allow us to see when the inputs are put in properly
% and how to do this for a full system.

% If an input to an FFT is real - the FFT produces complex results  - but
% with Hermitic properties - the negative frequency is equal to its complex
% conjugate P(-f) = (P(f))*. 
% 
% Taking the ifft of a Hermitic matrix yeilds real results.
%
% You can multiply Hermitic FFT's - which corresponds to a convolution in
% real space/time - to yield another Hermitic function in Fourier space. 
% When this is inverse Fourier transformed back to real space/time, the
% results are again completely real. 
%
% If any ifft result comes back as complex - something is wrong. This is to
% do with the boundary conditions. The input has to be padded with zeros to
% ensure that the periodicity is correct.
%

clear
clc

% Make the input matrixies
final = 10; % The final value in space/time
spacer = 0.01; % Space between adjacent cells
x = (0:0.01:final-spacer);
y = (0:0.01:final-spacer);

[X1,Y1] = meshgrid(x,y);

X = cos(2*pi*((X1*cos(pi/4))+Y1*sin(pi/4))) ;%  + cos(2*pi*5*XYgrid) ; % Creates an array that
% changes only in X. This will help to see.
caption  =('X = cos(2*pi*((X1*cos(pi/4))+Y1*sin(pi/4)))');


Y1 = Y1.*0; Y1(480:520,480:520) = 1; % Rectangular 'box' function. Could be anthing though. 

m = 2^nextpow2(length(x)); % Find the next power of two to pad to. Most efficient at this value.

fftX = fftn(X, [m,m]); % Taking the FFT of the first function (X).
ifftX = ifftn(fftX, [m,m]); % Taking the IFFT of the first function (X).

fftY = fftn(Y1, [m,m]); % Taking the FFT of the second function (Y).
ifftY = ifftn(fftY, [m,m]);% Taking the IFFT of the first function (Y).

fftconv = fftX.*fftY; % Multiplying the functions in Fourier space - convoluting in real space/time. 
ifftconv = ifftn(fftconv, [m,m]); % Inversing the results. 

% The result from this should be real. 



figure(17)

subplot(3,1,1)
imagesc(x,y,X)
title (caption);
xlabel 'x'
ylabel 'y'
colorbar

subplot(3,1,2)
imagesc(abs(fftshift(fftX)))
caxis([0,100])
title 'FFT results';
xlabel 'K space in x (arb units)'
ylabel 'K space in y (arb units)'
colorbar

subplot(3,1,3)
imagesc(x,y, ifftX(1:length(x),1:length(y)))
title 'Returned iFFT - all real';
xlabel 'x'
ylabel 'y'
colorbar


figure(18)

subplot(3,1,1)
imagesc(x,y,Y1)
title 'Input function - Dirac delta';
xlabel 'x'
ylabel 'y'
colorbar

subplot(3,1,2)
imagesc(abs(fftshift(fftY)))
caxis([0,100])
title 'FFT results';
xlabel 'K space in x (arb units)'
ylabel 'K space in y (arb units)'
colorbar

subplot(3,1,3)
imagesc(x,y, ifftY(1:length(x),1:length(y)))
title 'Returned iFFT - all real';
xlabel 'x'
ylabel 'y'
colorbar
%%
figure(19)

subplot(2,1,1)
imagesc(abs(fftconv))
title 'FFTx \otimes FFTy';
xlabel 'x'
ylabel 'y'
colorbar
caxis([0,10])

subplot(2,1,2)
imagesc(ifftshift(ifftconv))
title 'IFFT of both FFT''s multiplied';
xlabel 'K space in x (arb units)'
ylabel 'K space in y (arb units)'
colorbar
caxis([0,1000])

