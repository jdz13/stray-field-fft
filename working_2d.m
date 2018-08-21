
clc
clear

% Make up the variables to be convoluted.

mu0 = 4*pi*10^-7;

nx = 100; % Number in both positive and negative x direction
ny = 100; % +/- y


% The cell size needs to be defined to put the model into the proper scale.
CSx = 10^-3;
CSy = 10^-3;


CVol = CSx*CSy; % This allows for the proper moment distribution. 

% Also need to define the saturation magnetisation 
Ms = 10^6; % [A/m]

linex = (-nx:nx).*CSx;
liney = (-ny:ny).*CSy;

[X,Y] = meshgrid(linex,liney);
    
extra = 10^-10; % arbitrarily small number to give non-zeros
radialN = sqrt(X.^2 + Y.^2 ); % producing the distance from the centre
radialN = radialN+extra; % adding the extra gives non-zero for all values

% The 3D Green's function is G = 1/[4pi(r)]
% THIS CODE REQUIRES THE 2D GREENS FUNCTION, G = 1/[2pilog(r)].

Greens = 1./(2*pi().*log(radialN)); % Green's function for 2D
[gGx, gGy] = gradient(Greens); % The input has to be the gradient (2D)
maggradGreens = sqrt(gGx.^2 + gGy.^2); % Add the components (quatrature)

m = 2^nextpow2(length(linex)+length(liney)-1); % Find the next power of two to pad to. Most efficient at this value. Must be over double

% Now for a magnetic dirac delta matrix

Mag = zeros(2*nx+1,2*ny+1); % To produce the space invariant vector

% Add in the dirac delta function somewhere


spacer = 25;  % this controls the number of cells between adjacent Dirac delta functions

Mag(round(length(Mag(:,1))/2),round(length(Mag(:,1))/2)) = 1;
Mag(round(length(Mag(:,1))/2),round(length(Mag(:,1))/2)-1) = -1;
%Mag(round(length(Mag(:,1))/2)-spacer,round(length(Mag(:,1))/2),round(length(Mag(:,1))/2)) = 1;


Mag = Mag.*Ms.*CVol;


% Let's have a look at the three variables.


figure(4)
clf
subplot(3,3,1)
imagesc(linex,liney,Greens)
colorbar
title 'Green''s function G(r)'
subplot(3,3,2)
imagesc(linex,liney,maggradGreens)
colorbar
caxis([0,10^-3])
title ('Gradient of the Green''s function \nablaG(r)')
subplot(3,3,3)
imagesc(linex,liney,Mag)
colorbar
title 'Magnetic distribution (at their location)'


GreensFFT = (fftn(maggradGreens,[m,m]));

%obtain greens function gradients for all three orientations
GFFTx = (fftn(gGx,[m,m]));
GFFTy = (fftn(gGy,[m,m]));

MagFFT = (fftn(Mag,[m,m]));

MagxFFT = (fftn(Mag,[m,m]));
MagyFFT = (fftn(Mag,[m,m]));

subplot(3,3,5)
imagesc(abs(fftshift(GreensFFT)))
colorbar
caxis([0,0.2])
title 'Greens FFT'
subplot(3,3,6)
contour(abs(fftshift(MagFFT)))
colorbar
title 'Magnetic FFT'


CXY = convn(Mag,maggradGreens,'same').*mu0;
fXfY = MagFFT.*GreensFFT;

fX = GFFTx.*MagxFFT;
fY = GFFTy.*MagyFFT;


iFFT = (ifftn(fXfY)).*mu0;

iFFTx = (ifftn(fX)).*mu0;
iFFTy = (ifftn(fY)).*mu0;


subplot(3,3,7)
imagesc(linex,liney,(abs(fftshift(fXfY))))
colorbar
title 'FFT(X) x FFT(Y)'


in = [nx+1,ny+1];
en = [(3*nx+1),(3*ny+1)];

actual = iFFT(in(1):en(1),in(2):en(2));

subplot(3,3,8)
imagesc(linex,liney,(actual))
colorbar
caxis([-10^-11,10^-11])
colsca = get(colorbar, 'limits');
title 'iFFT(XY)'

subplot(3,3,9)
imagesc(linex,liney,((CXY)))
colorbar
caxis([-10^-11,10^-11])
title 'Result from convolution'

% Just as a bit of fun, label the whole plot in the empty area.
subplot(3,3,4)

th = title({'2D FFT vs Convolution',' ','with Gradient \nabla operator', ' ', 'JDZ'});
% get the position of the title
titlePos = get( th , 'position');
% change the x and y values
titlePos(1) = 0.4;
titlePos(2) = -0.1;

% update the position
set( th , 'position' , titlePos);

titleSiz = get(th, 'FontSize');
newtsize = 18;
set( th , 'FontSize', newtsize);
axis 'off'


%% 
% Plots to look at the field, both X and Y components as well as the total.
% This also has the Quiver plot on to visualise the field. Quiver plot is
% made using the X and Y components of the field calculated above. 

figure (10)
subplot(1,2,1)
imagesc(linex, liney, (iFFTx(in(1):en(1),in(2):en(2))))
colorbar
caxis([-10^-10,10^-10])
title 'Field X'

subplot(1,2,2)
imagesc(linex, liney, (iFFTy(in(1):en(1),in(2):en(2))))
colorbar 
caxis([-10^-10,10^-10])
title 'Field Y'

figure(11)
imagesc(linex, liney, sqrt((iFFTx(in(1):en(1),in(2):en(2))).^2+(iFFTy(in(1):en(1),in(2):en(2))).^2))
colorbar
caxis([0,10^-10])
title 'Total field - Dipole'
axis equal
hold on

removal12 = (abs(iFFTx) <= 10^-10); % has to be abs otherwise won't find the right values
removal13 = (abs(iFFTy) <= 10^-10);
removal11 = removal12.*removal13;
iFFTremx = iFFTx.*removal11;
iFFTremy = iFFTy.*removal11;

quiver(X, Y,(iFFTremx(in(1):en(1),in(2):en(2))),((iFFTremy(in(1):en(1),in(2):en(2)))))


%%
len = size(actual);
lines = actual(:,round(len(2)/2));
figure(12)
semilogy(liney,lines)
title 'Field at a distance plot'
ylabel 'Field (T)'
xlabel 'Distance (m)'