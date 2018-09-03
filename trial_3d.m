clc 
clear

%test line 2

% start to think about the grid. Start with thinking about the variables.
%------------------------------------------------------------------------

rangex = 10^-3; % time/space limit in x [s] or [m]
rangey = 10^-3; % time/space limit in y [s] or [m]
rangez = 10^-3; % time/space limit in z [s] or [m]

nx = 30; % number of cells to spread over in both positive and negative direction
ny = 30;
nz = 30;

extra = 10^-10;

linex = (-rangex:2*rangex/(2*nx):rangex)+extra;
liney = (-rangey:2*rangey/(2*ny):rangey)+extra;
linez = (-rangez:2*rangez/(2*nz):rangez)+extra;

[X,Y,Z] = meshgrid(linex,liney,linez);


mx = 2^nextpow2(3*length(X)+1);
my = 2^nextpow2(3*length(Y)+1);
mz = 2^nextpow2(3*length(Z)+1);

radialN = sqrt(X.^2 + Y.^2 + Z.^2);

% Amazing. Space is sorted. Now let's turn this into something useful. 
% The 2D Green's function is G = -1/2pi log(r)

Greens = log(radialN)./(-2*pi());

[gGx,gGy,gGz] = gradient(Greens);

GG = sqrt(gGx.^2 + gGy.^2 + gGz.^2);

% Now for a magnetic dirac delta matrix

Mag = zeros(2*nx+1,2*ny+1, 2*nz+1);

% Add in the dirac delta function somewhere

Mag(floor(length(Mag(:,1))/2-1),floor(length(Mag(:,1))/2),floor(length(Mag(:,1))/2)) = 1;
Mag(floor(length(Mag(:,1))/2+1),floor(length(Mag(:,1))/2),floor(length(Mag(:,1))/2)) = -1;

% Get on with the plotting of these so we can see what we're dealing with
% here.
%------------------------------------------------------------------------

Msat = 10^6; % [Am^2]
volume = (rangex/nx)*(rangey/ny)*(rangez/nz); % [m^3]
scaling = Msat*volume; % To be gactored into the Dirac delta functions 

Mag = Mag.*scaling;

MagFFT = fftn(Mag,[mx,my,mz]);

GreensFFT = fftn(GG,[mx,my,mz]);

both = MagFFT.*GreensFFT;

inverse =  ifftn(both,[mx,my,mz]);

in = [nx+1,ny+1,nz+1];
en = [(3*nx+1),(3*ny+1),(3*nz+1)];

real_part = inverse(in(1):en(1),in(2):en(2),in(3):en(3));

CXY = convn(Mag,GG, 'same');

n=1;

% now to start to think about plotting the two component vector.
% need each individual component. 

FFTgGx = fftn(gGx,[mx,my,mz]);
FFTgGy = fftn(gGy,[mx,my,mz]);
FFTgGz = fftn(gGz,[mx,my,mz]);

fftX = MagFFT.*FFTgGx;
fftY = MagFFT.*FFTgGy;
fftZ = MagFFT.*FFTgGz;

iFFTx = ifftn(fftX,[mx,my,mz]);
iFFTy = ifftn(fftY,[mx,my,mz]);
iFFTz = ifftn(fftZ,[mx,my,mz]);

Xcomp = iFFTx(in(1):en(1),in(2):en(2),in(3):en(3)); 
Ycomp = iFFTy(in(1):en(1),in(2):en(2),in(3):en(3)); 
Zcomp = iFFTz(in(1):en(1),in(2):en(2),in(3):en(3));
%%
figure(n)
clf

xslice = 0;
yslice = [];
zslice = [];

slice(X,Y,Z, real_part, xslice,yslice,zslice)
colorbar
caxis([-10^-10,10^-10])
hold on


%%

    filterx = abs(Xcomp)<= 10^-9;
    filtery = abs(Ycomp)<= 10^-9;
    filterz = abs(Zcomp)<= 10^-9;
    filter = filterx.*filtery.*filterz;

    Xcompfil = Xcomp.*filter;
    Ycompfil = Ycomp.*filter;
    Zcompfil = Zcomp.*filter;
    
    totfil = sqrt(Xcompfil.^2 + Ycompfil.^2 + Zcompfil.^2);
    
    
quiver3(X,Y,Z,Xcompfil,Ycompfil,Zcompfil);

%%
n = n+1;

realdat2D = real_part(:,:,round(length(X)/2)+1);
%realdat2D = reshape(realdat2D, size(X,2), size(X,3));
quivdatX = Xcompfil(:,:,round(length(X)/2)+1);
%quivdatX = reshape(quivdatX, size(X,2), size(X,3));
quivdatZ = Zcompfil(:,:,round(length(Z)/2)+1);
%quivdatZ = reshape(quivdatZ, size(Z,2), size(Z,3));
quivdatY = Ycompfil(:,:,round(length(Z)/2)+1);

quivx = X(:,:,round(length(X/2)));
quivz = Z(:,:,round(length(Z/2)));
quivy = Y(:,:,round(length(Z/2)));

figure(n)
clf
imagesc(linex,liney,realdat2D)
colorbar
caxis([-10^-10,10^-10])
hold on 
%quiver (quivx,quivdatX)
%quiver (quivz,quivdatZ)
%quiver (quivy,quivdatY)
quiver(quivx,quivy,quivdatX,quivdatY)
title ('2D plane of 3D output')

%%

% A tester to see the profile of the field.

n = n+1;

middlefield = zeros(1,size(real_part,3));

for o = 1:size(real_part,3)  
    go = real_part(:,:,o);
    go = reshape(go,size(real_part,1),size(real_part,2));
    
    middlefield(o) = go(round((size(go,1))/2),round((size(go,1))/2));
    
end 

figure(n)
plot(linex,middlefield)
axis([-10^-3,10^-3,-2*10^-10,0])
ylabel 'Field (T)'
xlabel ' Distance from the particle (m)'
title 'Field profile at a distance'