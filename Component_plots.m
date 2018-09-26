%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the components from the FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partno = 1;

figure(8)
clf
subplot(2,2,1)
slice(X,Y,Z,Xcomp./mu0, [],[],0)
caxis([-10^-4,10^-4])
polarmap
colorbar
title 'X component'
axis([-10^-3,10^-3,-10^-3,10^-3,-10^-3,10^-3])

subplot(2,2,2)
slice(X,Y,Z,Ycomp./mu0, 0,[],0)
caxis([-10^-4,10^-4])
polarmap
colorbar
title 'Y component'
axis([-10^-3,10^-3,-10^-3,10^-3,-10^-3,10^-3])

subplot(2,2,3)
slice(X,Y,Z,Zcomp./mu0, 0,[],[])
caxis([-10^-4,10^-4])
polarmap
colorbar
title 'Z component'
axis([-10^-3,10^-3,-10^-3,10^-3,-10^-3,10^-3])

subplot(2,2,4)
slice(X,Y,Z,sqrt(Xcomp.^2+Ycomp.^2+Zcomp.^2)./mu0, 0,0,[])
caxis([-10^-4,10^-4])
polarmap
colorbar
title 'tot B FFT'
axis([-10^-3,10^-3,-10^-3,10^-3,-10^-3,10^-3])