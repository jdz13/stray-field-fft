%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the components from the FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partno = 1;

figure(8)
clf
subplot(2,2,1)
slice(X,Y,Z,Xcomp, [],[],0)
caxis([-10^-10,10^-10])
polarmap
colorbar
title 'X component'

subplot(2,2,2)
slice(X,Y,Z,Ycomp, 0,[],0)
caxis([-10^-10,10^-10])
polarmap
colorbar
title 'Y component'

subplot(2,2,3)
slice(X,Y,Z,Zcomp, 0,[],[])
caxis([-10^-10,10^-10])
polarmap
colorbar
title 'Z component'

subplot(2,2,4)
slice(X,Y,Z,sqrt(Xcomp.^2+Ycomp.^2+Zcomp.^2), 0,0,[])
caxis([-10^-10,10^-10])
polarmap
colorbar
title 'tot B FFT'