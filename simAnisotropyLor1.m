% Normal powder average for anisotropic chemical shifts plus lorentzian
% LB, Kii in kHz frequency units!!!
% DW in micro seconds

function [spc,f] = simAnisotropyLor1(N,DW,TD,A,LB,Kxx,Kyy,Kzz,D,Klor,W)
%%
W = W*1e3;
Klor = Klor*1e3;

[spc f] = nmrpowderF(N,DW,TD,A,LB,Kxx,Kyy,Kzz);

Ysurface = trapz(f,spc);
Y = spc / Ysurface; % Normalized powder spectrum

% Y = zeros(size(f));
Y = Y + D*2/pi*W./(4*(f-Klor).^2+W^2);

spc = A*Y/max(Y);

