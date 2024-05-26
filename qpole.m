%==========================================================================
% Wrapper for simulating quadrupole NMR spectra
% Anton Potoènik @ IJS, F5,  19.02.2014 
%
% For Spin < 0, only central transition is computed
%
%  Kiso = (Kxx + Kyy + Kzz)/3    ...  Slichter book
%  Kani = (2*Kzz - Kxx - Kyy)/3
%  Kasy = Kxx - Kyy
%
% ! K tensor orientation is aligned with quadrupole tensor orientation!
%==========================================================================


function [spc,x] = qpole(N,fmax,NOP,Spin,A,LB,Kiso,Kani,Kasy,niQ,eta,DniQ,offset)
%%

[fs Ws] = qpoleHist2(N, Spin, Kiso, Kani, Kasy, niQ, eta, DniQ);

x = linspace(-fmax,fmax,NOP+2)';
spc = zeros(size(x));

for i=1:size(Ws)
   n = hist(fs(:,i), x)';
   spc = spc + Ws(i)*n;
end
x(1)=[]; x(end)=[];
spc(1) = []; spc(end) = [];

% Ispc = trapz(x,spc);

DW = 1/(x(2)-x(1))/NOP;
T = DW*(0:NOP-1)';

FID = fft(spc);
LBfun = exp(-T*LB)+exp(-(max(T)+DW-T)*LB);
FID = FID.*LBfun/2;

SPC = fft(FID)/NOP;
SPC = [SPC(end); SPC(end:-1:2)];
% Ispcp = trapz(x,real(SPC));
% SPC = Ispc/Ispcp*SPC;

% figure(1)
% plot(T,real(FID),'b',T,imag(FID),'r');
% figure(2)
% plot(x,real(SPC),'r',x,spc,'k');

spc = real(SPC)/max(real(SPC))*A + offset;

