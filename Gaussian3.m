%==========================================================================
% Gaussian function
% Anton Potoènik @ IJS, F5,  oct. 2009 - jun. 2012 
%
% LB, Kii in kHz frequency units!!!
% DW in micro seconds
% powder averagn method: zcw
%==========================================================================


function [spc,f] = Gaussian3(fmax,TD,A1,LB1,fc1,A2,LB2,fc2,A3,LB3,fc3)
%%
Y = zeros(TD,1);
f = linspace(-fmax,fmax,TD)';
   
    Y = Y + A1/sqrt(2*pi)/LB1.*exp(-(f-fc1).^2/2/LB1^2); % prefactor: 1/sqrt(2/pi)
    Y = Y + A2/sqrt(2*pi)/LB2.*exp(-(f-fc2).^2/2/LB2^2); % prefactor: 1/sqrt(2/pi)
    Y = Y + A3/sqrt(2*pi)/LB3.*exp(-(f-fc3).^2/2/LB3^2); % prefactor: 1/sqrt(2/pi)

%spc = abs(real(Y));
% spc = Y/max(Y);  %Normalize
% spc = A*spc;
spc = Y;
