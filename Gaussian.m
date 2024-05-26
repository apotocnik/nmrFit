%==========================================================================
% Gaussian function
% Anton Potoènik @ IJS, F5,  oct. 2009 - jun. 2012 
%
% LB, Kii in kHz frequency units!!!
% DW in micro seconds
% powder averagn method: zcw
%==========================================================================


function [spc,f] = Gaussian(fmax,TD,A,LB,fc)
%%
Y = zeros(TD,1);
f = linspace(-fmax,fmax,TD)';
   
    Y = Y + 1/sqrt(2/pi)/LB.*exp(-(f-fc).^2/2/LB^2); % prefactor: 1/sqrt(2/pi)

%spc = abs(real(Y));
spc = Y/max(Y);  %Normalize
spc = A*spc;
