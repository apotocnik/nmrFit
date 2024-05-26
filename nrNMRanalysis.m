%--------------------------------------------------------------------------
% nrNMRanalysis  -  Numerical analysis of NMR spectra
%
% Author: Anton Potocnik, F5, IJS
% Date:   31.10.2009 - 19.11.2012
% Arguments:
%       [M0, M1, M2] = nranalysis(f, Y, [bline_corr], [int_cutoff])
% Input:
%        f, Y        - FFT specter only REAL part!!! (kHZ)
%               new: - if Y is imag, analyse with Error Bar
%        bline_corr  - base line linear correction 
%                      [freq_span_left freq_span_right];
%                      the same units as freq
%        int_cutoff  - analize only above int_cutoff
% Output:
%        M0  - Intensity (area under specter integral)
%        M1  - center of specter
%        M2  - width of specter (= sqrt(M2))
%--------------------------------------------------------------------------

function [M0 M1 M2] = nrNMRanalysis(f, Y, varargin)
bline_corr = 0;
int_cutoff = 0; 
n = numel(varargin);
if n>0
    bline_corr = varargin{1};
end
if n>1
    int_cutoff = varargin{2};
end

if ~isreal(Y)
    err = imag(Y);
end
Y = real(Y);

%% Integral Base Line Linear Correction
% Convert freq values to index
if numel(bline_corr) == 2
   [C left_idx]  = min(abs(f-(f(1)+bline_corr(1))));
   [C right_idx] = min(abs(f-(f(end)-bline_corr(2))));
else
   [C left_idx]  = min(abs(f-(f(1)+bline_corr)));
   [C right_idx] = min(abs(f-(f(end)-bline_corr)));
end

if left_idx > 1 || right_idx < numel(f)     % Do nothing if bline_corr==0
    % Exctract data for linear fit
    x = f([1:left_idx right_idx:end]);
    z = Y([1:left_idx right_idx:end]);
    % Linear fit
    p = polyfit(x,z,1);

    % Correct base line
    Y = Y - polyval(p,f);
end

%% Integrate with cutoff
% if cutoff is empty dont cutoff
if isempty(int_cutoff)
    int_cutoff = -Inf;
end
y = Y;

if int_cutoff > -Inf
    y = y - int_cutoff;
    y(y<0) = 0;
end

%% M0, M1, M2 analysis

if exist('err','var')
    err(isnan(err)) = 1e10;
%     w = 1./err;
%     M0 = trapz(f,y.*w);
%     M1 = 1/M0*trapz(f,y.*f.*w);
%     M2 = sqrt(1/M0*trapz(f,y.*(f-M1).*(f-M1).*w));
    y = y - err;
    y(y<0) = 0;
    M0 = trapz(f,y);
    M1 = 1/M0*trapz(f,y.*f);
    M2 = sqrt(1/M0*trapz(f,y.*(f-M1).*(f-M1)));
else
%     y=abs(y);
    M0 = trapz(f,y);
    M1 = 1/M0*trapz(f,y.*f);
    M2 = sqrt(1/M0*trapz(f,y.*(f-M1).*(f-M1)));
end






    
    