%--------------------------------------------------------------------------
% nmr_sort  - sorts specters in nmr structure
%
% Version: 1.1
% Author: Anton Potocnik, F5, IJS
% Date:   04.11.2009
%       
% Input:    
%       nmr         
%       nmr.sort    temp,freq,dates
%       nmr.sortDir descend,ascend
%--------------------------------------------------------------------------


if ~isfield(nmr,'sortDir')
    nmr.sortDir = 'descend';  % 'ascend'
end

% SORT in ascend order
[v ix] = sort(nmr.(nmr.sort));


if strcmp(nmr.sortDir,'descend')    % problems with cell arrays
    ix = ix(numel(ix):-1:1);
end

nmr.data = reshape(nmr.data(ix),[],1);
nmr.fft = reshape(nmr.fft(ix),[],1);
nmr.temp = reshape(nmr.temp(ix),[],1);
nmr.freq = reshape(nmr.freq(ix),[],1);
nmr.dates = reshape(nmr.dates(ix),[],1);

% if isfield(nmr,'sim')
%     nmr.sim = reshape(nmr.sim(ix),[],1);
% end

clear v ix
    