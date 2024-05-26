%--------------------------------------------------------------------------
% nmr_delete  - delete specters in nmr structure
%
% Version: 2.0
% Author: Anton Potocnik, F5, IJS
% Date:   25.10.2009 - 5.1.2013
%       
% Arguments nmr = nmr_delete(epr,v)
% Input:    
%       nmr         
%       v           vector of indeces to delete
%--------------------------------------------------------------------------

function nmr = nmr_delete(nmr,v)

    try nmr.data(v)=[]; catch e; end;
    try nmr.temp(v)=[]; catch e; end;
    try nmr.freq(v)=[]; catch e; end;
    try nmr.dates(v)=[]; catch e; end;
    try nmr.fft(v)=[]; catch e; end;
    try nmr.sim(v)=[]; catch e; end;
    try nmr.nra(v)=[]; catch e; end;

    nmr.N = numel(nmr.data);
    
    