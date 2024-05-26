%--------------------------------------------------------------------------
%% NMR SIGNAL ANALYSIS Load data
% Version 2.1 - Loads almost everything
% from EPR_load
% Author: Anton Potocnik, F5, IJS
% Date:   24.10.2009 - 24.10.2009
% Input:  nmr.path 
%         [nmr.sort]    % 'temp' 'date'
% Output  nmr.data{}.exp .fname, .Temperature, .Frequency, .TD,
%                        .DW, .TAU, .NS, .NSC, .GAIN, 
%                        .D0, .D1, .D2,...
%                   .T   (Time [secunds])
%                   .Y   (Signal [Volts])
%            .temp()
%            .dates()
%            .freq()
%            .N         % numel(nmr.data)
%            .fft{}
%--------------------------------------------------------------------------

%% Get file list and its information
files = dir(nmr.path);
fpath = fileparts(nmr.path);

if isempty(files)
    error('Files not found!')
end

%% Get previous data
ind = 1;   % data counter, different than i if data allready exist
if isfield(nmr,'data')      % if data allready exists add new data
    ind = numel(nmr.data)+1;
    data = nmr.data;
    temp = nmr.temp;
    date = nmr.dates;
    freq = nmr.freq;
end

%% Load files
disp(' ');
disp('##########################################');
disp('Loading data from:');
disp(['  ' fpath]);

for i=1:numel(files)
    if files(i).isdir == true  % Skip Folders
        continue;
    end
    fname = fullfile(fpath,files(i).name);
    
    % Recognize file format
    desc = textread(fname,'%s',1,'bufsize',8190); 
    switch desc{1}
        case '[PARAMETERS]'              
            [T X Y exp] = sevenmrload(fname);      
        otherwise
            disp(['File not recognised!!' files(i).name]);
            return;
    end
    
    tmp = pwd;
    cd(fpath);
    data{ind}.fname = fullfile(pwd,files(i).name);
    cd(tmp);
    clear tmp
    data{ind}.date = files(i).date;
    data{ind}.T = T;
    data{ind}.signal = complex(X,-Y); 
    if ~isfield(exp,'Temperature')
        exp.Temperature = 295;     % Room temperature 22°C
    end
    data{ind}.exp = exp;
    
    temp(ind) = exp.Temperature;
    freq(ind) = exp.Freq;
    
    if isfield(exp,'DATEEND')
        date{ind} = [exp.DATEEND ' ' exp.TIMEEND];
    else
        date{ind} = files(i).date;
    end
    
    sprintf('%d\t%s\tT=%3.2fK\t f=%3.2fMHz',ind,files(i).name,temp(ind),freq(ind));
    ind = ind+1;
end

%% SORT and write to nmr structure
% if isfield(nmr,'sort')
%     switch nmr.sort
%         case 'dates'
%             [v ix] = sort(dates);
%         case 'freq'
%             [v ix] = sort(freq);
%         otherwise
%             [v ix] = sort(temp,'ascend');
%     end
% else 
%     [v ix] = sort(temp,'descend');
% end
% nmr.data = reshape(data(ix),[],1);
% nmr.temp = reshape(temp(ix),[],1);
% nmr.freq = reshape(freq(ix),[],1);
% nmr.dates = reshape(date(ix),[],1);
nmr.data = reshape(data,[],1);
nmr.temp = reshape(temp,[],1);
nmr.freq = reshape(freq,[],1);
nmr.dates = reshape(date,[],1);

nmr.N = numel(nmr.data);

clear temp freq date X Y T i ind files fpath v ix data desc exp fname fft m PHASE DE SHL LB spc f

