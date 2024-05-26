%--------------------------------------------------------------------------
% sevenmrload  - loads data files from our Sevenmr NMR spectrometers
%
% Author: Anton Potocnik, F5, IJS
% Date:   24.10.2009
% Arguments:
%       [T, Y, exp] = sevenmrload(filename)
%--------------------------------------------------------------------------

function [T, X, Y, exp] = sevenmrload(filename) 

fid = fopen(filename);
Cont=1;
while Cont,
    S = fgetl(fid); %disp(S);
	[Ukaz, VredS]	= strtok(S,'=');

    VredS2  = VredS(2:length(VredS));
    if length(VredS) > 1,
       Vred = str2double( VredS2 );
    else
       Vred = 0;
    end
        %disp([Ukaz '->' VredS]);
    if strcmp(Ukaz,'SF'),
        exp.Freq = Vred/1e6; %Ceprav je to v bistvu hitrost digitizerja
    elseif strcmp(Ukaz,'FR'),
        exp.Freq = Vred; %Ceprav je to v bistvu hitrost digitizerja
    elseif strcmp(Ukaz,'DW'),
        exp.DW = Vred;
    elseif strcmp(Ukaz,'NSC'),
        exp.NS = Vred;
    elseif strcmp(Ukaz,'D4'),
        exp.D4 = Vred;
    elseif strcmp(Ukaz,'D0'),
        exp.D0 = Vred;
    elseif strcmp(Ukaz,'D1'),
        exp.D1 = Vred;
    elseif strcmp(Ukaz,'TAU'),
        exp.TAU	= Vred;
    elseif  strcmp(Ukaz,'TD'),
        exp.TD	= Vred;      
    elseif  strcmp(Ukaz,'_ITC_R0'),
        exp.Temperature	= Vred;
    elseif  strcmp(Ukaz,'TIMESTA'),
        exp.TIMESTA	= VredS2;
    elseif  strcmp(Ukaz,'DATESTA'),
        exp.DATESTA	= VredS2;
    elseif  strcmp(Ukaz,'TIMEEND'),
        exp.TIMEEND	= VredS2;
    elseif  strcmp(Ukaz,'DATEEND'),
        exp.DATEEND	= VredS2;
    elseif  strcmp(Ukaz,'PPFILE'),
        exp.COMMENT	= VredS2;
    elseif  strcmp(Ukaz,'[DATA]'),
        Cont = 0;
        if exp.TD>0
            [NMR_FIDxy,count] = fscanf(fid, '%f', [2, exp.TD]);
            X = NMR_FIDxy(1,:);
            Y = NMR_FIDxy(2,:);
            
%             if count ~= 2*exp.TD 
%                 error('NMR .DAT reading error!!!');
%             end
            if count < 2*exp.TD 
                Y = [Y zeros(1,(exp.TD-numel(Y)))];
                X = [X zeros(1,(exp.TD-numel(X)))];
            end
            if count > 2*exp.TD 
                Y = Y(1:exp.TD);
                X = X(1:exp.TD);
            end
            
            if numel(X)==exp.TD %ker je zadnja pika vcasih pokvarjena, jo "zbrisem"
                X(exp.TD)= X(exp.TD-1);
                Y(exp.TD)= Y(exp.TD-1);
            end
        end
    else
        %disp([Ukaz, VredS]);
    end
		
end
fclose(fid);

[m n]=size(X);
if n>m
    X=X';
    Y=Y';
end
% 
%     desc = textscan(fid, '%s',1);                  % [PARAMETERS]
%     
%     % Check for correct file format
%     if strcmp(desc{1},'[PARAMETERS]') == 0
%         fclose(fid);
%         error(['Wrong file format: ' filename]);
%     end
%     
%     desc = textscan(fid, '%s',44,'delimiter','=');                  % read 
%     for i=1:2:44
%        eval(['exp.par.' desc{1}{i} '=' desc{1}{i+1} ';']) 
%     end
%     
%     %desc = textscan(fid, '%s',1);                  % [ADDITIONAL]
%      % search for [ADDITIONAL]
%     while(strcmp(desc{1},'[ADDITIONAL]') == 0)
%         desc = textscan(fid, '%s',1);
%     end
%     
%     desc = textscan(fid, '%s',20,'delimiter','=');                  % read
%     for i=1:2:20
%        eval(['exp.add.' desc{1}{i} '=''' desc{1}{i+1} ''';']) 
%     end
%     
%     %desc = textscan(fid, '%s',1);                  % [VARIABLES]
%      % search for [VARIABLES]
%     while(strcmp(desc{1},'[VARIABLES]') == 0)
%         desc = textscan(fid, '%s',1);
%     end
%     desc = textscan(fid, '%s',34,'delimiter','=');                  % read
%     for i=1:2:34
%        eval(['exp.var.v' desc{1}{i} '=''' desc{1}{i+1} ''';']) 
%     end
% 
%     % search for [DATA]
%     while(strcmp(desc{1},'[DATA]') == 0)
%         desc = textscan(fid, '%s',1);
%     end
% 
%     C = textscan(fid, '%s %s');

% fclose(fid);

    
% convert to double matrix
% A = cell2mat(cellfun(@(x) str2double(x),C,'UniformOutput',false));
% X = A(:,1);
% Y = A(:,2);

% replace "," with "." and convert to double for exp structure


%.Temperature, .Frequency, .TD,
%                        .DW, .TAU, .NS, .NSC, .GAIN, 
%                        .D0, .D1, .D2,...
% exp.Temperature = str2double(exp.var.v_ITC_R0);
% exp.Freq = str2double(exp.var.v_FR_SYNT);
% exp.TD = exp.par.TD;
% exp.DW = exp.par.DW;
% exp.TAU = exp.par.TAU;
% exp.NS = exp.par.NS;
% exp.NSC = exp.par.NSC;
% exp.GAIN = exp.par.GAIN;

T = (1:exp.TD)';%*exp.DW;
% % Correct amplitude
% Y = Y/exp.GAIN;
    
    
    