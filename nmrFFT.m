%--------------------------------------------------------------------------
% nmrFFT  - fourier transform of NMR fid signal with auto phase
%
% Author: Anton Potocnik, F5, IJS
% Date:   04.11.2009-06.06.2014
% Arguments:
%       [spc f PHS] = nmrFFT(fid,SHL,PHS,LB,DE,TD,DW,avrPoints) 
% INPUT:
%       fid ... fid signal (complex)
%       SHL ... shift left (points), auto if sqrt(-1), fit squ. if -sqrt(-1)
%       PHS ... phase (rad), auto if sqrt(-1)
%       LB  ... Real part:Lorentzian, Imag part: Gaussian broadening (kHz)
%       DE  ... double echo (0 fill with zeros|1 fill with left of echo)
%       TD  ... time domain (2^N!!!)
%       DW  ... dwell time (s)
%       avrPoints ...  number of averaging points
% OUTPUT:
%       spc ... specter (complex)
%       f   ... frequency axis (Hz)
%       PHS ... phase (from auto phase)
%       SHL ... shift left (from auto SHL)
%       fid ... changed fid
%--------------------------------------------------------------------------

function [spc f PHS SHL fid] = nmrFFT(fid,SHL,PHS,LB,DE,TD,DW,avrPoints) 

M = round((avrPoints-1)/2);


% =========================================================================
% --- SHL
% =========================================================================

if ~isreal(SHL)  % auto SHL
    [bla ind] = max(abs(fid)); 
    
    if imag(SHL) < 0  % Fit square function to the avrPoints to determine SHL
        if avrPoints > 2
            rng = (ind-M:ind+M)';
            rng(rng<1) = 1; 
            rng(rng>TD) = TD;
            sig = abs(fid(rng));
            p = polyfit(rng,sig,2);
            if p(1) < 0 && p(2) > 0
               ind = round(-p(2)/2/p(1));
            end
        else
           msgbox('AvrPoints < 2! Cannot fit square for SHL.'); 
           return
        end
    end
    SHL = ind;
end


% =========================================================================
% --- PHASE
% =========================================================================

if ~isreal(PHS) % Auto phase correction
    rng = (SHL-M:SHL+M)';
    rng(rng<1) = 1; 
    rng(rng>TD) = TD;
    
    PHS = -angle(mean(fid(rng))); 
end
fid = fid*exp(sqrt(-1)*PHS);  % Phase correction of the fid


% =========================================================================
% --- Double Echo
% =========================================================================

if DE == 1 % Fill with fid
    fid = [fid(SHL:end); fid(1:SHL-1)]; 
else
    fid = [fid(SHL:end); zeros(SHL-1,1)];
end


% =========================================================================
% --- Line Broadening
% =========================================================================

T = DW*(1:TD)'; % seconds
LBexp = real(LB);
LBgau = imag(LB);

if DE
    LBfunLor = exp(-T*LBexp) + exp((T-max(T)-DW)*LBexp);
    LBfunGau = exp(-(T*LBgau).*(T*LBgau)) + exp((T-max(T)-DW)*LBgau.*(T-max(T)-DW)*LBgau);
    LBfunLor = LBfunLor/max(LBfunLor);
    LBfunGau = LBfunGau/max(LBfunGau);
else
    LBfunLor = exp(-T*LBexp);
    LBfunGau = exp(-(T*LBgau).*(T*LBgau));
end

fid = fid.*LBfunLor.*LBfunGau;


% =========================================================================
% --- Fourier Transform 
% =========================================================================

if DE
    % It is already periodic!!!
else
    fid(1) = (fid(1)+fid(end))/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
spc = fftshift(fft(fid,TD));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs = 1/DW;
f = Fs*(-TD/2:TD/2-1)'/TD; %Hz



% %======================================================================
% % Check SPC peak phases and corrects if set
% %======================================================================
% 
% if SPCphase > 0
%     
%     for i=1:indN
%         [ind,cv] = searchclosest(Fall{i},PH_v);
%         phi0 = angle(SPCall{i}(ind));
%         phiall_spc(i,:) = [xPar(i),phi0];
%     end
% 
%     switch autoPH_SPC
%         case 1
%             PHImean = mean(phiall_spc(intersect(autoCorrInd,IND),2));
%             phi_spc = PHImean*ones(indN,1);
% 
%         case 2
%             phi_spc = phiall_spc(:,2);
%         
%         otherwise
%             error('Wrong autoPH_SPC option!')
%     end
% 
%     for i=1:indN
%         SPCall{i} = SPCall{i}*exp(-complex(0,1)*phi_spc(i));
%         phiall(i,3) = phiall(i,3) + phi_spc(i);
%     end
% end



