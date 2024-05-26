%--------------------------------------------------------------------------
% nmr_plot  - plots signal and fit from nmr structure
%
% Version: 2.0
% Author: Anton Potocnik, F5, IJS
% Date:   27.10.2008 - 25.01.2009
% Arguments:
%       nmr_plot(nmr,vect,[space],[with_fit],[normalize])
%--------------------------------------------------------------------------

function nmr_plot(nmr, v, varargin)
space = 0;
wf = 0;
norm_ = 1;

if nargin>2
    space = varargin{1};
end
if nargin>3
    wf = varargin{2};
end
if nargin>4
    norm_ = varargin{3};
end


figure;
% 
% axes('LineWidth',1.5,...
%     'FontSize',14,...
%     'FontName','arial');


k=0;
for i=v
    T = nmr.data{i}.T;
    X = nmr.data{i}.X;
    Y = nmr.data{i}.Y;
            
    if norm_ > 0
        mmin = min(Y);
        mmax = max(Y);
        norm = 1/(mmax-mmin);
        Y = (Y-Y(end))*norm;
    else
        norm = 1;
    end
    
    if wf==true  % With fit?
        f1 = nmr.sim.fits{i}.f;
        fY = f1(H);
        if norm_ > 0
            fY = (fY-fY(1))*norm;
        end
        plot(H,k*space+Y, H,k*space + fY,'LineWidth',1.0);
    else
        if isfield(nmr.data{i},'X')
            X = nmr.data{i}.X;
            if norm_ > 0
                X = (X-X(1))*norm;
            end
            plot(H,k*space+Y,H,k*space+X,'c','LineWidth',1.0);
        else
            plot(H,k*space+Y,'LineWidth',1.5);
        end
    end
    hold on
%     text(min(H),k*space+Y(1),sprintf('%3.0fK %3.1fGHz',nmr.temp(i), nmr.freq(i)),'VerticalAlignment','bottom','FontSize',12,'FontName','Arial');
    text(max(H),k*space+Y(end),sprintf('%3.0fK',nmr.temp(i)),'VerticalAlignment','bottom','FontSize',12,'FontName','Arial');
    k=k+1;
end
               
% specter = [nmr.data{v}.H nmr.data{v}.Y nmr.data{v}.H f1(nmr.data{v}.H)];
xlabel('H (G)','FontSize',14,'FontName','Arial')
ylabel('dP/dH (arb. units)','FontSize',14,'FontName','Arial')
title([nmr.material '  ' nmr.date],'FontSize',18,'FontName','Arial')
hold off
    
    