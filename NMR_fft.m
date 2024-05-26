%--------------------------------------------------------------------------
% NMR SIGNAL Fast Fourier Transform
% Version 1.0
%
% Author: Anton Potocnik, F5, IJS
% Date:   25.10.2009 - 25.10.2009
% Input:  nmr.fft.shl      ... shift left
%         nmr.fft.phs      ... shift phase
%         nmr.fft.chosen   ... vector of chosen signals
%         nmr.fft.plot     ... plot=1 noplot=0
%
% Output: nmr.fft.spc      ... specter
%         nmr.fft.phase    ... phase
%--------------------------------------------------------------------------
% TODO
% 
%--------------------------------------------------------------------------


disp(' ');
disp('##########################################');
disp('Fast Fourier Transform');


for i = nmr.fft.chosen

    T = nmr.data{i}.T;
    X = nmr.data{i}.X;
    Y = nmr.data{i}.Y;
    
    [Z A w xc] = nranalysis(H, Y, epr.nra.bline_corr, epr.nra.cutoff, epr.nra.w_method, epr.nra.xc_method);
    
    if nmr.fft.plot == 1
        plot(H,Z)
        legend(['T = ' num2str(epr.temp(i)) 'K']);
        %axis([1500 5000 0 20])
        Mov(i) = getframe;
    end
    
    
    disp(sprintf('T=%3.0fK\tA=%5.3e\tw=%3.3f\txc=%3.3f\t', [nmr.temp(i),A,w,xc]));
end

%% Save specter

nmr.fft.spc = spc;
nmr.fft.phase = phase;


% %% Plot
% if epr.nra.plot==1
%     plot_results(epr.nra.results_g,epr)
% end


%% Clear
clear i T X Y spc phase Mov
