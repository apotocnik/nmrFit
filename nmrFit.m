function varargout = nmrFit(varargin)
% nmrFIT M-file for nmrFit.fig
%      nmrFIT, by itself, creates a new nmrFIT or raises the existing
%      singleton*.
%
%      H = nmrFIT returns the handle to a new nmrFIT or the handle to
%      the existing singleton*.
%
%      nmrFIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in nmrFIT.M with the given input arguments.
%
%      nmrFIT('Property','Value',...) creates a new nmrFIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nmrFit_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nmrFit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nmrFit

% Last Modified by GUIDE v2.5 08-Jun-2014 22:39:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nmrFit_OpeningFcn, ...
                   'gui_OutputFcn',  @nmrFit_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before nmrFit is made visible.
function nmrFit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nmrFit (see VARARGIN)

% Choose default command line output for nmrFit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using nmrFit.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5,1));
end

if ~evalin('base','exist(''nmr'',''var'')')
    assignin('base','nmr',[]);
end

% Initialize popSumLib
sims = nmrsim_lib('');  % Get simmulations name cell array
set(handles.popSimLib, 'String',sims);

% UIWAIT makes nmrFit wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Set pointer values to guidata
pointers.X1 = 0;
pointers.X2 = 0;
pointers.Y1 = 0;
pointers.Y2 = 0;
set(handles.figure1,'UserData',pointers);   % Set pointers to figure1


% --- Outputs from this function are returned to the command line.
function varargout = nmrFit_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)

    printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)




% --- Executes on selection change in lbFiles.
function lbFiles_Callback(hObject, eventdata, handles)
    axes(handles.axes1);
    cla;

    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');

%     if numel(nmr.fft) < idx
%         butFFT_Callback(hObject, eventdata, handles);
%     end
    
    nmrplot(nmr,idx,handles);
    nmr_update(handles);




% --- Executes during object creation, after setting all properties.
function lbFiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% reloadFileListBox(handles);



% User-function: Reload entries in FileListBox
function reloadFileListBox(handles)
if evalin('base','exist(''nmr'',''var'')') == 1
    nmr = evalin('base','nmr');
    if isfield(nmr,'data')
        for i=1:numel(nmr.data)
            [pathstr, name, ext] = fileparts(nmr.data{i}.fname);
            shownames{i} = name;
        end
        set(handles.lbFiles, 'String', shownames);
    else
        set(handles.lbFiles, 'String', []);
    end
else
    set(handles.lbFiles, 'String', []);
end




% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)



% User function: plots data with fit
function nmrplot(nmr,idx,handles)
% nmr       ... nmr structure
% idx       ... specter index (could be a vector)
% handles   ... handles for check Fit and FFT/Spectrum

    chk_Sim = get(handles.chkPlotSim, 'Value');
    chk_Err = get(handles.chkErrBar, 'Value');
    
%     fft_spc = get(handles.chkFFTSpc,'Value');
    tbutSPC = get(handles.tbutSPC,'Value');
    tbutFID = get(handles.tbutFID,'Value');
    tbutmFID = get(handles.tbutmFID,'Value');
    

    %axes(handles.axes1); % It is making it very slow!!!
    %cla;    % Clear graph
    color = 'brgmcy';
    j=0;
    M={};

    for i=idx
        
        % Check datasets --------------------------------------------------
        isFID = false; 
        ismFID = false;
        if isElement(nmr,'data',i) 
            if numel(nmr.data{i}.signal) > 3 % Check for measured FID
                isFID = true;
            end
            if isfield(nmr.data{i},'signal_new') % Check for modified FID
                ismFID = true;
            end
        end
        isSPC = isElement(nmr,'fft',i); % Check for calculated SPC
        isSIM = isElement(nmr,'sim',i); % Check for calculated SIMulation
        
        % Check what datasets can be displayed ----------------------------
        if tbutmFID == 1 && ismFID == false
            disp(['No mFID available: ' num2str(i) '!'])
            tbutmFID = 0;
            tbutFID = 1;
        end
        
        if tbutFID == 1 && isFID == false
            disp(['No FID available: ' num2str(i) '!'])
            tbutFID = 0;
            tbutSPC = 1; % Probably only spectrum
        end
        
        if tbutSPC == 1 && isSPC == false && isSIM == false
            disp(['No SPC available: ' num2str(i) '!'])
            tbutSPC = 0;
            if isFID, 
                tbutFID = 1; 
            else
                msgbox('Problem! No FID, mFID, or SPC dataset found!');
            end
        end
        
        % Start drawing ---------------------------------------------------
        % FID .............................................................
        if tbutFID  
            DW = nmr.data{i}.exp.DW;
            T = nmr.data{i}.T*DW;
            fid = nmr.data{i}.signal;

            X = real(fid);
            Y = imag(fid);
            AB = sqrt(Y.*Y+X.*X);

            plot(T,X,color(mod(j,6)+1),T,Y,color(mod(j+1,6)+1),T,AB,'k');
            hold on

            if isElement(nmr,'fft',i)
                SHL = nmr.fft{i}.shl;
                v = ylim(handles.axes1);
                plot([SHL,SHL]*DW,v,'b-.')
                
                if isfield(nmr.fft{i},'autoSHLsquare')
                    avrPoints = nmr.fft{i}.avrPoints;
                    TD = nmr.fft{i}.TD;
                    if nmr.fft{i}.autoSHLsquare == 1 && avrPoints > 1 % added at the same time
                        
                        Bo = round((avrPoints-1)/2);
                        rng = (SHL-Bo:SHL+Bo)';
                        rng(rng<1) = 1; 
                        rng(rng>TD) = TD;
                        
                        sig = abs(fid(rng));
                        p = polyfit(rng,sig,2);
%                         if p(1) < 0 && p(2) > 0
                           xx = linspace(rng(1),rng(end),100);
                           yy = polyval(p,xx);
                           plot(xx*DW,yy,'-r','LineWidth',3)
%                         end
                    end
                end
            end

            xlabel('Time (s)');
            %ylabel('Signal (V)');
            M = {'real','imag','abs'};
            xlim([min(T) max(T)]);
        end
        
        % modified FID ....................................................
        if tbutmFID 
            fid = nmr.data{i}.signal_new;
            T = (1:numel(fid))*nmr.data{i}.exp.DW;

            X = real(fid);
            Y = imag(fid);
            AB = sqrt(Y.*Y+X.*X);

            plot(T,X,color(mod(j,6)+1),T,Y,color(mod(j+1,6)+1),T,AB,'k');
            hold on

            xlabel('Time (s)');
            %ylabel('Signal (V)');
            M = {'real','imag','abs'};
            xlim([min(T) max(T)]);
        end
        
        
        % FFT spectrum ....................................................
        if tbutSPC && isSPC  
            y = real(nmr.fft{i}.spc);
            xlimit = str2num(nmr.fft{i}.xlim);
            if nmr.fft{i}.isPPM
                x = nmr.fft{i}.f;
            else
                x = nmr.fft{i}.f/1000; % display kHz
            end
            % Get only data within xlimit
%                 [C idx_low] = min(abs(x-xlimit(1)));
%                 [C idx_high] = min(abs(x-xlimit(2)));
%                 sel = idx_low:idx_high;
            sel = 1:numel(x);

            if ~isreal(nmr.fft{i}.spc) && chk_Err
                z = imag(nmr.fft{i}.spc);
                errorbar(x(sel),y(sel),z(sel),color(mod(j,6)+1)) 
            else
                plot(x(sel),y(sel),color(mod(j,6)+1))
            end
            hold on
            
            xlim(xlimit)
            
            if nmr.fft{i}.isPPM
                xlabel('Frequency (ppm)')
            else
                xlabel('Frequency (kHz)')
            end
            M{j+1} =  [ num2str(nmr.temp(i)) 'K'];
        end
        
        % SIM spectrum ....................................................
        if chk_Sim && tbutSPC && isSIM  
            plot(nmr.sim{i}.f,nmr.sim{i}.spc,'k','LineWidth',2);
            hold on
            if isSPC
                if nmr.fft{i}.isPPM
                    xlabel('Frequency (ppm)')
                else
                    xlabel('Frequency (kHz)')
                end
            else
                xlabel('Frequency')
            end
        end
        

        j=j+1;
        
    end
    hold off
    
    if get(handles.chkAutoYrange,'Value') == 0
        YL = str2double(get(handles.edtYrange,'String'));
        ylim([-YL YL]);
    else
        ylim auto
    end
    
    grid on
    legend(M);
%     set(handles.chkFFTSpc,'Value',fft_spc);


    
% User function: check if there is element in the array within structure
function answer = isElement(struct,field,idx)
    answer = false;
    
    if isfield(struct,field)
        if numel(struct.(field)) >= idx
            if ~isempty(struct.(field){idx})
                answer = true;
            end
        end
    end

        
        
    

% User function: updates all nmr text boxes
function nmr_update(handles)
    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');

    % Update file list
    set(handles.txtSell,'String',['Selected = ' num2str(idx)]);
    if numel(idx)>1
        return
    end
    

    % Update NRA 
    try
        set(handles.txtM0,'String',['M0 = ' num2str(nmr.nra{idx}.results(2))]);
        set(handles.txtM1,'String',['M1 = ' num2str(nmr.nra{idx}.results(4))]);
        set(handles.txtM2,'String',['M2 = ' num2str(nmr.nra{idx}.results(6))]);  
    catch e 
        set(handles.txtM0,'String','M0 = ');
        set(handles.txtM1,'String','M1 = ');
        set(handles.txtM2,'String','M2 = ');
    end
    
    if get(handles.chkUpdate, 'Value') == 1
        % Update NRA
        try
            set(handles.edtBLCorr,'String',mat2str(nmr.nra{idx}.bline_corr));
            set(handles.edtCutOff,'String',num2str(nmr.nra{idx}.cutoff));
            set(handles.edtNRange,'String',nmr.nra{idx}.range); 
        catch e
            set(handles.edtBLCorr,'String','[0 0]');
            set(handles.edtCutOff,'String','-Inf');
            set(handles.edtNRange,'String','[-Inf Inf]');
        end
        
        % Update simulation functions 
        if isfield(nmr,'sim')
            i=idx(end);
            if i<=numel(nmr.sim) && ~isempty(nmr.sim{i})
                sim_lib = get(handles.popSimLib, 'String');
                for j=1:numel(sim_lib)
                   if strcmp(nmr.sim{idx}.sim_name,sim_lib{j})
                       set(handles.popSimLib, 'Value', j);
                   end
                end
                set(handles.uitSimPar, 'Data', nmr.sim{i}.coefs);
                set(handles.txtChi,'String',['Chi2 = ' num2str(nmr.sim{idx}.chi2)]);
                set(handles.edtSimRange,'String',nmr.sim{idx}.range);
            end
        end % Update simulation function

        % Update FFT parameters
        if isfield(nmr,'fft') 
            i=idx(end);
            if i<=numel(nmr.fft)
                set(handles.edtSHL,'String',num2str(nmr.fft{i}.shl));
                set(handles.chkSHLauto,'Value',nmr.fft{i}.autoshl);
                set(handles.edtXLim,'String',num2str(nmr.fft{i}.xlim));  
                set(handles.edtTD,'String',num2str(nmr.fft{i}.TD));  
                set(handles.edtPhase,'String',num2str(round(nmr.fft{i}.phs/pi*180)));  
                set(handles.chkDE,'Value',nmr.fft{i}.de);  
                set(handles.edtExpBroad,'String',num2str(real(nmr.fft{i}.lb)/1000));
                set(handles.edtGauBroad,'String',num2str(imag(nmr.fft{i}.lb)/1000));
                set(handles.chkPPM,'Value',nmr.fft{i}.isPPM);
                set(handles.chkNorm,'Value',nmr.fft{i}.isNorm);
                set(handles.edtREF,'String',num2str(nmr.fft{i}.REF)); 
                
                if isfield(nmr.fft{i},'autophs')
                    set(handles.chkAutoPhase,'Value',nmr.fft{i}.autophs);
                end
                if isfield(nmr.fft{i},'autoRange')
                    set(handles.chkAutoRange,'Value',nmr.fft{i}.autoRange);
                end
                if isfield(nmr.fft{i},'turn')
                    set(handles.chk380,'Value',nmr.fft{i}.turn);
                end
                if isfield(nmr.fft{i},'autoSHLsquare')
                    set(handles.chkSHLsquare,'Value',nmr.fft{i}.autoSHLsquare);
                end
                if isfield(nmr.fft{i},'avrPoints')
                    set(handles.edtAvrPoints,'String',num2str(nmr.fft{i}.avrPoints));
                end
            end
        end
        
    end
    
    % Update data description
    
    if isfield(nmr.data{idx}.exp,'DATESTA')
    
        txt0 = ['Temp: ' sprintf('% 4i',nmr.data{idx}.exp.Temperature) 'K'];
        txt0 = [txt0 '  Date: ' nmr.data{idx}.exp.DATESTA ' ' nmr.data{idx}.exp.TIMESTA];
        txt0 = [txt0 '  FREQ: ' sprintf('%.4f',nmr.data{idx}.exp.Freq) 'MHz'];

        txt1 = ['DW: ' sprintf('%3.2f',nmr.data{idx}.exp.DW*1e6) 'us'];
        txt1 = [txt1 '  TD: ' sprintf('%4d',nmr.data{idx}.exp.TD)];
        txt1 = [txt1 '  NSC: ' sprintf('%5d',nmr.data{idx}.exp.NS)];
        txt1 = [txt1 '  TAU: ' sprintf('%3.2f',nmr.data{idx}.exp.TAU*1e6) 'us'];
        txt1 = [txt1 '  D1: ' sprintf('%3.2f',nmr.data{idx}.exp.D1*1e6) 'us'];
        txt1 = [txt1 '  D0: ' sprintf('%3.6f',nmr.data{idx}.exp.D0) 's'];
    else
        txt0 = '';
        txt1 = '';
    end
    set(handles.txtDesc,'String',{txt0,txt1});
    set(handles.txtDesc,'FontName','FixedWidth')
        

function [X Y Xind] = readPointer(handles)

    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');
    fftspc = get(handles.tbutSPC,'Value');
    [X Y] = ginput(1);
    if fftspc == 1
%         f = nmr.fft{idx}.f;
%         kx = numel(f)/(f(2)-f(1))*1000;  %kHz
%         Xind = round(X*kx);
        Xind = [];
    else
        kx = 1/nmr.data{idx}.exp.DW;
        Xind = round(X*kx);
    end

    num2clip(X);
    
    X = num2str(X);
    Y = num2str(Y);
    Xind = num2str(Xind);
   


% --- Executes on button press in butPointer1.
function butPointer1_Callback(hObject, eventdata, handles)
    
    [X Y Xind] = readPointer(handles);
    

    if isempty(Xind) 
        set(handles.butPointer1, 'String',['X: ' X '  Y: ' Y])
        disp(['Pointer1: X = ' X '  Y = ' Y]);
    else
        set(handles.butPointer1, 'String',['X: ' X '(' Xind ') Y: ' Y])
        disp(['Pointer1: X = ' X '(' Xind ')  Y = ' Y]);
    end

    pointers = get(handles.figure1,'UserData');   % Get pointers from figure1
    pointers.X1 = str2double(X);
    pointers.Y1 = str2double(Y);
    set(handles.figure1,'UserData',pointers);     % Set pointers to figure1

    dX = abs(pointers.X1-pointers.X2);
    dY = abs(pointers.Y1-pointers.Y2);
    set(handles.txtdH, 'String',['dX: ' num2str(dX,'%f') '  dY: ' num2str(dY,'%f')]);


% --- Executes on button press in butPointer2.
function butPointer2_Callback(hObject, eventdata, handles)
    
    [X Y Xind] = readPointer(handles);
    
    if isempty(Xind) 
        set(handles.butPointer2, 'String',['X: ' X '  Y: ' Y])
        disp(['Pointer2: X: ' X '  Y: ' Y]);
    else
        set(handles.butPointer2, 'String',['X: ' X '(' Xind ')  Y: ' Y])
        disp(['Pointer2: X: ' X '(' Xind ')  Y: ' Y]);
    end
    pointers = get(handles.figure1,'UserData');   % Get pointers from figure1
    pointers.X2 = str2double(X);
    pointers.Y2 = str2double(Y);
    set(handles.figure1,'UserData',pointers);     % Set pointers to figure1

    dX = abs(pointers.X1-pointers.X2);
    dY = abs(pointers.Y1-pointers.Y2);
    set(handles.txtdH, 'String',['dX: ' num2str(dX,'%f') '  dY: ' num2str(dY,'%f')]);
    


% --- Executes on button press in butNext.
function butNext_Callback(hObject, eventdata, handles)
% hObject    handle to butNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)






% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on button press in butNAnalyse.
function butNAnalyse_Callback(hObject, eventdata, handles)

    blcorr = eval(get(handles.edtBLCorr,'String'));
%     chk_AnaErr = get(handles.chkAnaErr,'Value');
%     chk_Err = get(handles.chkErrBar,'Value');
    
    
    cutoff = str2double(get(handles.edtCutOff,'String'));

    idx = get(handles.lbFiles, 'Value');
    idx = idx(1); % If more than one selected
    
    nmr = evalin('base','nmr');

    if ~isfield(nmr,'fft') || idx>numel(nmr.fft)
        disp('No specter found!');
    end

    f = nmr.fft{idx}.f;
    if nmr.fft{idx}.isPPM == 0
        f = f/1000; %kHz
    end
    
%     if ~isreal(nmr.fft{idx}.spc) && chk_Err % For Plotting
%         err = imag(nmr.fft{idx}.spc);
%         plotErr = 1;
%     else
%         plotErr = 0;
%     end
%     if ~isreal(nmr.fft{idx}.spc) && chk_AnaErr % For Analysing
%         spc = nmr.fft{idx}.spc;
%         err = imag(spc);
%     else
        spc = real(nmr.fft{idx}.spc);
%     end    
    
    % Get only data within xlimit
    xlimit = str2num(nmr.fft{idx}.xlim);
    [C idx_low] = min(abs(f-xlimit(1)));
    [C idx_high] = min(abs(f-xlimit(2)));
    sel = idx_low:idx_high;
    f = f(sel);
    spc = spc(sel);
    
    range = get(handles.edtNRange, 'String');
    [f spc] = extrange(f,spc,range);

    [M0, M1, M2] = nrNMRanalysis(f, spc, blcorr, cutoff);
    
    spc = real(spc);
    set(handles.txtM0,'String',['M0  = ' num2str(M0)]);
    set(handles.txtM1,'String',['M1  = ' num2str(M1)]);
    set(handles.txtM2,'String',['M2 = ' num2str(2*M2)]); 
                % peak to peak for gaussian is 2*sM2!!!
                % 2*M2 contains 68% of gaussian signal (like 2/3)
%     if plotErr
%         errorbar(f, spc, err)
%     else
        plot(f, spc)
%     end
    hold on
    xl = ylim;
    plot([M1 M1], [xl(1), xl(2)],'k-..','Linewidth',1.5)
    plot([M1-M2 M1-M2], [xl(1), xl(2)],'r--','Linewidth',1)
    plot([M1+M2 M1+M2], [xl(1), xl(2)],'r--','Linewidth',1)
    plot([min(f) max(f)], [cutoff, cutoff],'g-','Linewidth',1)
    grid on
    hold off
    
    % Save results
    nmr.nra{idx}.cutoff = cutoff;
    nmr.nra{idx}.range = range;
    nmr.nra{idx}.bline_corr = blcorr;  % kHz
    nmr.nra{idx}.results = [nmr.temp(idx) M0 0 M1 0 2*M2 0];
    assignin('base','nmr',nmr);

    

    
% --- Executes on button press in butNRArun.
function butNRArun_Callback(hObject, eventdata, handles)
    
    idx = get(handles.lbFiles, 'Value');
    max_idx = numel(get(handles.lbFiles, 'String'));
    nmr = evalin('base','nmr');

    while idx <= max_idx
        set(handles.lbFiles, 'Value',idx);
        nmrplot(nmr,idx,handles);
        getframe;
        butNAnalyse_Callback(hObject, eventdata, handles);
        idx = idx+1;
    end

    
    
% --- Executes on button press in butNRAplot.
function butNRAplot_Callback(hObject, eventdata, handles)
    nmr = evalin('base','nmr');
    
    if ~isfield(nmr,'nra')
        errordlg('No NRA results to plot!','Error')
        return;
    end
    
    if evalin('base','isfield(nmr,''material'')')
        title_str = [nmr.material '  ' nmr.date];
    else
        title_str = [nmr.dates{1}];
    end
    
    results = []; i = 1;
    while (i <= nmr.N)
        try
            results = [results; nmr.nra{i}.results];
        catch e
        end
        i=i+1;
    end
    
    plot_nmresults(results,title_str);



% % --- Executes on button press in btnNRArun.
% function btnNRArun_Callback(hObject, eventdata, handles)
% 

% % --- Executes during object creation, after setting all properties.
% function popWmethod_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to popWmethod (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: popupmenu controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 



% % --- Executes during object creation, after setting all properties.
% function popXcmethod_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to popXcmethod (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: popupmenu controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% % --- Executes on button press in butDelete.
% function butDelete_Callback(hObject, eventdata, handles)




% --- Executes during object creation, after setting all properties.
function edtCutOff_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end





% --- Executes on button press in chkPlotSim.
function chkPlotSim_Callback(hObject, eventdata, handles)


% --- Executes on button press in butReload.
function butReload_Callback(hObject, eventdata, handles)

    nmr = evalin('base','nmr');
    if isfield(nmr,'data')
        for i=1:numel(nmr.data)
            [pathstr, name, ext] = fileparts(nmr.data{i}.fname);
            shownames{i} = name;
        end
        set(handles.lbFiles, 'String', shownames);
    else
        set(handles.lbFiles, 'String', []);
    end



% --- Executes on button press in buttNRA.
function buttNRA_Callback(hObject, eventdata, handles)
    set(handles.panNRA,'Visible','on');
    set(handles.panFit,'Visible','off');
    set(handles.panSim,'Visible','off');

% --- Executes on button press in buttFitting.
function buttFitting_Callback(hObject, eventdata, handles)
    set(handles.panFit,'Visible','on');
    set(handles.panNRA,'Visible','off');
    set(handles.panSim,'Visible','off');

% --- Executes on button press in buttSimul.
function buttSimul_Callback(hObject, eventdata, handles)
    set(handles.panSim,'Visible','on');
    set(handles.panNRA,'Visible','off');
    set(handles.panFit,'Visible','off');




% --------------------------------------------------------------------
function mnuLoad_Callback(hObject, eventdata, handles)

    [FileName,PathName] = uigetfile('*.*','MultiSelect','on');
    if ~isequal(FileName, 0)
        cd(PathName)
        nmr = evalin('base','nmr');
        if iscell(FileName)
            for i=1:numel(FileName)
                nmr.path = [PathName FileName{i}];
                NMR_load;
            end
        else
            nmr.path = [PathName FileName];
            NMR_load;
            prompt = {'Temperature:','Frequency:'};
            dlg_title = 'Load nmr data';
            num_lines = 1;
            t = nmr.data{nmr.N}.exp.Temperature;
            f = nmr.data{nmr.N}.exp.Freq;
            def = {num2str(t),num2str(f)};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            nmr.data{nmr.N}.exp.Temperature = str2double(answer{1});
            nmr.data{nmr.N}.exp.Freq = str2double(answer{2});
            nmr.temp(nmr.N) = str2double(answer{1});
            nmr.freq(nmr.N) = str2double(answer{2});
        end
%         nmrplot(nmr,nmr.N,handles);
%         nmr_update(handles);
        assignin('base','nmr',nmr);
        butReload_Callback(hObject, eventdata, handles);
        set(handles.lbFiles, 'Value',nmr.N);
        butFFT_Callback(hObject, eventdata, handles);
    end



    
    
% --------------------------------------------------------------------
function mnuLoadSpc_Callback(hObject, eventdata, handles)
    
    [FileName,PathName] = uigetfile('*.*','MultiSelect','on');
    if ~isequal(FileName, 0)
        fname = [PathName FileName];
        cd(PathName)
        nmr = evalin('base','nmr');
        if iscell(FileName)
            prompt = {'Temperature:','Frequency:','Xlim:','isPPM (1/0)'};
            dlg_title = 'Load nmr spectrum';
            num_lines = 1;

            def = {num2str(300),num2str(95.575),'[-200 600]','1'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            for i=1:numel(FileName)
                fname = [PathName FileName{i}];
                if isfield(nmr,'N'), ind = nmr.N + 1; 
                else ind = 1; end;
                exp.Temperature = str2double(answer{1});
                exp.Freq = str2double(answer{2});
                exp.DW = 1;
                nmr.data{ind}.T = [0 0 0];
                nmr.data{ind}.signal = [0 0 0];
                nmr.data{ind}.exp = exp;
                nmr.data{ind}.fname = fname;
                nmr.temp(ind) = exp.Temperature;
                nmr.freq(ind) = exp.Freq;
                nmr.dates{ind} = datestr(now);

                data = load(fname,'ascii');
                [f id] = sort(data(:,1));
                spc = data(id,2);
                
                
                if size(data,2)>2
                   spc = complex(spc,data(:,3)); 
                end

                
                nmr.fft{ind}.spc = spc;
                nmr.fft{ind}.TD = numel(f);
                nmr.fft{ind}.shl = 0;
                nmr.fft{ind}.phs = 0;
                nmr.fft{ind}.de = 0;
                nmr.fft{ind}.lb = 0;
                nmr.fft{ind}.xlim = answer{3};
                nmr.fft{ind}.isPPM = str2num(answer{4});
%                 if nmr.fft{ind}.isPPM == 0
%                     nmr.fft{ind}.f = f*1000;
%                 else
                nmr.fft{ind}.f = f;
%                 end
                nmr.N = ind; 
                nmr.fft{ind}.autoshl = 0;
                nmr.fft{ind}.isNorm = 0;
                nmr.fft{ind}.REF = exp.Freq;
            end
        else
            if isfield(nmr,'N'), ind = nmr.N + 1; 
            else ind = 1; end;
            prompt = {'Temperature:','Frequency:','Xlim:','isPPM (1/0)'};
            dlg_title = 'Load nmr spectrum';
            num_lines = 1;

            def = {num2str(300),num2str(95.575),'[-50 50]','1'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            exp.Temperature = str2double(answer{1});
            exp.Freq = str2double(answer{2});
            exp.DW = 1;
            nmr.data{ind}.T = [0 0 0];
            nmr.data{ind}.signal = [0 0 0];
            nmr.data{ind}.exp = exp;
            nmr.data{ind}.fname = fname;
            nmr.temp(ind) = str2double(answer{1});
            nmr.freq(ind) = str2double(answer{2});
            nmr.dates{ind} = datestr(now);
            
            data = load(fname,'ascii');
            [f id] = sort(data(:,1));
             spc = data(id,2);
            if size(data,2)>2
               spc = complex(spc,data(:,3)); 
            end
            
            nmr.fft{ind}.f = f;
            nmr.fft{ind}.spc = spc;
            nmr.fft{ind}.TD = numel(f);
            nmr.fft{ind}.shl = 0;
            nmr.fft{ind}.phs = 0;
            nmr.fft{ind}.de = 0;
            nmr.fft{ind}.lb = 0;
            nmr.fft{ind}.xlim = answer{3};
            nmr.fft{ind}.isPPM = str2num(answer{4});
%             if nmr.fft{ind}.isPPM == 0
%                 nmr.fft{ind}.f = f*1000;
%             else
            nmr.fft{ind}.f = f;
%             end
            nmr.N = ind;
            nmr.fft{ind}.autoshl = 0;
            nmr.fft{ind}.isNorm = 0;
            nmr.fft{ind}.REF = exp.Freq;
        end
        assignin('base','nmr',nmr);
        butReload_Callback(hObject, eventdata, handles);
        set(handles.lbFiles, 'Value',nmr.N);
        set(handles.chkFFTSpc,'Value',1);
        nmrplot(nmr,nmr.N,handles);
        nmr_update(handles);
    end
    
    

% --------------------------------------------------------------------
function mnuSampleD_Callback(hObject, eventdata, handles)

    nmr = evalin('base','nmr');
    prompt = {'Material:','Mass [mg]:','Date:'};
    dlg_title = 'Sample Description';
    num_lines = 1;
    if isfield(nmr,'material')
        def = {nmr.material,num2str(nmr.mass),nmr.date};
    else
        def = {'C60','0.0',datestr(now, 'yyyy-mm-dd')};
    end
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if ~isempty(answer) % If CANCEL answer is empty
        nmr.material = answer{1};
        nmr.mass = answer{2};
        nmr.date = answer{3};
    end

    assignin('base','nmr',nmr);


% --- Executes on button press in butEdit.
function butEdit_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function mnuSort_Callback(hObject, eventdata, handles)
% hObject    handle to mnuSort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nmr = evalin('base','nmr');
    
prompt = {'By (temp/freq/dates):','Direction (ascend/descend):'};
dlg_title = 'Sort data';
num_lines = 1;
sortby = 'temp';
sortdir = 'descend';

if isfield(nmr,'sort')
    sortby = nmr.sort;
end
if isfield(nmr,'sortDir')
    sortdir = nmr.sortDir;
end

def = {sortby,sortdir};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if numel(answer)==0
    return;
end

nmr.sort = answer{1};
nmr.sortDir = answer{2};
nmr_sort;

assignin('base','nmr',nmr);

butReload_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function toolSave_ClickedCallback(hObject, eventdata, handles)

    [FileName,PathName] = uiputfile('*.nfi');
    if ~isequal(FileName, 0)
        cd(PathName);
        nmr = evalin('base','nmr');
        nmr.runscript = get(handles.edtRunScript,'String');
        save([PathName FileName],'nmr');
    end


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)

    [FileName,PathName] = uigetfile('*.nfi');
    if ~isequal(FileName, 0)
        cd(PathName)
        load([PathName FileName],'-mat');

        if ~exist('nmr','var')
            errordlg('Wrong file format!');
            return
        end
        
        set(handles.lbFiles, 'Value',nmr.N);
        if isfield(nmr,'runscript')
            set(handles.edtRunScript,'String',nmr.runscript);
        end
        
        if isfield(nmr,'nra')
            if ~iscell(nmr.nra)   % recover data if old format
                if isfield(nmr.nra,'cutoff')
                    cutoff = nmr.nra.cutoff;
                else
                    cutoff = -Inf;
                end
                if isfield(nmr.nra,'range') 
                   range = nmr.nra.range; 
                else 
                   range = '[-Inf Inf]'; 
                end
                if isfield(nmr.nra,'bline_corr') 
                   bline_corr = nmr.nra.bline_corr; 
                else
                   bline_corr = [0 0];
                end
                
                if isfield(nmr.nra,'results') 
                    for i=1:size(nmr.nra.results,1)
                        if nmr.nra.results(i,1) ~= 0
                            nraa{i}.cutoff = cutoff;
                            nraa{i}.range = range;
                            nraa{i}.bline_corr = bline_corr;
                            nraa{i}.results = nmr.nra.results(i,:);
                        else
                            nraa{i} = [];
                        end
                    end
                    nmr.nra = nraa;
                end

            end
        end
        assignin('base','nmr',nmr);
        butReload_Callback(hObject, eventdata, handles);
        nmr_update(handles);
    end



% --------------------------------------------------------------------
function mnuExpSim_Callback(hObject, eventdata, handles)

    [FileName,PathName] = uiputfile('*.dat');
    if ~isequal(FileName, 0)
        cd(PathName);
        nmr = evalin('base','nmr');
        idx = get(handles.lbFiles, 'Value');
        export_fit([PathName FileName],nmr,idx);
    end
    msgbox(['Fit has been saved to ' [PathName FileName]],'Saved');


% --------------------------------------------------------------------
function mnuExpAFit_Callback(hObject, eventdata, handles)

    nmr = evalin('base','nmr');
    names = get(handles.lbFiles, 'String');
    FileName = ['sim' names{1} '.dat'];
    [FileName,PathName] = uiputfile('*.dat','Save to directory',FileName);
    if ~isequal(FileName, 0)
        ret = questdlg('All files will be saved in this directory with an appendix ''sim...'' .','Confirm');
        if ~strcmp(ret,'Yes')
            return
        end
        cd(PathName);

        for i=1:numel(names)
            FileName = ['sim' names{i} '.dat'];
            export_fit([PathName FileName],nmr,i);
        end
    end
    msgbox(['Fits have been saved to ' PathName],'Saved');



% % User function: export selected fit
% function export_fit(Filename, nmr, idx)
% 
%     f = nmr.sim{idx}.f';
%     spc = nmr.sim{idx}.spc';
%     F = nmr.sim.fits{idx}.f(H)';
% 
%     fid = fopen(Filename, 'wt');
%     fprintf(fid, '# nmrFit 2.0  exported data+fit file;  Anton Potocnik @ IJS F5\n');
%     fprintf(fid, '# ------------------------------------------------------------\n');
%     fprintf(fid, '# Original data file:\n');
%     fprintf(fid, '# %s\n',nmr.data{idx}.fname);
%     fprintf(fid, '# Fit function:\n');
%     fprintf(fid, '# %s\n',nmr.sim.fitfun);
%     fprintf(fid, '# Coefficients:\n');
%     tmp = fieldnames(nmr.sim.coef);
%     for i=1:size(tmp,1)
%         fprintf(fid, '#  %s\t',tmp{i});
%         fprintf(fid, '%f\t+-%f\n',[nmr.sim.results(idx,i*2);nmr.sim.results(idx,i*2+1)]);
%     end
%     fprintf(fid, '# Data:\n');
%     fprintf(fid, '# X\t \tY\t \tFit\n');
% 
%     fprintf(fid, '%e\t%e\t%e\n', [H; Y; F]);
%     fclose(fid);



% User function: export selected fit
function export_NRAres(Filename, nmr)

%     T = nmr.nra.results(:,1);
%     M0 = nmr.nra.results(:,2);
%     M1 = nmr.nra.results(:,4);
%     M2 = nmr.nra.results(:,6);
% 
%     fid = fopen(Filename, 'wt');
% %     fprintf(fid, '# nmrFit 2.0  exported NRA results;  Anton Potocnik @ IJS F5\n');
% %     fprintf(fid, '# ----------------------------------------------------------\n');
% %     fprintf(fid, '# Material/Mass/Date:\n');
% %     fprintf(fid, '# %s\n',[nmr.material '  ' num2str(nmr.mass) 'mg' '  ' nmr.date]);
% %     fprintf(fid, '# W method:\t%s\n',nmr.nra.w_method);
% %     fprintf(fid, '# Xc method:\t%s\n',nmr.nra.xc_method);
% %     fprintf(fid, '# Integration Cut Off:\t%s\n',num2str(nmr.nra.cutoff));
% % 
% %     fprintf(fid, '# Data:\n');
%     fprintf(fid, 'T\t \tA\t \tdH\t \tg\t \txc\n');
% 
%     fprintf(fid, '%f\t%f\t%f\t%f\t%f\n', [T'; M0'; M1'; M2']);
%     fclose(fid);

    results = []; i = 1;
    while (i <= nmr.N)
        try
            results = [results; nmr.nra{i}.results];
        catch e
        end
        i=i+1;
    end
    
    num2clip(results)
    save(Filename,'results','-ascii')



% User function: export selected fit
function export_Simres(Filename, nmr, handles)

    idx = get(handles.lbFiles, 'Value');
    max_idx = numel(get(handles.lbFiles, 'String'));

    Results = [];
    for j = 1:max_idx
        if isempty(nmr.sim{j})
           continue 
        end
        result = [nmr.temp(j)];
        for i = 1:size(nmr.sim{j}.coefs,1)
            result = [result nmr.sim{j}.coefs{i,2}];
        end
        Results = [Results; result];
    end

    num2clip(Results)
    save(Filename,'Results','-ascii')


% --------------------------------------------------------------------
function mnuNRAres_Callback(hObject, eventdata, handles)

    nmr = evalin('base','nmr');
    if ~isfield(nmr,'nra')
        errordlg('NRA results not available!');
        return
    end

    [FileName,PathName] = uiputfile('*.dat');
    if ~isequal(FileName, 0)
        cd(PathName);
        export_NRAres([PathName FileName],nmr);
    end
    msgbox(['NRA results have been copied and saved to ' [PathName FileName]],'Saved');


% --------------------------------------------------------------------
function mnuSimRes_Callback(hObject, eventdata, handles)

    nmr = evalin('base','nmr');
    if ~isfield(nmr,'sim')
        errordlg('Fit results not available!');
        return
    end

    [FileName,PathName] = uiputfile('*.dat');
    if ~isequal(FileName, 0)
        cd(PathName);
        export_Simres([PathName FileName],nmr,handles);
    end
    msgbox(['Results copied and saved to ' [PathName FileName] '.']);




% --------------------------------------------------------------------
function mnuNew_Callback(hObject, eventdata, handles)
    answ = questdlg('All data in base workspace will be deleted! Do you really want to proceed?','New Experiment');

    if strcmp(answ,'Yes')
        evalin('base','clear nmr');
        assignin('base','nmr',[]);
        mnuSampleD_Callback(hObject, eventdata, handles)
        butReload_Callback(hObject, eventdata, handles);
        axes(handles.axes1);
        cla;
    end


% --------------------------------------------------------------------
function conEdit_Callback(hObject, eventdata, handles)

    idx = get(handles.lbFiles, 'Value');
    names = get(handles.lbFiles, 'String');
    if numel(names)<1
        return;
    end

    nmr = evalin('base','nmr');

    prompt = {'Temperature:','Frequency:'};
    dlg_title = names{idx};
    num_lines = 1;
    t = nmr.temp(idx);
    f = nmr.freq(idx);

    def = {num2str(t),num2str(f)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if numel(answer)==0
        return;
    end

    nmr.data{idx}.exp.Temperature = str2double(answer{1});
    nmr.data{idx}.exp.Freq = str2double(answer{2});
    nmr.temp(idx) = str2double(answer{1});
    nmr.freq(idx) = str2double(answer{2});

    assignin('base','nmr',nmr);

    butReload_Callback(hObject, eventdata, handles);
    set(handles.lbFiles, 'Value',nmr.N);
    nmrplot(nmr,nmr.N,handles);
    nmr_update(handles);

    
% --------------------------------------------------------------------
function conDelete_Callback(hObject, eventdata, handles)

    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');

    if nmr.N < 1
        return; % Not selected
    end

    if strcmp(questdlg('Do you really want to delete this measurement?','Delete'), 'Yes') ~= 1;
        return
    end

    nmr = nmr_delete(nmr,idx);
    assignin('base','nmr',nmr);

    if nmr.N < 1
        set(handles.lbFiles, 'Value',0);
    else
        if idx > nmr.N 
            set(handles.lbFiles, 'Value',nmr.N);
        else
            set(handles.lbFiles, 'Value',idx);
        end
    end

    butReload_Callback(hObject, eventdata, handles);



% --- Executes during object creation, after setting all properties.
function edtSHL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtSHL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function edtPhase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butFFT.
function butFFT_Callback(hObject, eventdata, handles)

    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');

    signal = nmr.data{idx}.signal;
    DW = nmr.data{idx}.exp.DW;
    TD = nmr.data{idx}.exp.TD;
    
    SHL = str2double(get(handles.edtSHL,'String'));
    autoSHL = get(handles.chkSHLauto,'Value');
    DE = get(handles.chkDE,'Value'); % Fill with signal
    LBexp = str2double(get(handles.edtExpBroad,'String'))*1000; % kHz
    LBgau = str2double(get(handles.edtGauBroad,'String'))*1000; % kHz
    userTD = str2double(get(handles.edtTD,'String'));
    PHASE = str2double(get(handles.edtPhase,'String'))/180*pi;
    autoPhase = get(handles.chkAutoPhase,'Value');
    
    autoSHLsquare = get(handles.chkSHLsquare,'Value');
    avrPoints = str2double(get(handles.edtAvrPoints,'String'));
    
    isNorm = get(handles.chkNorm,'Value');
    is380 = get(handles.chk380,'Value');
    REF = str2double(get(handles.edtREF,'String')); % MHz
    isPPM = get(handles.chkPPM,'Value');
    
    Range = str2num(get(handles.edtXLim,'String'));
    autoRange = get(handles.chkAutoRange,'Value');
    
    if autoSHL, SHL = sqrt(-1); end;
    if autoSHLsquare, SHL = -SHL; end;
    if autoPhase, PHASE = sqrt(-1); end;
    
    if userTD > TD
        TD = userTD;
        signal = [signal; zeros(TD-numel(signal),1)];
    end
    
    LB = complex(LBexp, LBgau); % Re:Lorentzian, Im: Gaussian broadening
    
    if is380, signal = conj(signal); end % Reverse X axis
    
    
    % ================ perform full FFT analysis ==========================
    
    [spc f PHASE SHL signal_new] = nmrFFT(signal,SHL,PHASE,LB,DE,TD,DW,avrPoints);
    
    % =====================================================================
    
    if isPPM
        f = (f + (nmr.freq(idx)-REF)*1e6)/REF; % to ppm, REF is in MHz, f is in Hz
    end
    
    if autoRange
        Range = round(autoRangeFun(spc,f)/1e3);
        if Range(2)-Range(1) < 1
            Range = [-80 80];
        end
        if isPPM, set(handles.edtXLim,'String', ['[' num2str(Range*1e3) ']']); end;
        if ~isPPM, set(handles.edtXLim,'String', ['[' num2str(Range) ']']); end;
    end

    [ind] = find(f<Range(1)*1e3|f>Range(2)*1e3);
    f(ind)=[];
    spc(ind)=[];
    
    if isNorm
        spc = spc/max(real(spc));
    end
    
%     nmr.data{idx}.signal = nmr.data{idx}.signal*exp(sqrt(-1)*PHASE); % correct signal
    nmr.data{idx}.signal_new = signal_new;
    nmr.fft{idx}.f=f;
    nmr.fft{idx}.spc=spc;
    nmr.fft{idx}.shl=SHL;
    nmr.fft{idx}.autoshl=autoSHL;
    nmr.fft{idx}.autophs=autoPhase;
    nmr.fft{idx}.turn=is380;
    nmr.fft{idx}.phs=PHASE; %phs0+PHASE
    nmr.fft{idx}.de=DE;
    nmr.fft{idx}.lb=LB;
    nmr.fft{idx}.TD=TD;
    nmr.fft{idx}.REF=REF;
    nmr.fft{idx}.isPPM=isPPM;
    nmr.fft{idx}.isNorm=isNorm;
    nmr.fft{idx}.autoRange=autoRange;
    nmr.fft{idx}.xlim=get(handles.edtXLim,'String');
    
    nmr.fft{idx}.autoSHLsquare = autoSHLsquare;
    nmr.fft{idx}.avrPoints = avrPoints;
    
    set(handles.edtPhase,'String',round(nmr.fft{idx}.phs/pi*180)); % Show phase
    set(handles.edtSHL,'String',SHL); % Show SHL
    set(handles.sldPhase,'Value',100); % Reset phase slider
    if isfield(nmr.fft{idx},'phs0')
        nmr.fft{idx} = rmfield(nmr.fft{idx},'spc0');
        nmr.fft{idx} = rmfield(nmr.fft{idx},'phs0');
    end
        
    assignin('base','nmr',nmr);
    
    %     set(handles.chkFFTSpc,'Value',1); % Show spectrum
    set(handles.tbutSPC,'value',1);
    set(handles.tbutFID,'value',0);
    set(handles.tbutmFID,'value',0);
    
    nmrplot(nmr,idx,handles);
    nmr_update(handles);


function Range = autoRangeFun0(spc, f)
    spc = abs(real(spc));
    M0 = trapz(spc);
    M1 = trapz(spc.*f)/M0;
    M2 = trapz(spc.*(f-M1).*(f-M1))/M0;
    w = sqrt(M2);
    [m max_i] = min(abs(f-(w+M1)));
    [m min_i] = min(abs(f-(-w+M1)));
    Range = [f(min_i) f(max_i)];
    
function Range = autoRangeFun(spc, f)
    spc = abs(real(spc));
    f0 = f;
    [ma M1] = max(spc);
    in = spc>(ma/5);
    spc = spc(in);
    f = f(in);
    
    M0 = trapz(spc);
    M1 = trapz(spc.*f)/M0;
    M2 = trapz(spc.*(f-M1).*(f-M1))/M0;
    w = sqrt(M2);
    [m max_i] = min(abs(f0-(10*w+M1)));
    [m min_i] = min(abs(f0-(-10*w+M1)));
    Range = [f0(min_i) f0(max_i)];
    
        
    
    
% --- Executes on button press in butPhaseCorr.
function butPhaseCorr_Callback(hObject, eventdata, handles)




% --------------------------------------------------------------------
function toolNew_ClickedCallback(hObject, eventdata, handles)

    mnuNew_Callback(hObject, eventdata, handles)





% --- Executes during object creation, after setting all properties.
function edtXLim_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function text27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes on button press in tbutFFTSpectrum.
function tbutFFTSpectrum_Callback(hObject, eventdata, handles)
% hObject    handle to tbutFFTSpectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes during object creation, after setting all properties.
function tbutFFTSpectrum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tbutFFTSpectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function edtSHL_Callback(hObject, eventdata, handles)
% hObject    handle to edtSHL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtSHL as text
%        str2double(get(hObject,'String')) returns contents of edtSHL as a double



function edtXLim_Callback(hObject, eventdata, handles)
% hObject    handle to edtXLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtXLim as text
%        str2double(get(hObject,'String')) returns contents of edtXLim as a double



function edtPhase_Callback(hObject, eventdata, handles)
% hObject    handle to edtPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtPhase as text
%        str2double(get(hObject,'String')) returns contents of edtPhase as a double


% --- Executes during object creation, after setting all properties.
function butPhaseCorr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to butPhaseCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuEPR_Callback(hObject, eventdata, handles)
% hObject    handle to mnuEPR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuExport_Callback(hObject, eventdata, handles)
% hObject    handle to mnuExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ListBox_Callback(hObject, eventdata, handles)
% hObject    handle to ListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function edtExpBroad_Callback(hObject, eventdata, handles)
% hObject    handle to edtExpBroad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtExpBroad as text
%        str2double(get(hObject,'String')) returns contents of edtExpBroad as a double


% --- Executes during object creation, after setting all properties.
function edtExpBroad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtExpBroad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in butSignal.
function butSignal_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function butDelete_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function butReload_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function butEdit_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function butNext_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function chkUpdate_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function chkPlotSim_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function btnSaveNRAH1_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function btnSaveNRAdH_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function btnSaveNRAH2_CreateFcn(hObject, eventdata, handles)

% --- Otherwise, executes on mouse press in 5 pixel border or over lbFiles.
function lbFiles_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes on key press with focus on lbFiles and none of its controls.
function lbFiles_KeyPressFcn(hObject, eventdata, handles)



% --------------------------------------------------------------------
function cmnCpySpc_Callback(hObject, eventdata, handles)

    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');
    ma = [];

    for i=idx
        T = nmr.temp(i);
        if isfield(nmr,'fft')
            if i<= numel(nmr.fft)
                f = [0; nmr.fft{i}.f];
                Y = [T; real(nmr.fft{i}.spc)];
                ma = [ma f Y];
            else
                disp('FFT does not exist for this file!')
            end 
        else
            disp('FFT field not found!')
        end
    end
    num2clip(ma)



% --------------------------------------------------------------------
function cmnCpyFit_Callback(hObject, eventdata, handles)

    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');
    ma = [];

    for i=idx
        if isfield(nmr,'sim')
            if i<= numel(nmr.sim)
                f = nmr.sim{i}.f;
                Y = real(nmr.fft{i}.spc);
                spc = nmr.sim{i}.spc;
                ma = [ma f spc]; %[ma f Y spc]; 
            else
                disp('Sim does not exist for this file!')
            end 
        else
            disp('Sim field not found!')
        end
    end
    num2clip(ma)




% --------------------------------------------------------------------
function mnuNormAll_Callback(hObject, eventdata, handles)
    nmr = evalin('base','nmr');
    max_idx = numel(get(handles.lbFiles, 'String'));
    
    r = inputdlg({'Normalize between [0 1] (norm) or devide by max value (max)?','Range:'},'Normalize',1,{'max',['1:' num2str(max_idx)]});
    if strcmp(r{1},'norm')
        norm = 1;
    else
        norm = 0;
    end

    for idx=str2num(r{2})
        spc = nmr.fft{idx}.spc;
        if norm == 1
            nmr.fft{idx}.spc = normalize(spc);
        else
            ma = max(real(nmr.fft{idx}.spc));
            nmr.fft{idx}.spc = nmr.fft{idx}.spc/ma;
        end
        if ~isfield(nmr.fft{idx},'spcOld')
            nmr.fft{idx}.spcOld = spc;
        end
    end
    assignin('base','nmr',nmr);

    set(handles.lbFiles, 'Value',idx);

    nmrplot(nmr,idx,handles);
    msgbox('Data hase been normalized. To undo use Renormalize All.','Normalize')

    
% --------------------------------------------------------------------
function mnuReNormAll_Callback(hObject, eventdata, handles)

    nmr = evalin('base','nmr');
    max_idx = numel(get(handles.lbFiles, 'String'));
    
    r = inputdlg('Range:','Normalize',1,{['1:' num2str(max_idx)]});
    done = [];
    idx_range = str2num(r{1});
    for idx=idx_range
        if ~isfield(nmr.fft{idx},'spcOld')
            disp(['idx = ' num2str(idx) ' renormalization not available!'])
            continue
        end
        nmr.fft{idx}.spc = nmr.fft{idx}.spcOld;
        done = [done idx];
    end
    assignin('base','nmr',nmr);

    set(handles.lbFiles, 'Value',idx);
    nmrplot(nmr,idx,handles);
    
    if numel(done) == numel(idx_range)
        msgbox('All spectra have been renormalized.','Renormalize');
    elseif numel(done) == 0
        msgbox('None of spectra have been renormalized.','Renormalize');
    else
        msgbox('Some spectra have been renormalized.','Renormalize');
    end


    

function Y = normalize(Y)
ymax = max(real(Y));
ymin = min(real(Y));

Y = (Y-ymin)/(ymax-ymin);




% --- Executes on button press in chkUpdate.
function chkUpdate_Callback(hObject, eventdata, handles)



% --------------------------------------------------------------------
function toolOpen_ClickedCallback(hObject, eventdata, handles)

    OpenMenuItem_Callback(hObject, eventdata, handles);



function edtNRange_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edtNRange_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function edtTau_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edtTau_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function edtN_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edtN_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function edtwQ_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edtwQ_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end



function edtT2_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edtT2_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end




% --- Executes on button press in butSimRun.
function butSimRun_Callback(hObject, eventdata, handles)




function edtA_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edtA_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edtIter_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edtIter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in butCpySimRes.
function butCpySimRes_Callback(hObject, eventdata, handles)




function edtSimRange_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edtSimRange_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in chkAutoPhase.
function chkAutoPhase_Callback(hObject, eventdata, handles)



% --- Executes on button press in butRescale.
function butRescale_Callback(hObject, eventdata, handles)

    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');
    xlstr = get(handles.edtXLim,'String');
    nmr.fft{idx}.xlim=xlstr;
    
    % cut data to range
    rng = str2num(xlstr);
    if nmr.fft{idx}.isPPM==0
        rng = rng*1e3;
    end
    iii = nmr.fft{idx}.f < rng(1) | nmr.fft{idx}.f > rng(2);
    nmr.fft{idx}.f(iii) = [];
    nmr.fft{idx}.spc(iii) = [];
    
    nmrplot(nmr,idx,handles);
    assignin('base','nmr',nmr);



% --- Executes on button press in butRescaleAll.
function butRescaleAll_Callback(hObject, eventdata, handles)

    max_idx = numel(get(handles.lbFiles, 'String'));
    nmr = evalin('base','nmr');
    xlstr = get(handles.edtXLim,'String');
    rng = str2num(xlstr);
    for idx=1:max_idx
        nmr.fft{idx}.xlim=xlstr;
        % cut data to range
        rng1 = rng;
        if nmr.fft{idx}.isPPM==0
            rng1 = rng*1e3;
        end
        iii = nmr.fft{idx}.f < rng1(1) | nmr.fft{idx}.f > rng1(2);
        nmr.fft{idx}.f(iii) = [];
        nmr.fft{idx}.spc(iii) = [];       
    end
    idx = get(handles.lbFiles, 'Value');
    nmrplot(nmr,idx,handles);
    assignin('base','nmr',nmr);




% --- Executes on button press in butFFTRun.
function butFFTRun_Callback(hObject, eventdata, handles)

    idx = get(handles.lbFiles, 'Value');
    max_idx = numel(get(handles.lbFiles, 'String'));
    nmr = evalin('base','nmr');

    while idx <= max_idx
        set(handles.lbFiles, 'Value',idx);
        nmrplot(nmr,idx,handles);
        getframe;
        butFFT_Callback(hObject, eventdata, handles);
        idx = idx+1;
    end

    


% --- Executes on button press in butScan.
function butScan_Callback(hObject, eventdata, handles)



function edtTD_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edtTD_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in chkFFTSpc.
function chkFFTSpc_Callback(hObject, eventdata, handles)

    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');
    nmrplot(nmr,idx,handles);



% --- Executes on button press in chkPltSig.
function chkPltSig_Callback(hObject, eventdata, handles)

% --- Executes on button press in rbuSignal.
function rbuSignal_Callback(hObject, eventdata, handles)

    if (get(hObject,'Value') == get(hObject,'Max'))
        set(handles.rbuZeros,'Value',0);
    end

% --- Executes on button press in rbuZeros.
function rbuZeros_Callback(hObject, eventdata, handles)

    if (get(hObject,'Value') == get(hObject,'Max'))
     	set(handles.rbuSignal,'Value',0);
    end

% --------------------------------------------------------------------
function mnuEdit_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function mnuSubstract_Callback(hObject, eventdata, handles)

% --- Executes on selection change in popSimLib.
function popSimLib_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popSimLib_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in butSimSet.
function butSimSet_Callback(hObject, eventdata, handles)

    sim_num = get(handles.popSimLib, 'Value');
    sim_list = get(handles.popSimLib, 'String');
    %idx = get(handles.lbFiles, 'Value');
    %nmr = evalin('base','nmr');
    
    [simu params startVal] = nmrsim_lib(sim_list{sim_num});
    coefs=[];
    for j=1:numel(params)
        variable = strtok(params{j}); % only the first word is variable name
        fix = false;
        if strcmp(variable,'TD') || strcmp(variable,'Spin') || strcmp(variable,'DW') || strcmp(variable,'N') || strcmp(variable,'NOP')  || strcmp(variable,'fmax')
            fix = true;
        end
        vrstica = {params{j},startVal(j),fix,'',false};
        coefs = [coefs; vrstica];
        set(handles.uitSimPar, 'Data', coefs);
    end


% --- Executes on button press in butSimulate.
function butSimulate_Callback(hObject, eventdata, handles)
    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');
    sim_num = get(handles.popSimLib, 'Value');
    sim_list = get(handles.popSimLib, 'String');
    range = get(handles.edtSimRange, 'String'); 
    disp('Sim...')

    coefs = get(handles.uitSimPar, 'Data');
    [nmrsim params] = nmrsim_lib(sim_list{sim_num});

    for i = 1:size(coefs,1)
        name = strtok(coefs{i,1}); % only the first word is variable name
        value = coefs{i,2};
        eval([name '=' num2str(value) ';']);
    end
    tic
    [spc f] = eval(nmrsim);
    toc
    for i = 1:size(coefs,1)
        name = strtok(coefs{i,1}); % only the first word is variable name
        eval(['nmr.sim{idx}.' name '=' name ';']);
    end
    
    
    % Calculate Chi2
    spc0 = real(nmr.fft{idx}.spc);
    f0 = real(nmr.fft{idx}.f);
    chi2 = calculateChi2(spc,f,spc0,f0,range);
    
    
    [f spc] = extrange(f,spc,range);
    nmr.sim{idx}.spc = spc;
    nmr.sim{idx}.f = f;
    nmr.sim{idx}.coefs = coefs;
    nmr.sim{idx}.sim_name = sim_list{sim_num};
    nmr.sim{idx}.chi2 = chi2;
    nmr.sim{idx}.range = get(handles.edtSimRange, 'String');
    disp('done')
    
    set(handles.chkPlotSim,'value',1);
    assignin('base','nmr',nmr);
    nmrplot(nmr,idx,handles);
    nmr_update(handles);

% Calculate Chi2
function chi2 = calculateChi2(spc,f,spc0,f0,range)

    [f spc] = extrange(f,spc,range);
    
    if numel(spc) > numel(spc0)
        spc0 = interp1(f0,spc0,f,'linear');
        f0 = f;
    end
    if numel(spc) < numel(spc0)
        spc = interp1(f,spc,f0,'linear');
        f = f0;
    end
    rm_ind = isnan(spc);
    spc(rm_ind) = [];
    f(rm_ind) = [];
    spc0(rm_ind) = [];
    f0(rm_ind) = [];
    rm_ind = isnan(spc0);
    spc(rm_ind) = [];
    f(rm_ind) = [];
    spc0(rm_ind) = [];
    f0(rm_ind) = [];
    
    [f0 spc0] = extrange(f0,spc0,range);
    [f spc] = extrange(f,spc,range);
    [m,n]=size(spc);
    spc0=reshape(spc0,m,n);
    chi2 = sum((spc-spc0).*(spc-spc0));

    
    
% --- Executes on button press in butFitSim.
function butFitSim_Callback(hObject, eventdata, handles)
    
    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');
    sim_num = get(handles.popSimLib, 'Value');
    sim_list = get(handles.popSimLib, 'String');
    Nrepeat = str2double(get(handles.edtRep,'String'));
    
    
    for l=1:Nrepeat
    disp('Fit...')

    coefs = get(handles.uitSimPar, 'Data');
    [nmrsim params] = nmrsim_lib(sim_list{sim_num});
    sim_name = nmrsim(1:strfind(nmrsim,'(')-1);

    % Set parameters & fitting parameter string including values
    param_string = '';
    param_idx = 1;
    x0 = [];
    Link = {};
    link_idx = 1;
   
    for i = 1:size(coefs,1)
        name = strtok(coefs{i,1}); % only the first word is variable name
        value = coefs{i,2};
        eval([name '=' num2str(value) ';']);
        
        if coefs{i,3} ~= 1 % it is not fixed

            isLinked = 0;
            if ~isempty(coefs{i,4}) % it is linked
                for j=1:numel(Link)
                   if strcmp(Link{j}{1},coefs{i,4})  % if you found a link
                       isLinked = Link{j}{2};
                       break
                   end
                end
                if isLinked <= 0
                    Link{link_idx} = {coefs{i,4},param_idx};
                    link_idx = link_idx + 1;
                end
            end
            if coefs{i,5} == 1 % it is log
                if isLinked > 0
                    param_string = [param_string ', exp(x(' num2str(isLinked) '))'];
                    param_idx = param_idx - 1;
                else
                    param_string = [param_string ', exp(x(' num2str(param_idx) '))'];
                    x0(param_idx) = log(value);
                end
            else
                if isLinked > 0
                    param_string = [param_string ', x(' num2str(isLinked) ')'];
                    param_idx = param_idx - 1;
                else
                    param_string = [param_string ', x(' num2str(param_idx) ')'];
                    x0(param_idx) = value;
                end
            end
            
            param_idx = param_idx + 1;
        else
            param_string = [param_string ', ' name];
        end
        
    end
    param_string(1)=[]; % delete first tick
    
    NIter = str2double(get(handles.edtIter,'String'));
    range = get(handles.edtSimRange, 'String'); 

    spc0 = real(nmr.fft{idx}.spc);
    f0 = real(nmr.fft{idx}.f);
    
    [f0 isort] = sort(f0);
    spc0 = spc0(isort);
%     if nmr.fft{idx}.isPPM == 1
%         f0 = f0*1000;
%     end
    opts = optimset('TolX',1e-10,'TolFun',1e-10,'Display','notify','MaxIter',NIter,'PlotFcns',@optimplotfval);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    results = eval(['fminsearch(@(x) Mini(spc0,f0,range,sim_name,{' param_string  '}),x0,opts);']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Set fitted values
    
    param_idx = 1;
    for i = 1:size(coefs,1)
        isLinked = 0;
        if coefs{i,3} ~= 1 % it is not fixed
            if ~isempty(coefs{i,4}) % it is linked
                for j=1:numel(Link) % Search for linkage
                   if strcmp(Link{j}{1},coefs{i,4})  % if you found a link
                       isLinked = Link{j}{2};
                       break
                   end
                end
            end
            name = strtok(coefs{i,1}); % only the first word is variable name
            if coefs{i,5} == 1 % it is log
                if isLinked <= 0
                    eval([name '= exp(results(' num2str(param_idx) '));']);
                else
                    eval([name '= exp(results(' num2str(isLinked) '));']);
                    if isLinked ~= param_idx, param_idx = param_idx - 1; end; % not for the first one linked
                end
            else
                if isLinked <= 0
                    eval([name '= results(' num2str(param_idx) ');']);
                else
                    eval([name '= results(' num2str(isLinked) ');']);
                    if isLinked ~= param_idx, param_idx = param_idx - 1; end; % not for the first one linked
                end
            end
            eval(['coefs{i,2}=' name ';']);
            param_idx = param_idx + 1;
        end
    end
    
    % Evaluate simulation
    [spc f] = eval(nmrsim);

    % Store parameters and spectrum
    for i = 1:size(coefs,1)
        name = strtok(coefs{i,1}); % only the first word is variable name
        eval(['nmr.sim{idx}.' name '=' name ';']);
    end
%%%
    nmr.sim{idx}.coefs = coefs;
    nmr.sim{idx}.sim_name = sim_list{sim_num};
    
    % Calculate Chi2
    spc0 = real(nmr.fft{idx}.spc);
    f0 = real(nmr.fft{idx}.f);
    chi2 = calculateChi2(spc,f,spc0,f0,range);
    
    nmr.sim{idx}.chi2 = chi2;
    nmr.sim{idx}.range = get(handles.edtSimRange, 'String');
    nmr.sim{idx}.spc = spc;
    nmr.sim{idx}.f = f;

    % Save and plot
    set(handles.txtChi, 'String', ['Chi2 = ' num2str(chi2)]);
    axes(handles.axes1); % It is making it very slow!!!
    set(handles.uitSimPar, 'Data', coefs);
    assignin('base','nmr',nmr);
    nmrplot(nmr,idx,handles);
    nmr_update(handles);
    disp('Finneshed!');
    
    end


function minimize = Mini(spc0,f0,range,sim_name,params)
    params_str = [];
    for i=1:numel(params)
        params_str = [params_str ', ' num2str(params{i})]; 
    end
    params_str(1) = [];
    [spc,f] = eval([ sim_name '(' params_str ');']);
    [f spc] = extrange(f,spc,range);
    
    if numel(spc) > numel(spc0)
        spc0 = interp1(f0,spc0,f,'linear');
        f0 = f;
    end
    if numel(spc) < numel(spc0)
        spc = interp1(f,spc,f0,'linear');
        f = f0;
    end
    rm_ind = isnan(spc);
    spc(rm_ind) = [];
    f(rm_ind) = [];
    spc0(rm_ind) = [];
    f0(rm_ind) = [];
    rm_ind = isnan(spc0);
    spc(rm_ind) = [];
    f(rm_ind) = [];
    spc0(rm_ind) = [];
    f0(rm_ind) = [];
    
    [f0 spc0] = extrange(f0,spc0,range);
    [f spc] = extrange(f,spc,range);
    [m,n]=size(spc);
    spc0=reshape(spc0,m,n);
    %plot(fQ,spc,'b',fQ,spc0,'k')
    % xlim([-1e5,1e5])
    %grid on
    %getframe;
    minimize = sum((spc-spc0).*(spc-spc0));
%     disp(['' params_str])




% --------------------------------------------------------------------
function mnuTestGauss_Callback(hObject, eventdata, handles)
    
    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');

    prompt = {'A:','x0:','w:'};
    dlg_title = 'Test Gaussian. Current data will be destroyed!!!';
    num_lines = 1;
    def = {'1','0','3'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if ~isempty(answer) % If CANCEL answer is empty
        A = str2double(answer{1});
        x0 = str2double(answer{2});
        w = str2double(answer{3});
    else 
        return
    end

    f = linspace(-5.0*w,5.0*w,1024);
    spc = A/sqrt(2*pi)/w*exp(-(f-x0).^2/2/w/w);

    nmr.fft{idx}.f = f;
    nmr.fft{idx}.spc = spc;
    
    assignin('base','nmr',nmr);
    nmrplot(nmr,idx,handles);




% --------------------------------------------------------------------
function mnuShifzSpc_Callback(hObject, eventdata, handles)
    nmr = evalin('base','nmr');
    idx_ = get(handles.lbFiles, 'Value');

    answer = inputdlg('Shift spc to the right (Hz):','Shift',1,{'0'});
    if numel(answer)==0
        return;
    end
    shift = evalin('base',answer{1});
    if numel(shift) == 1
        shift(idx_) = shift*ones(size(idx_));
    elseif numel(shift) ~= numel(idx_)
        msgbox('numel(shift) ~= numel(idx_) !!!')
        return
    end
    
    for idx = idx_
        if ~isfield(nmr.data{idx},'f0')
            nmr.fft{idx}.f0 = nmr.fft{idx}.f;
            nmr.fft{idx}.spc0 = nmr.fft{idx}.spc;
        end
        nmr.fft{idx}.f = nmr.fft{idx}.f + shift(idx)*ones(size(nmr.fft{idx}.f));
    end

    assignin('base','nmr',nmr);
    nmrplot(nmr,idx,handles);
    msgbox(['Specter ' num2str(idx) ' has been shifted!'],'Shift Data');


% --------------------------------------------------------------------
function mnuMultiply_Callback(hObject, eventdata, handles)
    nmr = evalin('base','nmr');
    idx_ = get(handles.lbFiles, 'Value');

    A = inputdlg('Multiply specter by a factor:','Multiply',1,{'0'});
    if numel(A)==0
        return;
    end

    for idx = idx_
        if ~isfield(nmr.data{idx},'spc0')
            nmr.fft{idx}.spc0 = nmr.fft{idx}.spc;
            nmr.fft{idx}.f0 = nmr.fft{idx}.f;
        end
        nmr.fft{idx}.spc = nmr.fft{idx}.spc * str2num(A{1});
    end
    
    assignin('base','nmr',nmr);
    nmrplot(nmr,idx,handles);
    msgbox(['Specter ' num2str(idx) ' has been Multiplied!'],'Mulitpy Data');


% --------------------------------------------------------------------
function mnuSubSpc_Callback(hObject, eventdata, handles)
    nmr = evalin('base','nmr');
    idx = get(handles.lbFiles, 'Value');

    prompt = {'idx =','factor ='};
    dlg_title = 'Substract reference data';
    num_lines = 1;
    index = idx;
    def = {num2str(index), '1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if numel(answer)==0
        return;
    end
    index = str2double(answer{1});
    factor = str2double(answer{2});
    % Write some protection if spc doesnot exits!

    f = nmr.fft{idx}.f;
    spc = nmr.fft{idx}.spc;

    if ~isfield(nmr.fft{idx},'spc0')
        nmr.fft{idx}.f0 = f;
        nmr.fft{idx}.spc0 = spc;
    end

    refSpc = nmr.fft{index}.spc;
    nmr.fft{idx}.spc = spc - factor*refSpc;
    assignin('base','nmr',nmr);
    nmrplot(nmr,idx,handles);
    
    
    
    % --------------------------------------------------------------------
function mnuSubSim_Callback(hObject, eventdata, handles)
    nmr = evalin('base','nmr');
    idx = get(handles.lbFiles, 'Value');

    prompt = {'idx =','factor ='};
    dlg_title = 'Substract Simulation';
    num_lines = 1;
    index = idx;
    def = {num2str(index), '1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if numel(answer)==0
        return;
    end
    index = str2double(answer{1});
    factor = str2double(answer{2});
    % Write some protection if spc doesnot exits!

    f = nmr.fft{idx}.f;
    spc = nmr.fft{idx}.spc;

    if ~isfield(nmr.fft{idx},'spc0')
        nmr.fft{idx}.f0 = f;
        nmr.fft{idx}.spc0 = spc;
    end

    Sim = nmr.sim{index}.spc;
    nmr.fft{idx}.spc = spc - factor*Sim;
    assignin('base','nmr',nmr);
    nmrplot(nmr,idx,handles);


    
% --------------------------------------------------------------------
function mnuSubBL_Callback(hObject, eventdata, handles)
    
    nmr = evalin('base','nmr');
    s = 0;    y0 = 0;    q = 0;
    prompt = {'q = ','s =','y0 ='};
    dlg_title = 'Substract y = q*x*x + s*x + y0';
    num_lines = 1;
    def = {num2str(q),num2str(s),num2str(y0)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if numel(answer)==0
        return;
    end
    q = str2double(answer{1});
    s = str2double(answer{2});
    y0 = str2double(answer{3});
    idx = get(handles.lbFiles, 'Value');
    
    if ~isfield(nmr.data{idx},'spc0')
        nmr.fft{idx}.spc0 = nmr.fft{idx}.spc;
        nmr.fft{idx}.f0 = nmr.fft{idx}.f;
    end

    f = nmr.fft{idx}.f;
    spc =  nmr.fft{idx}.spc;
    nmr.fft{idx}.spc = spc  - q*f.*f - s*f - y0;
    assignin('base','nmr',nmr);
    nmrplot(nmr,idx,handles);
    

    
% --------------------------------------------------------------------
function mnuResturn_Callback(hObject, eventdata, handles)
    nmr = evalin('base','nmr');
    idx_ = get(handles.lbFiles, 'Value');

    for idx = idx_
        if ~isfield(nmr.fft{idx},'spc0')
            msgbox('Nothing to return!','Error');
            return;
        end

        nmr.fft{idx}.f = nmr.fft{idx}.f0;
        nmr.fft{idx}.spc = nmr.fft{idx}.spc0;
    end

    assignin('base','nmr',nmr);
    nmrplot(nmr,idx,handles);
    msgbox(['Specter ' num2str(idx) ' has been restored!'],'Return Specter');





% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
k=0;
if strcmp(eventdata.Key,'downarrow')
    k = +1;
elseif strcmp(eventdata.Key,'uparrow')
    k = -1;
elseif strcmp(eventdata.Key,'leftarrow')
    k = -1;
elseif strcmp(eventdata.Key,'rightarrow')
    k = +1;
end
if strcmp(eventdata.Modifier,'control')
    nmr = evalin('base','nmr');
    idx = get(handles.lbFiles, 'Value');
    max_idx = numel(get(handles.lbFiles, 'String'));

    idx = idx + k;
    if idx(1) == 0 || idx(end) > max_idx
        return
    end
    
    set(handles.lbFiles, 'Value',idx);
    nmrplot(nmr,idx,handles);

end


% --------------------------------------------------------------------
function parseTemp_Callback(hObject, eventdata, handles)
nmr = evalin('base','nmr');
max_idx = numel(get(handles.lbFiles, 'String'));
names = get(handles.lbFiles, 'String');

for idx=1:max_idx
    tmp = names{idx};
    start = union(strfind(tmp,'-'),strfind(tmp,'_'));
    start = union(start,strfind(tmp,'.'));
    stop = union(strfind(tmp,'K'),strfind(tmp,'p0K'));
    stop = union(stop,strfind(tmp,'p5K'));
    min = 1000;
    istart = 0;
    istop = 0;
    t = 300;
    for i=1:numel(start)
        for j=1:numel(stop)
            k = stop(j)-start(i);
            if k<=1
                continue;
            end
            if k<min
                istart = start(i)+1;
                istop = stop(j)-1;
                min = k;
            end
        end
    end
    if min < 6 % nikoli ne bo šlo preko 10000K
        t = str2double(tmp(istart:istop));
        if isnan(t)
            t = 300;
        end
    end
    
    nmr.data{idx}.exp.Temperature = t;
    nmr.temp(idx) = t;
end
assignin('base','nmr',nmr);

butReload_Callback(hObject, eventdata, handles);
set(handles.lbFiles, 'Value',idx);
msgbox('Parsing Temperatures Done!','Parse Temperature')



% --- Executes on button press in chk380.
function chk380_Callback(hObject, eventdata, handles)
% hObject    handle to chk380 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk380


% --- Executes on button press in chkSHLauto.
function chkSHLauto_Callback(hObject, eventdata, handles)
% hObject    handle to chkSHLauto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkSHLauto


% --- Executes on button press in chkDE.
function chkDE_Callback(hObject, eventdata, handles)
% hObject    handle to chkDE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkDE



function edtREF_Callback(hObject, eventdata, handles)
% hObject    handle to edtREF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtREF as text
%        str2double(get(hObject,'String')) returns contents of edtREF as a double


% --- Executes during object creation, after setting all properties.
function edtREF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtREF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkPPM.
function chkPPM_Callback(hObject, eventdata, handles)
% hObject    handle to chkPPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkPPM


% --- Executes on button press in chkNorm.
function chkNorm_Callback(hObject, eventdata, handles)
% hObject    handle to chkNorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkNorm


% --------------------------------------------------------------------
function mnuCpySimRes_Callback(hObject, eventdata, handles)

    nmr = evalin('base','nmr');
    if ~isfield(nmr,'sim')
        errordlg('Sim results not available!');
        return
    end    

    idx = get(handles.lbFiles, 'Value');
    max_idx = numel(get(handles.lbFiles, 'String'));
    
    prompt = {'Range:'};
    dlg_title = 'Select results to copy';
    num_lines = 1;
    def = {['[1:' num2str(max_idx) ']']};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if ~isempty(answer) % If CANCEL answer is empty
        range = str2num(answer{1});
    else 
        return
    end

    Results = [];
    for j = range
        if j > numel(nmr.sim) || isempty(nmr.sim{j})
           continue 
        end
        result = [nmr.temp(j)];
        for i = 1:size(nmr.sim{j}.coefs,1)
            result = [result nmr.sim{j}.coefs{i,2}];
        end
        Results = [Results; result nmr.sim{j}.chi2];
    end

    num2clip(Results)
    
    msgbox(['Results copied.']);


% --- Executes during object creation, after setting all properties.
function edtBLCorr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtBLCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtBLCorr_Callback(hObject, eventdata, handles)
% hObject    handle to edtBLCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtBLCorr as text
%        str2double(get(hObject,'String')) returns contents of edtBLCorr as a double



function edtRep_Callback(hObject, eventdata, handles)
% hObject    handle to edtRep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtRep as text
%        str2double(get(hObject,'String')) returns contents of edtRep as a double


% --- Executes during object creation, after setting all properties.
function edtRep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtRep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togRunSim.
function togRunSim_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of togRunSim

    if ~get(hObject,'Value')
    	return 
    end

    idx = get(handles.lbFiles, 'Value');
    max_idx = numel(get(handles.lbFiles, 'String'));
    nmr = evalin('base','nmr');
    
    prompt = {'Index range:'};
    dlg_title = 'Run simulation fit';
    num_lines = 1;
    def = {[num2str(idx) ':1:' num2str(max_idx)]};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if ~isempty(answer) % If CANCEL answer is empty
        idx_range = str2num(answer{1});
%         max_idx = str2double(answer{2});
%         step = str2double(answer{3});
    else 
        return
    end
    
%     step = sign(max_idx-idx)*abs(step);

    for i_counter = idx_range %idx:step:max_idx
        set(handles.lbFiles, 'Value',i_counter);
        % --- place for changing the parameters ---
        script = get(handles.edtRunScript,'String');
        for j_counter=1:numel(script)
           eval(script{j_counter}) 
        end
        % -----------------------------------------
        
        butFitSim_Callback(hObject, eventdata, handles);
        pause(0.1)
        if ~get(hObject,'Value')
           return 
        end
    end
    set(hObject,'Value',0)






function edtRunScript_Callback(hObject, eventdata, handles)
% hObject    handle to edtRunScript (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtRunScript as text
%        str2double(get(hObject,'String')) returns contents of edtRunScript as a double


% --- Executes during object creation, after setting all properties.
function edtRunScript_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtRunScript (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)



% --- Executes on button press in chkAutoRange.
function chkAutoRange_Callback(hObject, eventdata, handles)
% hObject    handle to chkAutoRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkAutoRange


% --- Executes on slider movement.
function sldPhase_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');
    
    if ~isfield(nmr.fft{idx},'phs0'), nmr.fft{idx}.phs0=nmr.fft{idx}.phs; end;
    if ~isfield(nmr.fft{idx},'spc0'), nmr.fft{idx}.spc0=nmr.fft{idx}.spc; end;
    nmr.fft{idx}.phs = nmr.fft{idx}.phs0 + pi*0.005*(get(hObject,'Value')-100);

    nmr.fft{idx}.spc = nmr.fft{idx}.spc0.*exp(sqrt(-1)*pi*0.005*(get(hObject,'Value')-100));
    set(handles.edtPhase,'String',num2str(round(nmr.fft{idx}.phs/pi*180)))
    set(handles.chkAutoPhase,'Value',0);
    nmr.fft{idx}.autophs = 0;
    assignin('base','nmr',nmr);
    nmrplot(nmr,idx,handles)


% --- Executes during object creation, after setting all properties.
function sldPhase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sldPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over sldPhase.
function sldPhase_ButtonDownFcn(hObject, eventdata, handles)




% --- Executes on button press in chkErrBar.
function chkErrBar_Callback(hObject, eventdata, handles)
% hObject    handle to chkErrBar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkErrBar


% --- Executes on button press in butAnaErr.
function butAnaErr_Callback(hObject, eventdata, handles)


% --- Executes on button press in chkAnaErr.
function chkAnaErr_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function mnuOffset_Callback(hObject, eventdata, handles)
    nmr = evalin('base','nmr');
    idx_ = get(handles.lbFiles, 'Value');

    A = inputdlg({'Offset specter by:'},'Offset',1,{'T'});
    if numel(A)==0
        return;
    end
    
    for idx = idx_
        if ~isfield(nmr.data{idx},'spc0')
            nmr.fft{idx}.spc0 = nmr.fft{idx}.spc;
            nmr.fft{idx}.f0 = nmr.fft{idx}.f;
        end
        if strcmp(A{1},'T')
            nmr.fft{idx}.spc = nmr.fft{idx}.spc + nmr.temp(idx);
        elseif strcmp(A{1},'-T')
            nmr.fft{idx}.spc = nmr.fft{idx}.spc - nmr.temp(idx);  
        else
            nmr.fft{idx}.spc = nmr.fft{idx}.spc + str2double(A{1});
        end
    end
    
    assignin('base','nmr',nmr);
    nmrplot(nmr,idx,handles);
    msgbox(['Specter ' num2str(idx) ' has been Offseted!'],'Offset Data');



function edtCutOff_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function mnuReload_Callback(hObject, eventdata, handles)
    butReload_Callback(hObject, eventdata, handles)
    
    % Initialize popSumLib
    sims = nmrsim_lib('');  % Get simmulations name cell array
    set(handles.popSimLib, 'String',sims);
    
    % Check if there are any new files
    reloadFileListBox(handles);


% --------------------------------------------------------------------
function mnuCpyNRAres_Callback(hObject, eventdata, handles)
    nmr = evalin('base','nmr');
    
    if ~isfield(nmr,'nra')
        errordlg('No NRA results to copy!','Error');
        return;
    end
    
    results = []; i = 1;
    while (i <= nmr.N)
        try
            results = [results; nmr.nra{i}.results];
        catch e
        end
        i=i+1;
    end

    num2clip(results)
    msgbox('NRA results have been copied to clipboard.')
    


% --------------------------------------------------------------------
function mnuCpySpc3D_Callback(hObject, eventdata, handles)
    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');
    ma = [];

    for i=idx
        T = nmr.temp(i);
        if isfield(nmr,'fft')
            if i<= numel(nmr.fft)
                f = [nmr.fft{i}.f];
                Y = [real(nmr.fft{i}.spc)];
                ma = [ma; f T*ones(size(Y)) Y];
            else
                disp('FFT does not exist for this file!')
            end 
        else
            disp('FFT field not found!')
        end
    end
    num2clip(ma)




% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over txt13C400MHz.
function txt13C400MHz_ButtonDownFcn(hObject, eventdata, handles)


function txt87Rb400MHz_ButtonDownFcn(hObject, eventdata, handles)


function tct13C380MHz_ButtonDownFcn(hObject, eventdata, handles)


function txt133Cs400MHz_ButtonDownFcn(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function txt133Cs400MHz_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function tct13C380MHz_CreateFcn(hObject, eventdata, handles)
     

% --- Executes during object creation, after setting all properties.
function txt87Rb400MHz_CreateFcn(hObject, eventdata, handles)




function edtYrange_Callback(hObject, eventdata, handles)
% hObject    handle to edtYrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtYrange as text
%        str2double(get(hObject,'String')) returns contents of edtYrange as a double


% --- Executes during object creation, after setting all properties.
function edtYrange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtYrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkAutoYrange.
function chkAutoYrange_Callback(hObject, eventdata, handles)
% hObject    handle to chkAutoYrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkAutoYrange



function edtAvrPoints_Callback(hObject, eventdata, handles)
% hObject    handle to edtAvrPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtAvrPoints as text
%        str2double(get(hObject,'String')) returns contents of edtAvrPoints as a double


% --- Executes during object creation, after setting all properties.
function edtAvrPoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtAvrPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkSPCphase.
function chkSPCphase_Callback(hObject, eventdata, handles)
% hObject    handle to chkSPCphase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkSPCphase


% --- Executes on button press in chkSHLsquare.
function chkSHLsquare_Callback(hObject, eventdata, handles)
% hObject    handle to chkSHLsquare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkSHLsquare



function edtGauBroad_Callback(hObject, eventdata, handles)
% hObject    handle to edtGauBroad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtGauBroad as text
%        str2double(get(hObject,'String')) returns contents of edtGauBroad as a double


% --- Executes during object creation, after setting all properties.
function edtGauBroad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtGauBroad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10


% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butREFset.
function butREFset_Callback(hObject, eventdata, handles)
% hObject    handle to butREFset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butSPCphase.
function butSPCphase_Callback(hObject, eventdata, handles)



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over txtPointer1.
function txtPointer1_ButtonDownFcn(hObject, eventdata, handles)



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over txtPointer2.
function txtPointer2_ButtonDownFcn(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function txtPointer2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPointer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in butPointer2.
function pushbutton65_Callback(hObject, eventdata, handles)
% hObject    handle to butPointer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in tbutSPC.
function tbutSPC_Callback(hObject, eventdata, handles)

    if get(hObject,'Value') == 1
        set(handles.tbutFID,'Value',0);
        set(handles.tbutmFID,'Value',0);
    else
        set(hObject,'Value',1);
    end

    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');
    nmrplot(nmr,idx,handles);


% --- Executes on button press in tbutFID.
function tbutFID_Callback(hObject, eventdata, handles)
    if get(hObject,'Value') == 1
        set(handles.tbutSPC,'Value',0);
        set(handles.tbutmFID,'Value',0);
    else
        set(hObject,'Value',1);
    end

    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');
    nmrplot(nmr,idx,handles);


% --- Executes on button press in tbutmFID.
function tbutmFID_Callback(hObject, eventdata, handles)
    if get(hObject,'Value') == 1
        set(handles.tbutFID,'Value',0);
        set(handles.tbutSPC,'Value',0);
    else
        set(hObject,'Value',1);
    end

    idx = get(handles.lbFiles, 'Value');
    nmr = evalin('base','nmr');
    nmrplot(nmr,idx,handles);





% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
