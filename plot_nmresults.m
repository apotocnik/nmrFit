%--------------------------------------------------------------------------
% plot_nmresults  - plots multiple results: M0(T), M1(T), M2(T)
%
% Author: Anton Potocnik, F5, IJS
% Date:   27.01.2009 - 31.10.2009
% Arguments:
%       h = plot_nmresults(results,[title])
% Input:    {results} ... [T, M0, dM0, M1, dM1, M2, dM2]   cell array of results
%           title       ... title string or epr structure
%--------------------------------------------------------------------------

function plot_nmresults(results, title_str)
if nargin<2
    title_str = 'Title  Date';
end

if isstruct(title_str)
    if isfield(title_str,'material') && isfield(title_str,'date')
        title_str = [title_str.material '  ' title_str.date '  '];
    end
end

if ~iscell(results)
    results = {results};
end


%--------------------------------------------------------------------------
font_size_title = 18;
font_size_labels = 14;
font_size_numbers = 14;
small_font_size_labels = 10;
small_font_size_numbers = 10;
%--------------------------------------------------------------------------

% Create figure
figure1 = figure('Position',[560,50,560,670]);

%% Create axes1 for M0
axes1 = axes('Parent',figure1,'YMinorTick','on','XMinorTick','on',...
    'Position',[0.16 0.677 0.77 0.24],...
    'LineWidth',1.5,...
    'FontSize',font_size_numbers,...
    'FontName','Arial');
box('on');
hold('all');


% Create plot
for i=1:numel(results)
    % Get variables
    T = results{i}(:,1);
    M0 = results{i}(:,2);
    dM0 = results{i}(:,3);
    
    errorbar(T,M0,dM0,'Parent',axes1,'MarkerFaceColor',colors(i),...
    'MarkerEdgeColor',colors(i),'Marker',markers(i),'LineStyle','none');
end
grid on

% Create ylabel
ylabel('Intensity (arb. units)','FontSize',font_size_labels,'FontName','Arial');

% Create title
title(title_str,'FontSize',font_size_title);

% % Create xlabel
% xlabel('\itT \rm(K)','FontSize',font_size_labels,'FontName','Arial');


%% Create axes2 for M1
axes2 = axes('Parent',figure1,'YMinorTick','on','XMinorTick','on',...
    'Position',[0.16 0.377 0.77 0.24],...
    'LineWidth',1.5,...
    'FontSize',font_size_numbers,...
    'FontName','Arial');
box('on');
hold('all');

% Create plot
for i=1:numel(results)
    % Get variables
    T = results{i}(:,1);
    M1 = results{i}(:,4);
    dM1 = results{i}(:,5);
    
    errorbar(T,M1,dM1,'Parent',axes2,'MarkerFaceColor',colors(i),...
    'MarkerEdgeColor',colors(i),'Marker',markers(i),'LineStyle','none');
end
grid on

% Create ylabel
ylabel('M1 (kHz)','FontSize',font_size_labels,'FontName','Arial');

% % Create xlabel
% xlabel('\itT \rm(K)','FontSize',font_size_labels,'FontName','Arial');


%% Create axes3 for M2
axes3 = axes('Parent',figure1,'YMinorTick','on','XMinorTick','on',...
    'Position',[0.16 0.077 0.77 0.24],...
    'LineWidth',1.5,...
    'FontSize',font_size_numbers,...
    'FontName','Arial');
box('on');
hold('all');

% Create plot
for i=1:numel(results)
    % Get variables
    T = results{i}(:,1);
    M2 = results{i}(:,6);
    dM2 = results{i}(:,7);
    
    errorbar(T,M2,dM2,'Parent',axes3,'MarkerFaceColor',colors(i),...
    'MarkerEdgeColor',colors(i),'Marker',markers(i),'LineStyle','none');
end
grid on

% Create ylabel
ylabel('M2 (kHz)','FontSize',font_size_labels,'FontName','Arial');

% Create xlabel
xlabel('\itT \rm(K)','FontSize',font_size_labels,'FontName','Arial');




