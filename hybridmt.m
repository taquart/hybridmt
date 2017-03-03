function [Solution, History] = hybridmt(Input, varargin)
%HYBRIDMT Perform refinement of moment tensors using hybridMT technique.
%   Use HYBRIDMT to perform refinement of a group of moment tensor
%   solutions forming a tight cluster using Hybrid Moment Tensor technique.
%   The core implementation of moment tensor refinement techniques is based on
%   Andersen, L. M. (2001), "A relative moment tensor inversion technique
%   applied to seismicity induced by mining", PhD Thesis, Univ. of the
%   Witwatersrand, Johannesburg.
%
%   Solution = HYBRIDMT(Input) performs refinement of seismic moment tensors
%   provided in Input cell array using default parameters. The refined moment
%   tensors are returned in Solution cell array.
%
%   Solution = HYBRIDMT(Input, 'ParamName', ParamValue, ...) allows to specify
%   additional options in a form of ParamName-ParamValue.
%
%   See also FOCIMT
%
%   part of hybridMT package 
%   <a href="matlab:open('html/doc_hybridmt.html')">Reference page for hybridmt</a>

%   Copyright 2015-2016 Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%                       Patricia Martinez-Garzon <patricia@gfz-potsdam.de>
%
%   $Revision: 1.0.5 $  $Date: 2016.07.08 $

% Parse input parameters.
p = inputParser;
p.addRequired('Input', @(x) ischar(x) );
p.addParamValue('Iterations', 40, @(x) isscalar(x) && x > 0);
p.addParamValue('Weight', 0.1, @(x) isscalar(x) && x > 0);
p.addParamValue('Solution', 'full', @(x)any(strcmpi(x,{'full','deviatoric','dc'}))); %#ok<*NVREPL>
p.addParamValue('ProjectDir', 'results', @(x) ischar(x) );
p.addParamValue('Snapshots', 'on', @(x)any(strcmpi(x,{'off','on'}))); %#ok<*NVREPL>
p.addParamValue('Display', 'on', @(x)any(strcmpi(x,{'off','on'}))); %#ok<*NVREPL>
p.addParamValue('RatioLimit', 10, @(x) isscalar(x) && x >= 1);
p.addParamValue('TestData', 'off', @(x)any(strcmpi(x,{'off','on'})));
p.addParamValue('IgnoreStation', cell(0), @(x) iscell(x)); % passthrough to focimt.m
p.addParamValue('CorrectStation', cell(0), @(x) iscell(x));
p.addParamValue('VelocityModel', [], @(x) ismatrix(x)); % passthrough to focimt.m
p.addParamValue('FigureFormat', 'PNG', @(x)any(strcmpi(x,{'PNG','PS','SVG','PDF'})));

p.addParamValue('AuxiliaryFigures', 'on', @(x)any(strcmpi(x,{'off','on'}))); %#ok<*NVREPL>

p.parse(Input,varargin{:});

% Interpret input parameters.
n_iter = p.Results.Iterations; % number of iterations.
weight = p.Results.Weight; % weighting for displacement correction procedure
mt_solution_type = p.Results.Solution; % type of solution to consider.
output_dir = p.Results.ProjectDir; % Output directory name.
snapshots = strcmp(p.Results.Snapshots,'on');
ratio_limit = p.Results.RatioLimit;
display_progress = strcmp(p.Results.Display,'on');
test_data = strcmp(p.Results.TestData,'on');
StIgnored = p.Results.IgnoreStation;
CoStation = p.Results.CorrectStation;
picformat = lower(p.Results.FigureFormat);
min_npol = 10;
min_events = 20;
pick_ratio = 0.5;
polarity_match_ratio = 0.8;

%---- Setup properties for fociMT.
cmdline = {};
cmdline2 = {};
if ~isempty(StIgnored)
  cmdline = [cmdline 'IgnoreStation' {StIgnored}];
end
if snapshots
  cmdline = [cmdline 'BeachBallFormat' 'PNG'];
end
if ~isempty(p.Results.VelocityModel)
  cmdline = [cmdline 'VelocityModel' p.Results.VelocityModel];
  cmdline2 = [cmdline2 'VelocityModel' p.Results.VelocityModel];
end
% if ~isempty(CoStation)
%   cmdline = [cmdline 'CorrectStation' {CoStation}];
% end

disp('Preparing data for hybridMT inversion');

% Calculate first MT solution.
[Solution, Input] = focimt(Input,cmdline2{:}); % Calculate initial solution in order to get basic information.

% Prepare matrices for hybrid moment tensor by analysis of the MT solution.
n_events = numel(Solution); % Number of events.

if n_events < min_events
  fprintf('Low number of events taken for hybrid moment tensor inversion (%d)\n',n_events);
elseif n_events <= 1/ratio_limit
  error('Cannot perform hybrid moment tensor inversion with a single event.');
end

% Check for repetive IDs
IDListTemp = getsolution(Solution,mt_solution_type,'event_id');
if numel(IDListTemp) ~= numel(unique(IDListTemp))
  error('Input phase data contains event(s) with same event id(s).');
end

% Check how many stations we have.
Stations = cell(0);
P_MATCH = [];
BadEvents = cell(0);
BADEVENTSN = [];
k = 1;
for i=1:n_events
  event_id = Solution{i}.event_id;
  switch mt_solution_type
    case 'full'
      Sol = Solution{i}.full;
    case 'deviatoric'
      Sol = Solution{i}.deviatoric;
    case 'dc'
      Sol = Solution{i}.dc;
  end
  n = length(Sol.UTH);
  if n < min_npol
    BadEvents{k} = event_id;
    BADEVENTSN(k) = n; %#ok<AGROW>
    k = k + 1;
  end
  Stations = [Stations Sol.Station]; %#ok<AGROW>
  P_MATCH = [P_MATCH sign(Sol.UMEASURED.*Sol.UTH) == 1]; %#ok<AGROW> % count how many polarities match.
end

% "Stations" contain unique names of stations from all events.
Stations_unique = sort(unique(Stations));

% remove ignored stations.
if numel(StIgnored)
  for i=1:numel(StIgnored)
    I = ~strcmpi(Stations_unique,StIgnored{i});
    Stations_unique = Stations_unique(I);
  end
end

n_stations = numel(Stations_unique);
N_PICKS = zeros(1,n_stations);
N_MATCH = zeros(1,n_stations);
for i=1:n_stations
  N_PICKS(i) = sum(strcmpi(Stations,Stations_unique{i}));
  N_MATCH(i) = sum(P_MATCH(strcmpi(Stations,Stations_unique{i})));
end
Stations = Stations_unique;

%---- Display basic statistics.

% Event statistics.
fprintf(' Number of events: %d\n',n_events);
if numel(BadEvents)
  fprintf(' There are events with limited number of polarities (<%d):\n',min_npol);
  for i=1:numel(BadEvents);
    fprintf(' %s (N=%d)',BadEvents{i},BADEVENTSN(i));
    if mod(i,5) == 0
      fprintf('\n');
    end
  end
  fprintf('\n Consider removing them from inversion\r\n');
end
% Station statistics.
fprintf(' Number of stations: %d\n',n_stations);
ST = 1:n_stations;
idx_picks = [];
idx_match = [];
fprintf(' -----------------------------------------------------------------------------------\n');
while ~isempty(ST)
  if length(ST) >= 10
    ST_TEMP = ST(1:10);
    ST = ST(11:end);
  else
    ST_TEMP = ST(1:end);
    ST = [];
  end
  fprintf(' Station:     ');
  for i=ST_TEMP
    fprintf(' %6s',Stations{i});
  end
  fprintf('\n');
  fprintf(' %% picked:    ');
  for i=ST_TEMP
    v = N_PICKS(i)/n_events;
    if v < pick_ratio
      idx_picks = [idx_picks i]; %#ok<AGROW>
    end
    fprintf(' %6d',round(v*100));
  end
  fprintf('\n');
  fprintf(' %% p.matched: ');
  for i=ST_TEMP
    v = N_MATCH(i)/N_PICKS(i);
    if v < polarity_match_ratio
      idx_match = [idx_match i]; %#ok<AGROW>
    end
    fprintf(' %6d',round(v*100));
  end
  fprintf('\n -----------------------------------------------------------------------------------\n');
end
if ~isempty(idx_picks)
  fprintf('\n Warning: There are stations where amplitudes were picked rarely (<%1.0f%% of events):\n',pick_ratio*100);
  for i=idx_picks
    fprintf(' %s (%1.0f%%)',Stations{i},N_PICKS(i)/n_events*100);
  end
  fprintf('\n');
  fprintf(' Consider ignoring these stations in inversion (use ''IgnoreStation'' property)\n');
end
if ~isempty(idx_match)
  fprintf('\n Warning: There are stations where polarities frequently do not match (<%1.0f%% picks):\n',polarity_match_ratio*100);
  for i=idx_match
    fprintf(' %s (%1.0f%%)',Stations{i},N_MATCH(i)/N_PICKS(i)*100);
  end
  fprintf('\n');
  fprintf(' Consider ignoring these stations inversion (use ''IgnoreStation'' property)\n');
end
fprintf('\n');

if test_data
  return;
end

% Prepare history data.
History.RMS_ERROR = nan(n_events, n_iter);
History.U_RATIO_MEDIAN = nan(n_stations, n_iter); % changes in median U_RATIO with iteration.
History.U       = ones(n_stations,1); % current U multiplier.
History.U_CHANGES = ones(n_stations,n_iter); % cumulative history of U changes
History.U_RATIO = cell(1,n_iter); % changes of U_RATIO with iteration
History.POLARITY_VALID = ones(n_stations, n_iter); % changes in number of valid polarities with iteration
History.POLARITY_VALID_TOTAL = nan(1,n_iter);

if exist(output_dir,'dir')
  error('Project directory "%s" already exists. Please delete it manually.',output_dir);
end

% Correct ONCE input amplitudes.
if ~isempty(CoStation)
  for i=1:2:numel(CoStation)
    station = CoStation{i};
    factor = CoStation{i+1};
    for j=1:numel(Input)
      I = strcmpi(Input{j}.Station,station);
      if sum(I)
        Input{j}.OMEGA(I) = Input{j}.OMEGA(I) * factor;
      end
    end
  end
end

%---- Main loop of the program.
rmsg = '';
disp('Performing hybridMT inversion');
for iteration = 1:n_iter
  
  % Display progress information.
  progress = 100*iteration/n_iter;
  msg = sprintf(' Iteration: %3d/%3d |%s| (%1.1f%%)',iteration,n_iter,[repmat('o',1,floor(progress/5)) repmat('-',1,ceil((100-progress)/5))],progress);
  fprintf('%s',[rmsg, msg]); rmsg = repmat(sprintf('\b'), 1, length(msg));
  
  % Calculate moment tensor solutions for the whole dataset.
  if isempty(cmdline)
    [Solution, Input] = focimt(Input);
  else
    [Solution, Input] = focimt(Input,cmdline{:});
  end
  
  % Transfer figures and output data project folder.
  if ~exist(output_dir,'dir')
    mkdir(output_dir);
    if snapshots || strcmp(p.Results.AuxiliaryFigures,'on')
      for i=1:numel(Solution)
        mkdir(output_dir,Solution{i}.event_id);
      end
    end
    if snapshots
      for i=1:numel(Solution)
        mkdir(sprintf('%s/%s',output_dir,Solution{i}.event_id),'iter');
      end
    end
  end
  
  if snapshots % remove or move beach balls to appriopriate folders
    for i=1:numel(Solution)
      file = Solution{i}.event_id;
      switch mt_solution_type
        case 'full'
          movefileifexist([file '-full.png'],sprintf('%s/%s/iter/%03d.png',output_dir,file,iteration));
          deleteifexist([file '-deviatoric.png']);
          deleteifexist([file '-dc.png']);
        case 'deviatoric'
          movefileifexist([file '-deviatoric.png'],sprintf('%s/%s/iter/%03d.png',output_dir,file,iteration));
          deleteifexist([file '-full.png']);
          deleteifexist([file '-dc.png']);
        case 'dc'
          movefileifexist([file '-dc.png'],sprintf('%s/%s/iter/%03d.png',output_dir,file,iteration));
          deleteifexist([file '-full.png']);
          deleteifexist([file '-deviatoric.png']);
      end
    end
  end
  
  
  % This matrixes will hold theoretical and measured displacements.
  U_THEORETICAL   = nan(n_events, n_stations);
  U_MEASURED = nan(n_events, n_stations);
  
  for i=1:n_events
    
    % Get appropriate moment tensor solution.
    switch mt_solution_type
      case 'full'
        Sol = Solution{i}.full;
      case 'deviatoric'
        Sol = Solution{i}.deviatoric;
      case 'dc'
        Sol = Solution{i}.dc;
    end
    
    % Get station and displacement information
    Sta = Sol.Station;
    U_TH_TEMP = Sol.UTH;
    U_MEAS_TEMP = Sol.UMEASURED;
    
    % As for some stations there may be no displacement data, we have to
    % put data in appropriate places of the matrix.
    Sta_ID = nan(size(Sta));
    for j=1:numel(Sta)
      Sta_ID(j) = find(strcmp(Stations,Sta{j}));
    end
    if sum(isnan(Sta_ID))
      error('!!!!');
    end
    
    % Transfer displacements into appropriate place of matrix.
    U_THEORETICAL(i,Sta_ID) = U_TH_TEMP;
    U_MEASURED(i,Sta_ID) = U_MEAS_TEMP;
    
  end
  
  % Calculate ratio of amplitudes.
  U_RATIO = U_THEORETICAL ./ U_MEASURED;
  
  % Matrices for storing median displacement ratios and RMS errors.
  U_RATIO_MEDIAN = nan(1,size(U_RATIO,2));
  POLARITY_DISC  = nan(1,size(U_RATIO,2));
  RMS_ERROR = nan(size(U_RATIO,1),1);
  
  % Calculate displacement ratios for each station.
  for i=1:size(U_RATIO,2)
    I = ~isnan(U_RATIO(:,i)) & abs(U_RATIO(:,i)) < ratio_limit & abs(U_RATIO(:,i)) > 1/ratio_limit;
    U_RATIO_MEDIAN(i) = median(abs(U_RATIO(I,i)));
    I = ~isnan(U_RATIO(:,i));
    POLARITY_DISC(i) = sum(sign(U_RATIO(I,i)) > 0) / sum(I); % calculates how many theoretical and observe polarities match.
  end
  I = ~isnan(U_RATIO);
  polarity_valid_total = sum(sign(U_RATIO(I)) > 0) / sum(I(:));
  
  
  % Calculate RMS error for each event.
  for i=1:size(U_RATIO,1)
    I = ~isnan(U_RATIO(i,:));
    %RMS_ERROR(i) = sqrt( sum( (U_THEORETICAL(i,I) - U_MEASURED(i,I)).^2 ) / sum(I));
    RMS_ERROR(i) = sqrt( sum( (U_THEORETICAL(i,I) - U_MEASURED(i,I)).^2 ) ) / sqrt(sum( U_MEASURED(i,I).^2 ));
  end
  
  % Store information in history structure.
  History.RMS_ERROR(:,iteration) = RMS_ERROR;
  History.U_RATIO_MEDIAN(:,iteration) = U_RATIO_MEDIAN';
  History.U_RATIO{iteration} = U_RATIO;
  History.POLARITY_VALID(:,iteration) = POLARITY_DISC';
  History.POLARITY_VALID_TOTAL(iteration) = polarity_valid_total;
  
  % Update original amplitudes as proposed in Anderson's paper.
  WEIGHT = weight * ones(size(U_RATIO_MEDIAN));
  A = WEIGHT.*(U_RATIO_MEDIAN - 1);  % Correction factor.
  
  % Correct original displacement values from input data.
  for i=1:numel(Input)
    Sta_event = Input{i}.Station;
    U_event = Input{i}.OMEGA;
    for j=1:numel(U_event)
      k = find(strcmp(Sta_event{j},Stations));
      if ~isempty(k)
        U_event(j) = U_event(j) + A(k)*U_event(j);
      end
    end
    Input{i}.OMEGA = U_event; % Overwrite input moments with updated ones.
  end % loop for each event in particular iteration
  
  History.U = History.U + A'.*History.U; % History of moment corrections.
  History.U_CHANGES(:,iteration) = History.U;
  
  %---- Save output data.
  save(sprintf('%s/solution_%03d.mat',output_dir,iteration),'Solution');
  
  %---- Plot figures.
  if display_progress
    figure(1);
    subplot(2,3,1);
    plot(History.U_RATIO_MEDIAN');
    set(gca,'XLim',[0 n_iter+1]);
    ylabel('Amplitude ratio history');
    xlabel('iteration #');
    
    subplot(2,3,2);
    plot(History.RMS_ERROR');
    set(gca,'XLim',[0 n_iter+1]);
    ylabel('MT RMS error');
    xlabel('iteration #');
    
    subplot(2,3,3);
    imagesc(History.U_CHANGES);
    set(gca,'XLim',[0 n_iter+1]);
    ylabel('Relative changes to U');
    xlabel('iteration #');
    
    subplot(2,3,4);
    barh(POLARITY_DISC*100);
    set(gca,'YLim',[0 n_stations+1],'XLim',[0 100]);
    set(gca,'YTick',1:n_stations,'YTickLabel',Stations);
    grid on; box on;
    xlabel('Polarity match [%]');
    ylabel('Station');
    
    subplot(2,3,5);
    imagesc(log10(abs(U_RATIO)));
    set(gca,'XLim',[0 n_stations+1]);
    grid on;
    ylabel('Event #');
    xlabel('Station');
    
    subplot(2,3,6);
    barh(History.U');
    set(gca,'YLim',[0 n_stations+1]);
    set(gca,'YTick',1:n_stations,'YTickLabel',Stations);
    grid on; box on;
    xlabel('M_0 correction factor');
    ylabel('Station');
  end
  
  pause(0.1);
  
end
fprintf('\n');
disp(' hybridMT inversion finished successfully');

% Save output variables to project directory.
save(sprintf('%s/Solution.mat',output_dir),'Solution');
save(sprintf('%s/History.mat',output_dir),'History');

%---- Generate and save optimization overview.
disp('Generating summary figures');

% Create output directory for figures.
if ~exist([output_dir '/!iterations'],'dir')
  mkdir(output_dir,'!iterations');
end

% Generate figures.
f = figure('Visible','off');
plot(History.U_RATIO_MEDIAN','Color','k');
set(gca,'XLim',[0 n_iter+1],'YScale','log');
grid on; box on;
ylabel('Amplitude ratio history');
xlabel('iteration #');
saveas(gcf,sprintf('%s/!iterations/!amplitude_ratio_history.%s',output_dir,picformat));
close(f);

f = figure('Visible','off');
hold on;
plot(History.POLARITY_VALID','Color','k');
plot(History.POLARITY_VALID_TOTAL,'Color','r','LineWidth',2);
hold off;
set(gca,'XLim',[0 n_iter+1]);
grid on; box on;
ylabel('Polarity match history');
xlabel('iteration #');
saveas(gcf,sprintf('%s/!iterations/!polarity_match_history.%s',output_dir,picformat));
close(f);

f = figure('Visible','off');
plot(History.RMS_ERROR','Color','k');
set(gca,'XLim',[0 n_iter+1]);
grid on; box on;
ylabel('MT RMS error');
xlabel('iteration #');
saveas(gcf,sprintf('%s/!iterations/!rms_error_history.%s',output_dir,picformat));
close(f);

f = figure('Visible','off');
imagesc(History.U_CHANGES);
set(gca,'XLim',[0 n_iter+1]);
set(gca,'YTick',1:n_stations,'YTickLabel',Stations);
grid on; box on; colorbar;
xlabel('iteration #');
ylabel('Station');
title('Moment correction factor history');
saveas(gcf,sprintf('%s/!iterations/!amplitude_correction_history.%s',output_dir,picformat));
close(f);

f = figure('Visible','off');
barh(POLARITY_DISC*100);
set(gca,'YLim',[0 n_stations+1],'XLim',[0 100]);
set(gca,'YTick',1:n_stations,'YTickLabel',Stations);
grid on; box on;
xlabel('Polarity match [%]');
ylabel('Station');
saveas(gcf,sprintf('%s/!iterations/!polarity_match.%s',output_dir,picformat));
close(f);

f = figure('Visible','off');
barh(History.U');
set(gca,'YLim',[0 n_stations+1]);
set(gca,'YTick',1:n_stations,'YTickLabel',Stations);
grid on; box on;
xlabel('Moment correction factor');
ylabel('Station');
saveas(gcf,sprintf('%s/!iterations/!amplitude_correction.%s',output_dir,picformat));
close(f);

% Plot additional figures.
if strcmp(p.Results.AuxiliaryFigures,'on')
  hybridmt_show(output_dir,mt_solution_type,picformat);
end

disp('hybridMT finished successfully');

function movefileifexist(filein,fileout)

if exist(filein,'file')
  movefile(filein,fileout);
end

function deleteifexist(filein)

if exist(filein,'file')
  delete(filein);
end

