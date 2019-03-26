function hybridmt_show(projectdir, mt_solution_type,picformat)
%HYBRIDMT_SHOW Internal function to generate figures from hybridMT refinement.
%
%   part of hybridMT package
%   <a href="matlab:open('html/doc_hybridmt_show.html')">Reference page for hybridmt_show</a>

%   Copyright 2015 Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%                  Patricia Martinez-Garzon <patricia@gfz-potsdam.de>
%
%   $Revision: 1.0.5 $  $Date: 2017.01.16 $

d = dir([projectdir '/']);

n_events = 0;
n_iter = 0;
EventDirectories = cell(0);
EventID = cell(0);
SolutionFiles = cell(0);
for i=1:numel(d)
  if d(i).isdir && ~strcmp(d(i).name,'.') && ~strcmp(d(i).name,'..') && ~strcmp(d(i).name,'!iterations')
    n_events = n_events + 1;
    EventDirectories{n_events} = sprintf('%s/%s',projectdir,d(i).name);
    EventID{n_events} = d(i).name;
  end
  if ~isempty(strfind(d(i).name,'solution_')) && ~isempty(strfind(d(i).name,'.mat')) && d(i).isdir == 0
    n_iter = n_iter + 1;
    SolutionFiles{n_iter} = sprintf('%s/%s',projectdir,d(i).name);
  end
end

fprintf('Project directory: %s\n',projectdir);
fprintf('Number of events detected: %d\n',n_events);
fprintf('Number of iterations detected: %d\n',n_iter);

% Sort solution files (just in case)
SolutionFiles = sort(SolutionFiles);

% Prepare iteration matrix.
Iteration = cell(1,n_iter);

rmsg = '';
disp('Reading solution files');
for i=1:numel(SolutionFiles)
  % Display progress.
  progress = 100*i/numel(SolutionFiles);
  msg = sprintf(' File progress: %s |%s| (%1.1f%%)',SolutionFiles{i},[repmat('o',1,floor(progress/5)) repmat('-',1,ceil((100-progress)/5))],progress);
  fprintf('%s',[rmsg, msg]); rmsg = repmat(sprintf('\b'), 1, length(msg));
  
  % Load output file.
  Iteration{i} = load(SolutionFiles{i},'Solution');
  pause(0.01);
end
fprintf('\n');

% Prepare output event cell array.
Event = cell(1,n_events);
for i=1:n_events
  Event{i}.event_id = EventID{i};
  Event{i}.path = EventDirectories{i};
  Event{i}.solution_type = mt_solution_type;
  Event{i}.ISO = nan(n_iter,1);
  Event{i}.CLVD = nan(n_iter,1);
  Event{i}.DC = nan(n_iter,1);
  Event{i}.M0 = nan(n_iter,1);
  Event{i}.MT = nan(n_iter,1);
  Event{i}.M0ERRMAX = nan(n_iter,1);
  Event{i}.RMS = nan(n_iter,1);
  Event{i}.MW = nan(n_iter,1);
  Event{i}.P = nan(n_iter,2);
  Event{i}.T = nan(n_iter,2);
  Event{i}.B = nan(n_iter,2);
  Event{i}.F1 = nan(n_iter,3);
  Event{i}.F2 = nan(n_iter,3);
  Event{i}.MXX = nan(n_iter,6);
end

%---- Loop through iterations and gather relevant information.

for i=1:n_iter
  
  % Loop through events in particular iteration (should be always the same
  % number of events therein)
  
  if numel(Iteration{i}.Solution) ~= n_events
    error('Not enough events in iteration %d. Procedure aborted.',i);
  end
  
  for e=1:numel(Iteration{i}.Solution)
    event_id = Iteration{i}.Solution{e}.event_id;
    
    switch mt_solution_type
      case 'full'
        Solution = Iteration{i}.Solution{e}.full;
      case 'deviatoric'
        Solution = Iteration{i}.Solution{e}.deviatoric;
      case 'dc'
        Solution = Iteration{i}.Solution{e}.dc;
    end
    
    j = find(strcmpi(EventID, event_id)); % find index of a particular event.
    if isempty(j)
      error('Cannot find specific event_id ("%s") in the directory list. Procedure aborted.',event_id);
    end
    
    % i - iteration index (1...n_iter)
    % j - event index (1...n_events)
    Event{j}.ISO(i,:) = [Solution.ISO ];
    Event{j}.CLVD(i,:) = [Solution.CLVD ];
    Event{j}.DC(i,:) = [Solution.DC];
    Event{j}.F1(i,:) = [Solution.F1 ];
    Event{j}.F2(i,:) = [Solution.F2 ];
    Event{j}.P(i,:) = [Solution.P ];
    Event{j}.T(i,:) = [Solution.T ];
    Event{j}.B(i,:) = [Solution.B ];
    Event{j}.M0(i,:) = [Solution.M0 ];
    Event{j}.MT(i,:) = [Solution.MT ];
    Event{j}.M0ERRMAX(i,:) = [Solution.M0ERRMAX ];
    Event{j}.RMS(i,:) = [Solution.RMSERROR ];
    Event{j}.MW(i,:) = [Solution.MW ];
    Event{j}.MXX(i,:) = [Solution.MXX ];
    
  end
end

%---- Generate figures for each iteration.

if ~exist([projectdir '/!iterations'],'dir')
  mkdir(projectdir,'!iterations');
end

rmsg = '';
disp('Generating summary figures (iterations)');
for i=1:n_iter
  try
    titletext = sprintf('Iteration: %d',i);
    
    progress = 100*i/n_iter;
    msg = sprintf(' Figure progress: %3d/%3d |%s| (%1.1f%%)',i,n_iter,[repmat('o',1,floor(progress/5)) repmat('-',1,ceil((100-progress)/5))],progress);
    fprintf('%s',[rmsg, msg]); rmsg = repmat(sprintf('\b'), 1, length(msg));
    
    P = nan(n_events,2);
    T = nan(n_events,2);
    B = nan(n_events,2);
    MXX = nan(n_events,6);
    for j=1:n_events
      P(j,:) = Event{j}.P(i,:);
      T(j,:) = Event{j}.T(i,:);
      B(j,:) = Event{j}.B(i,:);
      MXX(j,:) = Event{j}.MXX(i,:);
    end
    
    % Plot P/T/B axes.
    f = figure('Visible','off');
    hold on;
    drawstereonet('Projection','wullf');
    [X,Y] = drawstereonet(P(:,1), 90 - P(:,2),'Projection','wullf');
    plot(X,Y,'ko','Marker','o','MarkerSize',5,'MarkerFaceColor','b');
    [X,Y] = drawstereonet(T(:,1), 90 - T(:,2),'Projection','wullf');
    plot(X,Y,'ko','Marker','o','MarkerSize',5,'MarkerFaceColor','r');
    [X,Y] = drawstereonet(B(:,1), 90 - B(:,2),'Projection','wullf');
    plot(X,Y,'ko','Marker','o','MarkerSize',5,'MarkerFaceColor',[0.7 0.7 0.7]);
    hold off;
    axis equal;
    set(gca,'Visible','off');
    text(0,1.05,titletext,'HorizontalAlignment','center','Interpreter','none','VerticalAlignment','bottom');
    saveas(gcf,sprintf('%s/!iterations/ptb-stereonet-%03d.%s',projectdir,i,picformat));
    close(f);
    
    % Hudson plot.
    f = figure('Visible','off');
    [U,V] = drawhudsonnet(MXX);
    hold on;
    plot(U,V,'ok','MarkerSize',5,'MarkerFaceColor','r');
    hold off;
    set(gca,'Visible','off');
    saveas(gcf,sprintf('%s/!iterations/hudson-stereonet-%03d.%s',projectdir,i,picformat));
    close(f);
    
    pause(0.01);
  catch
    warning('hybridmt_show:unexpected_exception','Unexpected exception while preparing the figure plots');
  end
end
fprintf('\n');


%---- Generate figures for each event.
pause(1.0);

rmsg = '';
disp('Generating summary figures (events)');
for i=1:n_events
  
  titletext = sprintf('%s (%s MT solution)',Event{i}.event_id,Event{i}.solution_type);
  
  % Progress bar.
  progress = 100*i/n_events;
  msg = sprintf(' Figure progress: %3d/%3d |%s| (%1.1f%%)',i,n_events,[repmat('o',1,floor(progress/5)) repmat('-',1,ceil((100-progress)/5))],progress);
  fprintf('%s',[rmsg, msg]); rmsg = repmat(sprintf('\b'), 1, length(msg));
  
  % Moment tensor decomposition changes.
  f = figure('Visible','off');
  hold on;
  plot(Event{i}.ISO,'r-','Marker','.');
  plot(Event{i}.CLVD,'g-','Marker','.');
  plot(Event{i}.DC,'b-','Marker','.');
  hold off;
  grid on; box on;
  set(gca,'YLim',[-100 100]);
  legend('Isotropic','CLVD','DC');
  xlabel('Iteration');
  ylabel('% Moment tensor components');
  title(titletext,'Interpreter','none');
  saveas(gcf,[Event{i}.path '/' Event{i}.event_id '-decomposition.' picformat]);
  close(f);
  
  % P/T/B axes directions.
  f = figure('Visible','off');
  subplot(2,1,1);
  hold on;
  plot(Event{i}.P(:,1),'b-','Marker','.');
  plot(Event{i}.T(:,1),'r-','Marker','.');
  plot(Event{i}.B(:,1),'k-','Marker','.');
  hold off;
  grid on; box on;
  set(gca,'YLim',[0 360]);
  legend('P trend','T trend','B trend');
  ylabel('P/T/B Trend');
  title(titletext,'Interpreter','none');
  subplot(2,1,2);
  hold on;
  plot(Event{i}.P(:,2),'b-','Marker','.');
  plot(Event{i}.T(:,2),'r-','Marker','.');
  plot(Event{i}.B(:,2),'k-','Marker','.');
  hold off;
  grid on; box on;
  set(gca,'YLim',[0 90]);
  legend('P plunge','T plunge','B plunge');
  xlabel('Iteration');
  ylabel('P/T/B Plunge');
  saveas(gcf,[Event{i}.path '/' Event{i}.event_id '-ptb_axes.' picformat]);
  close(f);
  
  % Error dropdown.
  f = figure('Visible','off');
  hold on;
  plot(Event{i}.RMS,Event{i}.M0ERRMAX,'k-','Marker','.');
  plot(Event{i}.RMS(end),Event{i}.M0ERRMAX(end),'ko');
  hold off;
  grid on; box on;
  set(gca,'YScale','log');
  xlabel('RMS Error');
  ylabel('Maximum error of MT component [Nm]');
  title(titletext,'Interpreter','none');
  saveas(gcf,[Event{i}.path '/' Event{i}.event_id '-errors.' picformat]);
  close(f);
  
  % Seismic moment information.
  f = figure('Visible','off');
  hold on;
  plot(Event{i}.M0,'k-','Marker','.');
  plot(Event{i}.MT,'r-','Marker','.');
  plot(Event{i}.M0ERRMAX,'b-','Marker','.');
  hold off;
  grid on; box on;
  set(gca,'YScale','log');
  xlabel('Iteration');
  ylabel('Seismic moment [Nm]');
  title(titletext,'Interpreter','none');
  saveas(gcf,[Event{i}.path '/' Event{i}.event_id '-moment.' picformat]);
  close(f);
  
  % P/T/B axis directions on stereonet plot.
  f = figure('Visible','off');
  hold on;
  drawstereonet('Projection','wullf');
  [X,Y] = drawstereonet(Event{i}.P(:,1), 90 - Event{i}.P(:,2),'Projection','wullf');
  plot(X,Y,'b-','Marker','.');
  plot(X(end),Y(end),'ko','Marker','o','MarkerSize',5,'MarkerFaceColor','r');
  text(X(1),Y(1),'P');
  [X,Y] = drawstereonet(Event{i}.T(:,1), 90 - Event{i}.T(:,2),'Projection','wullf');
  plot(X,Y,'r-','Marker','.');
  plot(X(end),Y(end),'ko','Marker','o','MarkerSize',5,'MarkerFaceColor','r');
  text(X(1),Y(1),'T');
  [X,Y] = drawstereonet(Event{i}.B(:,1), 90 - Event{i}.B(:,2),'Projection','wullf');
  plot(X,Y,'k-','Marker','.');
  plot(X(end),Y(end),'ko','Marker','o','MarkerSize',5,'MarkerFaceColor','k');
  text(X(1),Y(1),'B');
  hold off;
  axis equal;
  set(gca,'Visible','off');
  text(0,1.05,titletext,'HorizontalAlignment','center','Interpreter','none','VerticalAlignment','bottom');
  saveas(gcf,[Event{i}.path '/' Event{i}.event_id '-ptb_axes-stereonet.' picformat]);
  close(f);
  
  % Hudson plot of moment tensor decomposition.
  try
    f = figure('Visible','off');
    [U,V] = drawhudsonnet(Event{i}.MXX);
    hold on;
    plot(U,V,'-k','Marker','.');
    plot(U(end),V(end),'ok','MarkerSize',5,'MarkerFaceColor','r');
    hold off;
    set(gca,'Visible','off');
    saveas(gcf,[Event{i}.path '/' Event{i}.event_id '-hudson_plot.' picformat]);
    close(f);
  catch
    warning('hybridmt_show:unexpected_exception','Unexpected exception while preparing Hudson''s plot');
  end
  
  pause(0.01);
end
fprintf('\n');
