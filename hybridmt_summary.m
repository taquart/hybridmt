function hybridmt_summary( projectdir, mt_solution_type )
%   Copyright 2015 Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%                  Patricia Martinez-Garzon <patricia@gfz-potsdam.de>
%
%   $Revision: 1.0.0 $  $Date: 2019.03.28 $

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

fprintf('hybridMT Performance information\n');

RMS_RATIO = nan(size(Event));
for i=1:numel(Event)
  RMS_RATIO(i) = Event{i}.RMS(end) - Event{i}.RMS(1);
end
[~,I] = sort(RMS_RATIO,'descend');

j = 1;
disp('  # EIDX              ID  RMS1  RMSF  dRMS   RMS* IIDX*  ISOF/CLVDF/DCF  ISO*/CLVD*/DC*');
for i=I
  k = find(Event{i}.RMS == min(Event{i}.RMS),1,'last');
  
  fprintf('%3d %4d %15s %5.2f %5.2f ', ...
    j, i, ...
    Event{i}.event_id, ...
    Event{i}.RMS(1), ...
    Event{i}.RMS(end) ...
  );
  drms = Event{i}.RMS(end)-Event{i}.RMS(1);
  if drms > 0
    cprintf('red','%5.2f  ', drms);
  else
    cprintf('white','%5.2f  ', drms);
  end

  cprintf('white','%5.2f %5d  %4.0f/%5.0f/%3.0f  %4.0f/%5.0f/%3.0f\n', ...,
    Event{i}.RMS(k),...
    k, ...
    Event{i}.ISO(end), ...
    Event{i}.CLVD(end), ...
    Event{i}.DC(end), ...
    Event{i}.ISO(k), ...
    Event{i}.CLVD(k), ...
    Event{i}.DC(k) ...
  );
j = j + 1;
end
  


end

