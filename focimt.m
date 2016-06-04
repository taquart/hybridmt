function [Solution, Input, Params] = focimt(INPUT, varargin)
%FOCIMT Perform seismic full moment tensor inversion using fociMT software.
%   Use FOCIMT to perform the seismic moment tensor inversion using fociMT
%   software via MATLAB wrapper.
%
%   part of hybridMT package
%   <a href="matlab:open('html/doc_focimt.html')">Reference page for focimt</a>

%   Copyright 2015-2016 Grzegorz Kwiatek <taquart@gmail.com>
%                       Patricia Martinez-Garzon <patricia@gfz-potsdam.de>
%
%   $Revision: 1.0.12 $  $Date: 2016.06.04 $

% Parse input parameters.
p = inputParser;
p.addRequired('INPUT', @(x) ischar(x) || iscell(x));
p.addParamValue('Jacknife', 'off', @(x)any(strcmpi(x,{'on','off'}))); %#ok<*NVREPL>
p.addParamValue('Verbose', 'off', @(x)any(strcmpi(x,{'on','off'}))); %#ok<*NVREPL>
p.addParamValue('Norm', 'L2', @(x)any(strcmpi(x,{'L1','L2'})));
p.addParamValue('BeachBallSize', 300, @(x) isscalar(x) && x > 0);
p.addParamValue('BeachBallFormat', 'NONE', @(x)any(strcmpi(x,{'PNG','PS','NONE','SVG','PDF'})));
p.addParamValue('Projection', 'schmidt', @(x)any(strcmpi(x,{'wullf','schmidt'})));
p.addParamValue('Hemisphere', 'lower', @(x)any(strcmpi(x,{'lower','upper'})));
p.addParamValue('IgnoreStation', cell(0), @(x) iscell(x) || ischar(x));
p.addParamValue('VelocityModel', [], @(x) isnumeric(x) && size(x,2) == 2);
p.addParamValue('MinimumPhases', 8, @(x) isscalar(x) && x > 6);
p.addParamValue('Bootstrap', [], @(x ) all(size(x) == [1 2]) || all(size(x) == [1 3]) || all(size(x) == [1 4]));
p.addParamValue('CorrectStation', cell(0), @(x) iscell(x));
p.addParamValue('ProjectDir', '', @(x) ischar(x) );
p.addParamValue('PlotCross', 'on', @(x)any(strcmpi(x,{'on','off'})));
p.addParamValue('PlotStations', 'on', @(x)any(strcmpi(x,{'on','off'})));
p.addParamValue('PlotAxes', 'on', @(x)any(strcmpi(x,{'on','off'})));
p.addParamValue('PlotDC', 'on', @(x)any(strcmpi(x,{'on','off'})));
%
p.addParamValue('Solutions', 'FTD', @(x) ischar(x));
p.addParamValue('Decomposition', 'JostHerrmann', @(x)any(strcmpi(x,{'JostHerrmann','Vavrycuk'})));
p.addParamValue('NormalFaultColor', [], @(x)all(size(x) == [1 3]) || all(size(x) == [1 4]));
p.addParamValue('StrikeSlipFaultColor', [], @(x)all(size(x) == [1 3]) || all(size(x) == [1 4]));
p.addParamValue('ThrustFaultColor', [], @(x)all(size(x) == [1 3]) || all(size(x) == [1 4]));
p.addParamValue('DoubleCoupleColor', [], @(x)all(size(x) == [1 3]) || all(size(x) == [1 4]));
p.addParamValue('TShadingColor', [], @(x)all(size(x) == [1 3]) || all(size(x) == [1 4]));
p.addParamValue('PShadingColor', [], @(x)all(size(x) == [1 3]) || all(size(x) == [1 4]));
p.addParamValue('PlusColor', [], @(x)all(size(x) == [1 3]) || all(size(x) == [1 4]));
p.addParamValue('MinusColor', [], @(x)all(size(x) == [1 3]) || all(size(x) == [1 4]));
p.addParamValue('LabelColor', [], @(x)all(size(x) == [1 3]) || all(size(x) == [1 4]));

p.parse(INPUT,varargin{:});

% Interpret input parameters.
Params = p.Results;
ball = '-b ';
if strcmpi(p.Results.PlotStations,'on')
  ball = [ball 'S'];
end
if strcmpi(p.Results.PlotAxes,'on')
  ball = [ball 'A'];
end
if strcmpi(p.Results.PlotCross,'on')
  ball = [ball 'C'];
end
if strcmpi(p.Results.PlotDC,'on')
  ball = [ball 'D'];
end

% Interpratation of colors
colors_text = '';
colorvars = {'NormalFaultColor','StrikeSlipFaultColor','ThrustFaultColor', ...
  'DoubleCoupleColor','TShadingColor','PShadingColor','PlusColor','MinusColor',...
  'LabelColor'};
colorpars = {'n','s','r','d','t','p','+','-','l'};
for i=1:numel(colorvars)
  if ~isempty(eval(['p.Results.' colorvars{i}]))
    colors_text = [colors_text prncolor(colorpars{i},eval(['p.Results.' colorvars{i}]))]; %#ok<AGROW>
  end
end

if strcmpi(p.Results.Decomposition,'JostHerrmann')
  decomposition = 'D';
else
  decomposition = 'Y';
end

solutions = p.Results.Solutions;

jacknife = '';
if strcmpi(p.Results.Jacknife,'on')
  jacknife = '-j';
end
bootstrap = '';
if ~isempty(p.Results.Bootstrap)
  B = p.Results.Bootstrap;
  if numel(B) >= 2
    bootstrap = ['-rp ' num2str(B(1)) '/' num2str(B(2))];
  end
  if numel(B) >= 3
    bootstrap = [bootstrap ' -rr ' num2str(B(1)) '/' num2str(B(3))];
  end
  if numel(B) >= 4
    bootstrap = [bootstrap ' -ra ' num2str(B(1)) '/' num2str(B(4))];
  end
end
normfunc = p.Results.Norm;
bbsize = ['-z ' num2str(p.Results.BeachBallSize)];
bbformat = ['-t ' p.Results.BeachBallFormat];
if strcmpi(p.Results.Projection,'wullf')
  bbprojection = 'W';
else
  bbprojection = 'S';
end
if strcmpi(p.Results.Hemisphere,'lower')
  bbhemisphere = 'L';
else
  bbhemisphere = 'U';
end
bbproj = ['-p ' bbprojection bbhemisphere];
verbose = strcmpi(p.Results.Verbose,'on');
StIgnored = p.Results.IgnoreStation;
if ischar(StIgnored)
  StIgnored = {StIgnored};
end
VMODEL = p.Results.VelocityModel;
CoStation = p.Results.CorrectStation;
min_phases = p.Results.MinimumPhases;
projectdir = p.Results.ProjectDir;

if verbose
  tic;
end


% Check if project directory exists.
if ~isempty(projectdir)
  if exist(projectdir,'dir')
    error(['Project directory ' projectdir ' already exist.']);
  else
    mkdir(projectdir);
  end
end

% Load input file.
if verbose, disp('Preparing input file'); end
filename = '';
if ischar(INPUT)
  if exist(INPUT,'file')
    filename = INPUT;
    % Determine if the file is in raw or velocity 1d format.
    fid = fopen(filename);
    header = fgetl(fid);
    fclose(fid);
    EDATA = sscanf(header,'%*s %f %f %f %f %f');
    if numel(EDATA) == 5
      Input = readvel1d(filename, true, min_phases);
    else
      Input = readraw(filename, true, min_phases);
    end
  else
    error('File %s does not exist.',INPUT);
  end
elseif iscell(INPUT);
  Input = INPUT;
  clear INPUT;
else
  error('Input data should be either a valid file name of or a cell array in fociMT format.');
end

% Check if ID field is repeating...
IDNames = cell(numel(Input),1);
for i=1:numel(Input)
  IDNames{i} = Input{i}.event_id;
end
if numel(unique(IDNames)) ~= numel(IDNames)
  error('Input file contains repetitive ID numbers.');
end

% Remove unnecessary stations on request
if numel(StIgnored) > 0
  for i=1:numel(StIgnored)
    station_ignored = StIgnored{i};
    for j=1:numel(Input)
      I = ~strcmpi(Input{j}.Station,station_ignored);
      Input{j}.n_phases = sum(I);
      Input{j}.Station = Input{j}.Station(I);
      Input{j}.Component = Input{j}.Component(I);
      Input{j}.Phase = Input{j}.Phase(I);
      Input{j}.OMEGA = Input{j}.OMEGA(I);
      if strcmpi(Input{j}.format, 'raw')
        Input{j}.AZIMUTH = Input{j}.AZIMUTH(I);
        Input{j}.AOI = Input{j}.AOI(I);
        Input{j}.TAKEOFF = Input{j}.TAKEOFF(I);
        Input{j}.V = Input{j}.V(I);
        Input{j}.R = Input{j}.R(I);
        Input{j}.DENSITY = Input{j}.DENSITY(I);
      else
        Input{j}.S_NORTHING = Input{j}.S_NORTHING(I);
        Input{j}.S_EASTING = Input{j}.S_EASTING(I);
        Input{j}.S_Z = Input{j}.S_Z(I);
      end
    end
  end
end

% Correct input amplitudes on request.
if ~isempty(CoStation)
  for i=1:2:numel(CoStation)
    station = CoStation{i};
    factor = CoStation{i+1};
    for j=1:numel(Input)
      I = strcmpi(Input{j}.Station,station);
      Input{j}.OMEGA(I) = Input{j}.OMEGA(I) * factor;
    end
  end
end

% Prepare input data for moment tensor inversion.
if isempty(projectdir)
  temp = ['./' char(randi(26,1,8)+96)];
else
  temp = ['./' projectdir '/' char(randi(26,1,8)+96)];
end
focimt_deletefile({[temp '.txt'],[temp '-full.asc'],[temp '-deviatoric.asc'],[temp '-dc.asc'], ...
  [temp '-full-u.asc'],[temp '-deviatoric-u.asc'],[temp '-dc-u.asc'],[temp '_vmodel.txt']});

% If velocity model is used, create velocity model file for fociMT
% executable.
if strcmpi(Input{1}.format,'vel1d')
  if isempty(VMODEL)
    error(['Input file is in 1D velocity model format, but the velocity model ' ...
      'has not been specified. Use ''VelocityModel'' property to specify the velocity model matrix']);
  end
  fid = fopen([temp '_vmodel.txt'],'w');
  fprintf(fid,'%d\r\n',size(VMODEL,1));
  fprintf(fid,'%1.2f\r\n',VMODEL(:,1));
  fprintf(fid,'%1.2f\r\n',VMODEL(:,2));
  fclose(fid);
  vmodel = ['-m ' temp '_vmodel.txt'];
else
  vmodel = '';
end

commandline = ['-i ' temp '.txt ' colors_text ' -d ' decomposition 'WAFTUMVE -o ' temp ' -s ' solutions ' -n ' normfunc ' ' bbsize ' ' bbformat ' ' bbproj ' ' jacknife ' ' bootstrap ' ' vmodel ' ' ball];

% Prepare input file for focimt application.
writeinput([temp '.txt'], Input);

if verbose
  disp(['Executing fociMT (' commandline ')']);
end

% Execute focimt program.
focimt_path = mfilename('fullpath');
focimt_path = fileparts(focimt_path);
status = system(['"' focimt_path '/focimt.exe" ' commandline]);
if status ~= 0
  warning('FOCIMT:exec_status','fociMT binary exited with status: %d', status);
elseif verbose
  fprintf('fociMT finished with exit code: %d\n',status);
end

% Read output files.
Solution = cell(1);
if ~isempty(strfind(solutions,'F'))
  if verbose, disp('Reading full MT solution output file'); end
  Solution = readsolution([temp '-full'], 'full', Solution,true);
end
if ~isempty(strfind(solutions,'T'))
  if verbose, disp('Reading deviatoric MT solution output file'); end
  Solution = readsolution([temp '-deviatoric'], 'deviatoric', Solution,true);
end
if ~isempty(strfind(solutions,'D'))
  if verbose, disp('Reading double-couple MT solution output file'); end
  Solution = readsolution([temp '-dc'], 'dc', Solution,true);
end

% Delete unnecessary files.
if verbose, disp('Deleting unnecessary files'); end
focimt_deletefile({[temp '.txt'],[temp '-full.asc'],[temp '-deviatoric.asc'],[temp '-dc.asc'], ...
  [temp '-full-u.asc'],[temp '-deviatoric-u.asc'],[temp '-dc-u.asc'],[temp '_vmodel.txt']});

if exist(filename,'file') && exist(projectdir,'dir');
  copyfile(filename,projectdir);
end

% Save input parameters
if verbose, disp('Saving inversion parameters to Params.mat'); end
if exist(projectdir,'dir')
  save(['./' projectdir '/Params.mat'],'Params');
  save(['./' projectdir '/Solution.mat'],'Solution');
else
  save('Params.mat','Params');
  save('Solution.mat','Solution');
end

if verbose
  toc;
end

%==============================================================================
%==============================================================================
%==============================================================================
function focimt_deletefile(Files)

for i=1:numel(Files)
  filename = Files{i};
  if exist(filename,'file');
    delete(filename);
  end
  if exist(filename,'file');
    warning('FOCIMT:file_locked',['Cannot delete ' filename '. Executing "fclose all" command']);
    fclose all;
    delete(filename);
  end
  if exist(filename,'file');
    error(['Cannot delete file (' filename '. File locked outside of MATLAB?']);
  end
end

function colorvector = prncolor(ctype,colorvector)

if all(size(colorvector) == [1 3])
  colorvector = sprintf('-c%s %1.3f/%1.3f/%1.3f ',ctype,colorvector);
elseif all(size(colorvector) == [1 4])
  colorvector = sprintf('-c%s %1.3f/%1.3f/%1.3f/%1.3f ',ctype,colorvector);
else
  error('Color vector must be either of [1,3] or [1,4] size [R/G/B[/A]]');
end


