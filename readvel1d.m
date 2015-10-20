function Input = readvel1d(filename, varargin)
%READVEL1D Read fociMT input ASCII file in 1D velocity model format.
%   Use READVEL1D(filename) to read ASCII fociMT input file containing 
%   data for seismic moment tensor inversion.
%
%   part of hybridMT package
%   <a href="matlab:open('html/doc_readvel1d.html')">Reference page for readvel1d</a>

%   Copyright 2015 Grzegorz Kwiatek <taquart@gmail.com>
%                  Patricia Martinez-Garzon <patricia@gfz-potsdam.de>
%
%   $Revision: 1.0.1 $  $Date: 2015.09.11 $

if nargin == 3
  matrixmode = varargin{1};
  min_phases = varargin{2};
elseif nargin == 1
  matrixmode = true;
  min_phases = 0;
else
  error('Wrong number of input parameters');
end

if ~exist(filename, 'file')
  error(['File ' filename ' does not exist.']);
end

[fid,errmsg] = fopen(filename,'r');
if fid == -1
  fclose(fid);
  error(errmsg);
end

Input = cell(1);
j = 1;
EVENT_ID = {};
while 1
  event_id = fscanf(fid,'%s',1);
  n = fscanf(fid,'%d',1);
  e_northing = fscanf(fid,'%f',1);
  e_easting = fscanf(fid,'%f',1);
  e_z = fscanf(fid,'%f',1);
  density = fscanf(fid,'%f',1);
  
  if(isempty(event_id) || isempty(n))
    break;
  end
  if n < min_phases
    textscan(fid,'%s %s %s %f %f %f %f %f %f %f',n);
    continue;
  end
  
  if ~isempty(find(strcmp(EVENT_ID,event_id)))
    warning('Following event with repeating ID: %s will be ignored,',event_id);
    textscan(fid,'%s %s %s %f %f %f %f %f %f %f',n);
    continue;
  else
    EVENT_ID = {EVENT_ID event_id}; %#ok<AGROW>
  end
  
  Input{j}.event_id = event_id;
  Input{j}.format = 'vel1d';
  Input{j}.n_phases = n;
  Input{j}.e_northing = e_northing;
  Input{j}.e_easting = e_easting;
  Input{j}.e_z = e_z;
  Input{j}.density = density;
  Input{j}.matrixmode = matrixmode;
  Data = textscan(fid,'%s %s %s  %f  %f %f %f',Input{j}.n_phases);
  if matrixmode
    Input{j}.Station = Data{1};
    Input{j}.Component = Data{2};
    Input{j}.Phase = Data{3};
    Input{j}.OMEGA = Data{4};
    Input{j}.S_NORTHING = Data{5};
    Input{j}.S_EASTING = Data{6};
    Input{j}.S_Z = Data{7};
  else
    for i=1:n
      Input{j}.Phase(i).station = Data{1}{i};
      Input{j}.Phase(i).component = Data{2}{i};
      Input{j}.Phase(i).phase = Data{3}{i};
      Input{j}.Phase(i).omega = Data{4}(i);
      Input{j}.Phase(i).s_northing = Data{5}(i);
      Input{j}.Phase(i).s_easting = Data{6}(i);
      Input{j}.Phase(i).s_z = Data{7}(i);
    end
  end
  j = j + 1;
end
fclose(fid);



