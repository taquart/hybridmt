function Input = readraw(filename, varargin)
%READRAW Read fociMT input ASCII file in RAW format.
%   Use READRAW(filename) to read fociMT ASCII input file in RAW format 
%   containing data for seismic moment tensor inversion.
%
%   part of hybridMT package
%   <a href="matlab:open('html/doc_readraw.html')">Reference page for readraw</a>

%   Copyright 2015 Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%   $Revision: 1.0.3 $  $Date: 2015.10.02 $

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
  if(isempty(event_id) || isempty(n))
    break;
  end
  if n < min_phases
    textscan(fid,'%s %s %s %f %f %f %f %f %f %f',n);
    continue;
  end
  if ~isempty(find(strcmp(EVENT_ID,event_id))) %#ok<EFIND>
    warning('Event with repeating ID: %s will be ignored,',event_id);
    textscan(fid,'%s %s %s %f %f %f %f %f %f %f',n);
    continue;
  else
    EVENT_ID = {EVENT_ID event_id}; %#ok<AGROW>
  end
  
  Input{j}.event_id = event_id;
  Input{j}.n_phases = n;
  Input{j}.format = 'raw';
  Input{j}.matrixmode = matrixmode;
  
  Data = textscan(fid,'%s %s %s %f %f %f %f %f %f %f',Input{j}.n_phases);
  if matrixmode
    Input{j}.Station = Data{1};
    Input{j}.Component = Data{2};
    Input{j}.Phase = Data{3};
    Input{j}.OMEGA = Data{4};
    Input{j}.AZIMUTH = Data{5};
    Input{j}.AOI = Data{6};
    Input{j}.TAKEOFF = Data{7};
    Input{j}.V = Data{8};
    Input{j}.R = Data{9};
    Input{j}.DENSITY = Data{10};
  else
    for i=1:n
      Input{j}.Phase(i).station = Data{1}{i};
      Input{j}.Phase(i).component = Data{2}{i};
      Input{j}.Phase(i).phase = Data{3}{i};
      Input{j}.Phase(i).omega = Data{4}(i);
      Input{j}.Phase(i).azimuth = Data{5}(i);
      Input{j}.Phase(i).aoi = Data{6}(i);
      Input{j}.Phase(i).takeoff = Data{7}(i);
      Input{j}.Phase(i).v = Data{8}(i);
      Input{j}.Phase(i).r = Data{9}(i);
      Input{j}.Phase(i).density = Data{10}(i);
    end
  end
  j = j + 1;
end
fclose(fid);



