function varargout = getsolution(Res,solution,varargin)
%GETSOLUTION Extract parameter(s) from moment tensor solution cell array
%   Use GETSOLUTION to convenietly extract seismic moment tensor solution
%   results from multiple seismic events.
%
%   part of hybridMT package
%   <a href="matlab:open('html/doc_getsolution.html')">Reference page for getsolution</a>

%   Copyright 2015 Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%                  Patricia Martinez-Garzon <patricia@gfz-potsdam.de>
%
%   $Revision: 1.0.3 $  $Date: 2015.09.16 $

% Count first the number of elements
for i=1:numel(varargin)
  varname = varargin{i};
  n = 0;
  for j=1:numel(Res)
    if isfield(eval(['Res{j}.' solution]),varname)
      n = n + size(eval(['Res{j}.' solution '.' varname]),1);
    else
      n = n + size(eval(['Res{j}.' varname]),1);
    end
  end
  
  if isfield(eval(['Res{1}.' solution]),varname)
    m = size(eval(['Res{1}.' solution '.' varname]),2);
    if iscell(eval(['Res{1}.' solution '.' varname]))
      TEMP = cell(n,1);
      cella = true;
    else
      TEMP = nan(n,m);
      cella = false;
    end
  else
    m = size(eval(['Res{1}.' varname]),2);
    if iscell(eval(['Res{1}.' varname])) || ischar(eval(['Res{1}.' varname]))
      TEMP = cell(n,1);
      cella = true;
    else
      TEMP = nan(n,m);
      cella = false;
    end
  end
  
  n = 0;
  for j=1:numel(Res)
    if isfield(eval(['Res{j}.' solution]),varname)
      len = size(eval(['Res{j}.' solution '.' varname]),1);
      if cella
        TEMP{n+1} = eval(['Res{j}.' solution '.' varname]);
      else
        TEMP(n+1,:) = eval(['Res{j}.' solution '.' varname]);
      end
    else
      len = size(eval(['Res{j}.' varname]),1);
      if cella
        TEMP{n+1} = eval(['Res{j}.' varname]);
      else
        TEMP(n+1,:) = eval(['Res{j}.' varname]);
      end
    end
    n = n + len;
  end
  varargout{i} = TEMP;
end


