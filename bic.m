function h = bic(Solution, varargin)
%BIC Calculate BIC criterion for events from a particular Solution structure.
%   Copyright 2019 Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%
%   $Revision: 1.0.1 $  $Date: 2019.04.08 $

p = inputParser;
p.addRequired('Solution', @(x) iscell(x));
p.addParamValue('Sol1', 'full', @(x)any(strcmpi(x,{'full','deviatoric','dc'}))); %#ok<*NVREPL>
p.addParamValue('Sol2', 'dc', @(x)any(strcmpi(x,{'full','deviatoric','dc'}))); %#ok<*NVREPL>
p.addParamValue('Colorbar', 'RMSERROR', @(x)any(strcmpi(x,{'RMSERROR','ISO','CLVD','DC','MW'}))); %#ok<*NVREPL>
p.parse(Solution, varargin{:});

% Interpret input parameters.
sol1 = p.Results.Sol1;
sol2 = p.Results.Sol2;
cbar = p.Results.Colorbar;

if strcmpi(sol1,'full')
  n1 = 6;
elseif strcmpi(sol1,'deviatoric')
  n1 = 5;
else
  n1 = 4;
end

if strcmpi(sol2,'full')
  n2 = 6;
elseif strcmpi(sol2,'deviatoric')
  n2 = 5;
else
  n2 = 4;
end

eval(['RMS=getsolution(Solution,''full'',''' cbar ''');']);
X = getsolution(Solution,sol1,'RMSERROR');
Y = getsolution(Solution,sol2,'RMSERROR');
X_N = nan(size(X));
Y_N = nan(size(X));
for i=1:numel(X_N)
  eval(['X_N(i)=numel(Solution{i}.' sol1 '.UTH);']);
  eval(['Y_N(i)=numel(Solution{i}.' sol2 '.UTH);']);
end
BIC_X = (n1+1)*log(X_N) + X_N.*log(X./X_N);
BIC_Y = (n2+1)*log(Y_N) + Y_N.*log(Y./Y_N);

XLIM = [BIC_X; BIC_Y];
XLIM = [min(XLIM) max(XLIM)];

h = figure;
hold on;
scatter(BIC_X,BIC_Y,36,RMS,'filled','MarkerEdgeColor','k');
hold off;
xlabel(['BIC ' upper(sol1) ' moment tensors']);
ylabel(['BIC ' upper(sol2) ' moment tensors']);
box on; grid on;
axis equal;
h = colorbar; set(get(h,'title'),'string',cbar);
set(gca,'XLim',XLIM);
set(gca,'YLim',XLIM);
line(XLIM',XLIM','Color','k');
