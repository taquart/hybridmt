function [X, Y] = drawstereonet(varargin)
%DRAWSTEREONET Create stereonet plot.
%   Use DRAWSTEREONET to generate either equal area or equal-angle 
%   stereonet plot and project points onto it.
%
%   SYNTAX
%  
%   DRAWSTEREONET(...) generates stereonet plot using default values. Use 
%   'ParamName'-ParamValue arguments to adjust plot characteristics.
%
%   [X,Y] = DRAWSTEREONET(AZIMUTH, TAKEOFF, ...) projects matrices AZIMUTH
%   and TAKEOFF (of the same size) and returns X and Y matrices of the same
%   size. X and Y can be then used to plot points onto stereonet.
%
%   EXAMPLES
%
%   drawstereonet('GridColor','r','GridStepAzimuth',10) produces stereonet 
%   with red gridlines that has an azimuthal span of gridlines of 10 degrees.
%
%   drawstereonet([10 20 30 40],[50 60 70 80],'Projection','schmidt') produces
%   stereonet plot and immediately projects four points onto created stereonet
%   using Schmidt (equal-area) projection.
%
%   [X,Y] = drawstereonet([10 20 30 40],[50 60 70 80]) just projects input data
%   using default projection and returns [X,Y] vector. 

%   Copyright 2015 Patricia Martinez-Garzon <patricia@gfz-potsdam.de>
%                  Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%                  
%   $Revision: 1.0.2 $  $Date: 2015.09.11 $

% Parse input parameters.
p = inputParser;
p.addOptional('AZM', [], @(x) ismatrix(x) && isnumeric(x) );
p.addOptional('TKO', [], @(x) ismatrix(x) && isnumeric(x) );
p.addParamValue('Projection', 'schmidt', @(x)any(strcmpi(x,{'schmidt','wullf'}))); %#ok<*NVREPL>
p.addParamValue('Scale', 1, @(x) isscalar(x) && x > 0);
p.addParamValue('Grid', 'on', @(x)any(strcmpi(x,{'on','off'}))); %#ok<*NVREPL>
p.addParamValue('GridColor', 'k', @(x) ischar(x) || isvector(x));
p.addParamValue('GridStepAzimuth', 30, @(x) isscalar(x) && x > 0 && x < 360) ;
p.addParamValue('GridStepPlunge', 15, @(x) isscalar(x) && x > 0 && x < 90) ;
p.addParamValue('Location', [0 0], @(x) isvector(x) );
p.parse(varargin{:});

% Interpret input parameters.
scale_factor = p.Results.Scale;
grid_color = p.Results.GridColor;
grid_step_azimuth = p.Results.GridStepAzimuth;
grid_step_plunge = p.Results.GridStepPlunge;
projection = p.Results.Projection;
X0 = p.Results.Location(:,1);
Y0 = p.Results.Location(:,2);

AZM = p.Results.AZM;
TKO = p.Results.TKO;

% Convert points and exit if necessary.
if ~isempty(AZM) || ~isempty(TKO)
  AZM = AZM * pi/180;
  TKO = TKO * pi/180;

  I = TKO > pi/2;
  AZM(I) = AZM(I) + pi;
  TKO(I) = pi - TKO(I);
  if strcmpi(projection, 'schmidt')
    R = sqrt(2)*sin(TKO/2); % Schmidt (Lambert, equal-area)
  else
    R = tan(TKO/2);  % Wulff projection (Stereographic, equal-angle)
  end
  X = scale_factor*R.*sin(AZM);
  Y = scale_factor*R.*cos(AZM);
else
  X = [];
  Y = [];
end
if nargout > 0
  return;
end
  
% Plot stereonet and project points if available and not returned.
for i=1:length(X0)
  x0 = X0(i);
  y0 = Y0(i);
  
  a = (0:4:360)'*pi/180;
  
  % Plot background as patch object.
  Fvc.VERTICES = [scale_factor*cos(a) + x0 scale_factor*sin(a) + y0];
  Fvc.FACES = 1:size(Fvc.VERTICES,1);
  patch(Fvc,'FaceColor','w','EdgeColor','none','LineWidth',2);
  
  % Plot grid lines if necessary according to the projection.
  if strcmpi(p.Results.Grid,'on')
    for GAZM = (0:grid_step_azimuth:360)*pi/180
      GTKO = [0 90] * pi / 180;
      [GX,GY] = project_stereonet(GAZM*180/pi, GTKO*180/pi,'Projection',projection,'Scale',scale_factor);
      line(GX+ x0,GY+ y0,'Color',grid_color);
    end
    for GTKO = (0:grid_step_plunge:90)*pi/180
      GAZM = (0:2:360) * pi / 180;
      [GX,GY] = project_stereonet(GAZM*180/pi, GTKO*180/pi,'Projection',projection,'Scale',scale_factor);
      line(GX'+ x0,GY'+ y0,'Color',grid_color);
    end
  end
  
  patch(Fvc,'FaceColor','none','EdgeColor','k','LineWidth',2);
  set(gca,'Visible','off'); axis equal;
  set(gcf,'Color','w');
end

if nargout == 0 && ~isempty(X) && ~isempty(Y) 
  hold on;
  plot(X,Y,'ok','MarkerFaceColor','r');
  hold off;
end

