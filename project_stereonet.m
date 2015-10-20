function [X, Y, R] = project_stereonet(AZM, TKO, varargin)
%PROJECT_STEREONET Project data onto stereonet plot.
%   Use PROJECT_STEREONET to project input AZIMUTH and TAKEOFF data onto the
%   stereonet. The value of AZIMUTH is measured from North towards East. The
%   value of TAKEOFF angle is measured from BOTTOM towards TOP. By default, the
%   equal area projection is used to project the points on stereonet.
%
%   part of hybridMT package 
%   <a href="matlab:open('html/doc_project_stereonet.html')">Reference page for project_stereonet</a>

%   Copyright 2015 Patricia Martinez-Garzon <patricia@gfz-potsdam.de>
%   			         Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%
%   $Revision: 1.0.2 $  $Date: 2015.09.11 $

% Parse input data.
p = inputParser;
p.addRequired('AZM', @(x) ismatrix(x) );
p.addRequired('TKO', @(x) ismatrix(x) );
p.addParamValue('Projection', 'schmidt', @(x)any(strcmpi(x,{'schmidt','wullf'}))); %#ok<*NVREPL>
p.addParamValue('Scale', 1, @(x) isscalar(x) && x > 0);
p.parse(AZM,TKO,varargin{:});

% Interpret parsed input data.
scale_factor = p.Results.Scale;
stereo_projection = p.Results.Projection;
AZM = AZM * pi/180;
TKO = TKO * pi/180;

I = TKO > pi/2;
AZM(I) = AZM(I) + pi;
TKO(I) = pi - TKO(I);
if strcmpi(stereo_projection, 'schmidt')
  R = sqrt(2)*sin(TKO/2); % Schmidt (Lambert, equal-area)
else
  R = tan(TKO/2);  % Wulff projection (Stereographic, equal-angle)
end
X = scale_factor*R.*sin(AZM);
Y = scale_factor*R.*cos(AZM);

