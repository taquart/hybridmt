function genmt_raw(STRIKE, DIP, RAKE, varargin)
%GETMT_RAW Generate synthetic event/phase data for fociMT in RAW ASCII format.
%	Use GETMT_RAW to generate a synthetic amplitude data for the purpose of 
%   seismic moment inversion. The point source is assumed and the radiation 
%   formula of shear-tensile source model are used to calculate the expected
%   amplitudes at specified AZIMUTHs and TAKEOFF angles from the source.
%
%   SYNTAX
%
%   GENMT_RAW(STRIKE,DIP,RAKE) generates synthetic event/phase data for 
%   double-couple event(s) with focal mechanism specified by STRIKE, DIP and
%   RAKE vectors (the vectors must be of the same size and the angles are 
%   specified in degrees in standard seismological convention).
%
%   GENMT_RAW(STRIKE,DIP,RAKE,ANGLE) allows to specify the tensile angle for
%   every fault plane. The negative tensile angles correspond to compaction (with 
%   -90deg corresponding to MODE-I closing) whereas posivitive angles correspond 
%   to tensile opening (with +90deg corresponds to MODE-I opening). The value 0deg
%   corresponds to pure shear motion. 
%
%   GENMT_RAW(STRIKE, DIP, RAKE, ...) allows to specify additional parameters.
%
%   PARAMETERS
% 
%   'Distance'
%
%   part of hybridMT package
%   <a href="matlab:open('html/doc_genmt_raw.html')">Reference page for genmt_raw</a>

%   Copyright 2015 Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%                  Patricia Martinez-Garzon <patricia@gfz-potsdam.de>
%
%   $Revision: 1.0.2 $  $Date: 2015.09.11 $

%---- Parse input parameters.
p = inputParser;
p.addRequired('STRIKE', @(x) isvector(x) && all(x>=0 & x<360.0));
p.addRequired('DIP', @(x) isvector(x) && all(x>=0 & x<=90));
p.addRequired('RAKE', @(x) isvector(x) && all(x>=-180 & x<=180.0));
p.addOptional('ANGLE', nan, @(x) isvector(x) && all(x>=-90 & x<=90.0));
p.addParamValue('Distance', 3000, @(x) isscalar(x) && x > 0); %#ok<*NVREPL>
p.addParamValue('Velocity', 5000, @(x) isscalar(x) && x > 0);
p.addParamValue('StationBias', [], @(x) ismatrix(x));
p.addParamValue('Density', 2700, @(x) isscalar(x) && x > 0);
p.addParamValue('PoissonRatio', 0.25, @(x) isscalar(x) && x > 0);
p.addParamValue('MomentMagnitude', pi, @(x) isscalar(x));
p.addParamValue('Takeoff', [20 150 50 100], @(x) isvector(x) && all(x >= 0 & x <= 180));
p.addParamValue('Azimuth', [20 80 140 200 260 320], @(x) isvector(x) && all(x >= 0 & x < 360));
p.addParamValue('FileName', 'default.txt', @(x) ischar(x) );
p.addParamValue('PicksLost',0, @(x) isscalar(x) && x >= 0.0 && x < 1.0);

p.parse(STRIKE, DIP, RAKE, varargin{:});

% Assume basic parameters.
filename = p.Results.FileName;
if isempty(filename)
  commonfile = false;
else
  commonfile = true;
  if exist(filename,'file');
    delete(filename);
  end
end
dist = p.Results.Distance;
vel = p.Results.Velocity;
den = p.Results.Density;
M0 = 10.^(3*(p.Results.MomentMagnitude+6.07)/2); 
TAKEOFF = p.Results.Takeoff;
AZM     = p.Results.Azimuth;
ANGLE = p.Results.ANGLE;
if isnan(ANGLE)
  ANGLE = zeros(size(STRIKE));
elseif isscalar(ANGLE)
  ANGLE = ANGLE * ones(size(STRIKE));
end
sigma = p.Results.PoissonRatio;
picklost = p.Results.PicksLost;
BIAS = p.Results.StationBias;

for j=1:length(STRIKE)
  
  % Generate stations around the focal sphere.
  [TAKEOFFM,AZMM] = meshgrid(TAKEOFF,AZM);
  AOIM = 180.0 - TAKEOFFM;
  
  % Produce radiation pattern (pure shear faulting)
  RP = rpgen(STRIKE(j), DIP(j), RAKE(j), ANGLE(j), sigma, TAKEOFFM, AZMM);
  
  % Modify amplitude as it would be recorded at vertical sensor.
  RP_ZZ = RP .* cos(AOIM * pi / 180);
  
  % Rescale amplitude to seismic moment.
  AMPLITUDE = RP_ZZ ./ 4 / pi / den / vel ^ 3 / dist * M0;
  
  % Plotting over focal sphere.
  %   [X,Y,Z] = sph2cart(pi / 2 - AZMM * pi/180, pi / 2 - AOIM * pi/180, dist);
  %   hold on;
  %   i = find(RP_ZZ >=0 ); plot3(X(i),Y(i),Z(i),'or');
  %   i = find(RP_ZZ < 0 ); plot3(X(i),Y(i),Z(i),'ob');
  %   hold off;
  %   for i=1:length(X(:))
  %     text(X(i),Y(i),Z(i),num2str(i),'FontSize',8,'HorizontalAlignment','center','Verticalalignment','bottom');
  %     text(X(i),Y(i),Z(i),['I' num2str(AOIM(i)) ' T' num2str(TAKEOFFM(i)) ' A' num2str(AZMM(i))],'FontSize',8,'HorizontalAlignment','center','Verticalalignment','top');
  %   end
  %   axis vis3d;
  %   xlabel('x [m]');
  %   ylabel('y [m]');
  %   zlabel('z [m]');
  %   box on;
  %   grid on;
  
  % Output input file for foci-mt in standard format.
  if commonfile
    % There will be one file which contains all input data.
    fileid = sprintf('%03d_%03d_%+03d',STRIKE(j),DIP(j),RAKE(j));
  else
    % There will be one file generated pear strike/dip/rake triplet.
    filename = sprintf('%03d_%03d_%+03d.txt',STRIKE(j),DIP(j),RAKE(j));
    fileid = sprintf('%03d_%03d_%+03d',STRIKE(j),DIP(j),RAKE(j));
    if exist(filename, 'file')
      delete(filename);
    end
  end
  
  n_stations = length(TAKEOFFM(:));
  n_lost = floor(picklost * n_stations);
  ACC = ones(size(TAKEOFFM(:)));
  while n_lost > 0
    k = randi(n_stations);
    if ACC(k) == 1
      ACC(k) = 0;
      n_lost = n_lost - 1;
    end
  end
  
  % Actual number of stations.
  fid = fopen(filename,'a');
  fprintf(fid,'%s %d\r\n',fileid,sum(ACC));
  for i = 1:n_stations
    if(ACC(i))
      for k=1:size(BIAS,1)
        if i == BIAS(k,1)
          AMPLITUDE(i) = AMPLITUDE(i) * BIAS(k,2);
        end
      end
      fprintf(fid,'S%02d Z P %22.14g  %6.1f %6.1f %6.1f %5.0f %7.2f %5.0f\r\n', ...
        i, AMPLITUDE(i), AZMM(i), AOIM(i), TAKEOFFM(i), vel, dist, den);
    end
  end
  fclose(fid);
  
end
