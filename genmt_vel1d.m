function genmt_vel1d(STRIKE, DIP, RAKE, varargin)
%GENMT_VEL1D Generate synthetic event/phase data for fociMT in 1D velocity model ASCII format
%
%   part of hybridMT package
%   <a href="matlab:open('html/doc_genmt_vel1d.html')">Reference page for genmt_vel1d</a>

%   Copyright 2015 Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%                  Patricia Martinez-Garzon <patricia@gfz-potsdam.de>
%
%   $Revision: 1.0.1 $  $Date: 2015.06.01 $

%---- Parse input parameters.
p = inputParser;
p.addRequired('STRIKE', @(x) isvector(x) && all(x>=0 & x<360.0));
p.addRequired('DIP', @(x) isvector(x) && all(x>=0 & x<=90));
p.addRequired('RAKE', @(x) isvector(x) && all(x>=-180 & x<=180.0));
p.addOptional('ANGLE', nan, @(x) isvector(x) && all(x>=-90 & x<=90.0));
p.addParamValue('PoissonRatio', 0.25, @(x) isscalar(x) && x > 0); %#ok<*NVREPL>
p.addParamValue('Density', 2700, @(x) isscalar(x) && x > 0);
p.addParamValue('MomentMagnitude', pi, @(x) isscalar(x));
p.addParamValue('FileFormat', 'vel1D', @(x)any(strcmpi(x,{'raw','vel1D'})));
p.addParamValue('StationBias', [], @(x) ismatrix(x));
p.addParamValue('FileName', '', @(x) ischar(x) );
p.addParamValue('PicksLost',0, @(x) isscalar(x) && x >= 0.0 && x < 1.0);
p.addParamValue('Northing', meshgrid(-10000:5000:10000), @(x) ismatrix(x) );
p.addParamValue('Easting', meshgrid(-10000:5000:10000)', @(x) ismatrix(x) );
p.addParamValue('EventDepth', 1500, @(x) isvector(x) && all(x >= 0 & x < 6371000));
p.addParamValue('VelocityModel', [0.00 3.00 8.00 20.00 22.00; 4.10 5.47 5.75 6.02 7.90]', @(x) ismatrix(x));
p.parse(STRIKE, DIP, RAKE, varargin{:});

% Assume basic parameters.
filename = p.Results.FileName;
if isempty(filename)
  commonfile = false;
else
  commonfile = true;
  if exist(filename,'file')
    delete(filename);
  end
end

% Velocity model:
VMODEL = p.Results.VelocityModel;

% Station locations (column vector)
S_NORTHING = p.Results.Northing(:)';
S_EASTING  = p.Results.Easting(:)';
S_DEPTH    = zeros(size(S_EASTING));

% Event parameters.
STRIKE = STRIKE(:);
DIP = DIP(:);
RAKE = RAKE(:);
ANGLE = p.Results.ANGLE(:);
if isnan(ANGLE)
  ANGLE = zeros(size(STRIKE));
elseif isscalar(ANGLE)
  ANGLE = ANGLE * ones(size(STRIKE));
end
sigma = p.Results.PoissonRatio;
den = p.Results.Density;
M0 = 10.^(3*(p.Results.MomentMagnitude+6.07)/2);
E_DEPTH = -p.Results.EventDepth * ones(size(STRIKE));
E_NORTHING = 0 * ones(size(STRIKE));
E_EASTING = 0 * ones(size(STRIKE));

% Other parameters
if strcmpi(p.Results.FileFormat,'raw')
  format1d = false;
else
  format1d = true;
end

picklost = p.Results.PicksLost;
BIAS = p.Results.StationBias;

fclose all;
delete('input.txt');
delete('output.txt');

for j=1:numel(STRIKE)
  
  % Trace seismic ray parameters and calculate takeoff angles, azimuth,
  % distance and angle of incidence using 1D velocity model.
  
  e_depth = E_DEPTH(j);
  e_northing = E_NORTHING(j);
  e_easting = E_EASTING(j);
  
  H = S_DEPTH - e_depth;   % depth differences [m]
  RXY = sqrt( (S_NORTHING - e_northing).^2 + (S_EASTING - e_easting).^2); % epicentral distance [m]
  AZM = mod(360 + 180 * atan2(S_EASTING - e_easting, S_NORTHING - e_northing) / pi, 360)';
  
  %---- Perform ray tracing using 1D velocity model.
  
  % Generate input file for TT1D ray-tracing program.
  fid = fopen('input.txt','w');
  fprintf(fid,'%d\r\n',size(VMODEL,1));
  fprintf(fid,'%1.2f\r\n',VMODEL(:,1));
  fprintf(fid,'%1.2f\r\n',VMODEL(:,2));
  fprintf(fid,'DATA\r\n');
  fprintf(fid,'%f %f %f\r\n',[zeros(size(H)); H; RXY]/1000);
  fclose(fid);
  
  % Call ray-tracing program (fociMT undocumented feature).
  system('focimt.exe -m input.txt -o output.txt');
  
  % Get the output!
  DATA = load('output.txt');
  DATA = DATA(1:end-1,:); %%%% Temporary!!!
  TKO = DATA(:,6);
  AOI = DATA(:,7);
  R   = DATA(:,9) * 1000;
  
  % Estimate velocity in the source.
  velocity = 0;
  for i=size(VMODEL,1):-1:1
    if e_depth <= -VMODEL(i,1) * 1000
      velocity = VMODEL(i,2) * 1000;
      break;
    end
  end
  
  % Produce radiation pattern (pure shear faulting)
  RP = rpgen(STRIKE(j), DIP(j), RAKE(j), ANGLE(j), sigma, TKO, AZM);
  
  % Modify amplitude as it would be recorded at vertical sensor.
  RP_ZZ = RP .* cos(AOI * pi / 180);
  
  % Rescale amplitude to seismic moment.
  AMPLITUDE = RP_ZZ ./ 4 / pi / den / velocity ^ 3 ./ R * M0;
  
  % Plotting over focal sphere.
  %     [X,Y,Z] = sph2cart(pi / 2 - AZM * pi/180, pi / 2 - AOI * pi/180, 1);
  %     hold on;
  %     i = find(RP_ZZ >=0 ); plot3(X(i),Y(i),Z(i),'or');
  %     i = find(RP_ZZ < 0 ); plot3(X(i),Y(i),Z(i),'ob');
  %     hold off;
  %     for i=1:length(X(:))
  %       text(X(i),Y(i),Z(i),num2str(i),'FontSize',8,'HorizontalAlignment','center','Verticalalignment','bottom');
  %       text(X(i),Y(i),Z(i),['I' num2str(AOI(i)) ' T' num2str(TKO(i)) ' A' num2str(AZM(i))],'FontSize',8,'HorizontalAlignment','center','Verticalalignment','top');
  %     end
  %     axis vis3d;
  %     xlabel('x [m]');
  %     ylabel('y [m]');
  %     zlabel('z [m]');
  %     box on;
  %     grid on;
  %   pause;
  
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
  
  n_stations = length(TKO(:));
  n_lost = floor(picklost * n_stations);
  ACC = ones(size(TKO(:)));
  while n_lost > 0
    k = randi(n_stations);
    if ACC(k) == 1
      ACC(k) = 0;
      n_lost = n_lost - 1;
    end
  end
  
  if format1d
    % Generate file in 1D velocity model format.
    fid = fopen(filename,'a');
    fprintf(fid,'%s %d %1.0f %1.0f %1.0f %1.0f\r\n',fileid,sum(ACC), e_northing, e_easting, e_depth, den);
    for i = 1:n_stations
      if(ACC(i))
        for k=1:size(BIAS,1)
          if i == BIAS(k,1)
            AMPLITUDE(i) = AMPLITUDE(i) * BIAS(k,2);
          end
        end
        fprintf(fid,'S%02d Z P %22.14g  %14.0f %14.0f %14.0f\r\n', ...
          i, AMPLITUDE(i), S_NORTHING(i), S_EASTING(i), S_DEPTH(i));
      end
    end
    fclose(fid);
  else
    % Generate file in raw format
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
          i, AMPLITUDE(i), AZM(i), AOI(i), TKO(i), velocity, R(i), den);
      end
    end
    fclose(fid);
  end
end

fclose all;
delete('input.txt');
delete('output.txt');


