function writeinput(filename, Input)
%WRITEINPUT Converts ASCII input file in RAW or 1D velocity format to input cell array
%   Internal function for create input file read by fociMT application.
%
%   part of hybridMT package 
%   <a href="matlab:open('html/doc_writeinput.html')">Reference page for writeinput</a>

%   Copyright 2015 Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%                  Patricia Martinez-Garzon <patricia@gfz-potsdam.de>
%
%   $Revision: 1.0.2 $  $Date: 2015.09.11 $

if exist(filename,'file');
  delete(filename);
end

[fid,errmsg] = fopen(filename,'w');
if fid == -1
  fclose(fid);
  error(errmsg);
end

for j=1:length(Input)
  if strcmpi(Input{j}.format,'raw')
    if Input{j}.matrixmode == false
      fprintf(fid,'%s %d\r\n',Input{j}.event_id,Input{j}.n_phases);
      for i=1:Input{j}.n_phases
        fprintf(fid,'%6s %2s %2s %22.14g  %6.1f %6.1f %6.1f %5.0f %7.2f %5.0f\r\n', Input{j}.Phase(i).station, ...
          Input{j}.Phase(i).component,Input{j}.Phase(i).phase, Input{j}.Phase(i).omega, ...
          Input{j}.Phase(i).azimuth, Input{j}.Phase(i).aoi, Input{j}.Phase(i).takeoff, Input{j}.Phase(i).v,...
          Input{j}.Phase(i).r, Input{j}.Phase(i).density);
      end
    else
      fprintf(fid,'%s %d\r\n',Input{j}.event_id,Input{j}.n_phases);
      for i=1:Input{j}.n_phases
        fprintf(fid,'%6s %2s %2s %22.14g  %6.1f %6.1f %6.1f %5.0f %7.2f %5.0f\r\n', Input{j}.Station{i}, ...
          Input{j}.Component{i},Input{j}.Phase{i}, Input{j}.OMEGA(i), ...
          Input{j}.AZIMUTH(i), Input{j}.AOI(i), Input{j}.TAKEOFF(i), Input{j}.V(i),...
          Input{j}.R(i), Input{j}.DENSITY(i));
      end
    end
  elseif strcmpi(Input{j}.format,'vel1d')
    if Input{j}.matrixmode == false
      fprintf(fid,'%s %d %1.0f %1.0f %1.0f %1.0f\r\n',Input{j}.event_id,Input{j}.n_phases, ...
        Input{j}.e_northing,Input{j}.e_easting, Input{j}.e_z, Input{j}.density);
      for i=1:Input{j}.n_phases
        fprintf(fid,'%6s %2s %2s %22.14g  %14.0f %14.0f %14.0f\r\n', Input{j}.Phase(i).station, ...
          Input{j}.Phase(i).component,Input{j}.Phase(i).phase, Input{j}.Phase(i).omega, ...
          Input{j}.Phase(i).s_northing, Input{j}.Phase(i).s_easting, Input{j}.Phase(i).s_z);
      end
    else
      fprintf(fid,'%s %d %1.0f %1.0f %1.0f %1.0f\r\n',Input{j}.event_id,Input{j}.n_phases, ...
        Input{j}.e_northing,Input{j}.e_easting, Input{j}.e_z, Input{j}.density);
      for i=1:Input{j}.n_phases
        fprintf(fid,'%6s %2s %2s %22.14g  %14.0f %14.0f %14.0f\r\n', Input{j}.Station{i}, ...
          Input{j}.Component{i},Input{j}.Phase{i}, Input{j}.OMEGA(i), ...
          Input{j}.S_NORTHING(i), Input{j}.S_EASTING(i), Input{j}.S_Z(i));
      end
    end
  else
    error('Input file format not recognized.');
  end
end
fclose(fid);

