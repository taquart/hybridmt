function Solution = readsolution(filename, solution, Solution, matrixmode)
%READSOLUTION Read seismic moment tensor solution file created by fociMT
%
%   part of hybridMT package
%   <a href="matlab:open('html/doc_readsolution.html')">Reference page for readsolution</a>

%   Copyright 2015 Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%                  Patricia Martinez-Garzon <patricia@gfz-potsdam.de>
%
%   $Revision: 1.0.3 $  $Date: 2015.10.02 $

filename1 = [filename '.asc'];
filename2 = [filename '-u.asc'];

[fid,errmsg] = fopen(filename1,'r');

if fid == -1
  fclose(fid);
  error(errmsg);
end

%---- Read solution file.

j = 1;
while 1
  evstr = ['Solution{j}.' solution '(k).'];
  k = 1;
  
  event_id = fscanf(fid,'%s',1);
  n = fscanf(fid,'%d',1);
  if(isempty(event_id) || isempty(n))
    break;
  end
  
  Solution{j}.event_id = event_id;
  Solution{j}.n_trials = n;
  Solution{j}.calculation_dt = now;
  
  if matrixmode
    Data = textscan(fid,['%s %d  %f %f %f  %f %f %f %f  %f %f %f %f %f %f ' ...
      '%f %f %f %f %f %f  %s  %f %f %f %f %f %f  %f %f %f %f %f %f  %f'],n); %#ok<NASGU>
%     switch solution
%       case 'full'
%         Solution{j}.full.Type=Data{1};
%         Solution{j}.full.STATION_ID=Data{2};
%         Solution{j}.full.ISO=Data{3};
%         Solution{j}.full.CLVD=Data{4};
%         Solution{j}.full.DC=Data{5};
%         Solution{j}.full.M0=Data{6};
%         Solution{j}.full.MT=Data{7};
%         Solution{j}.full.M0ERRMAX=Data{8};
%         Solution{j}.full.MW=Data{9};
%         Solution{j}.full.P=[Data{10:11}];
%         Solution{j}.full.T=[Data{12:13}];
%         Solution{j}.full.B=[Data{14:15}];
%         Solution{j}.full.F1=[Data{16:18}];
%         Solution{j}.full.F2=[Data{19:21}];
%         Solution{j}.full.Fault=Data{22};
%         Solution{j}.full.MXX=[Data{23:28}];
%         Solution{j}.full.CXX=[Data{29:34}];
%         Solution{j}.full.RMSERROR=[Data{35}];
%       case 'deviatoric'
%         Solution{j}.deviatoric.Type=Data{1};
%         Solution{j}.deviatoric.STATION_ID=Data{2};
%         Solution{j}.deviatoric.ISO=Data{3};
%         Solution{j}.deviatoric.CLVD=Data{4};
%         Solution{j}.deviatoric.DC=Data{5};
%         Solution{j}.deviatoric.M0=Data{6};
%         Solution{j}.deviatoric.MT=Data{7};
%         Solution{j}.deviatoric.M0ERRMAX=Data{8};
%         Solution{j}.deviatoric.MW=Data{9};
%         Solution{j}.deviatoric.P=[Data{10:11}];
%         Solution{j}.deviatoric.T=[Data{12:13}];
%         Solution{j}.deviatoric.B=[Data{14:15}];
%         Solution{j}.deviatoric.F1=[Data{16:18}];
%         Solution{j}.deviatoric.F2=[Data{19:21}];
%         Solution{j}.deviatoric.Fault=Data{22};
%         Solution{j}.deviatoric.MXX=[Data{23:28}];
%         Solution{j}.deviatoric.CXX=[Data{29:34}];
%         Solution{j}.deviatoric.RMSERROR=[Data{35}];
%       case 'dc'
%         Solution{j}.dc.Type=Data{1};
%         Solution{j}.dc.STATION_ID=Data{2};
%         Solution{j}.dc.ISO=Data{3};
%         Solution{j}.dc.CLVD=Data{4};
%         Solution{j}.dc.DC=Data{5};
%         Solution{j}.dc.M0=Data{6};
%         Solution{j}.dc.MT=Data{7};
%         Solution{j}.dc.M0ERRMAX=Data{8};
%         Solution{j}.dc.MW=Data{9};
%         Solution{j}.dc.P=[Data{10:11}];
%         Solution{j}.dc.T=[Data{12:13}];
%         Solution{j}.dc.B=[Data{14:15}];
%         Solution{j}.dc.F1=[Data{16:18}];
%         Solution{j}.dc.F2=[Data{19:21}];
%         Solution{j}.dc.Fault=Data{22};
%         Solution{j}.dc.MXX=[Data{23:28}];
%         Solution{j}.dc.CXX=[Data{29:34}];
%         Solution{j}.dc.RMSERROR=[Data{35}];
%     end
        eval([ evstr 'Type=Data{1};']);
        eval([ evstr 'STATION_ID=Data{2};']);
        eval([ evstr 'ISO=Data{3};']);
        eval([ evstr 'CLVD=Data{4};']);
        eval([ evstr 'DC=Data{5};']);
        eval([ evstr 'M0=Data{6};']);
        eval([ evstr 'MT=Data{7};']);
        eval([ evstr 'M0ERRMAX=Data{8};']);
        eval([ evstr 'MW=Data{9};']);
        eval([ evstr 'P=[Data{10:11}];']);
        eval([ evstr 'T=[Data{12:13}];']);
        eval([ evstr 'B=[Data{14:15}];']);
        eval([ evstr 'F1=[Data{16:18}];']);
        eval([ evstr 'F2=[Data{19:21}];']);
        eval([ evstr 'Fault=Data{22};']);
        eval([ evstr 'MXX=[Data{23:28}];']);
        eval([ evstr 'CXX=[Data{29:34}];']);
        eval([ evstr 'RMSERROR=[Data{35}];']);
  else
    for i=1:n
      %                            -D------  -W---------  -A---------------
      %  -F---------------  -T
      Data = textscan(fid,['%s %d  %f %f %f  %f %f %f %f  %f %f %f %f %f %f ' ...
        '%f %f %f %f %f %f  %s  %f %f %f %f %f %f  %f %f %f %f %f %f'],1); %#ok<NASGU>
      eval([ evstr 'Type=Data{1}{:};']);
      eval([ evstr 'station_id=Data{2};']);
      eval([ evstr 'iso=Data{3};']);
      eval([ evstr 'clvd=Data{4};']);
      eval([ evstr 'dc=Data{5};']);
      eval([ evstr 'm0=Data{6};']);
      eval([ evstr 'mT=Data{7};']);
      eval([ evstr 'm0errmax=Data{8};']);
      eval([ evstr 'mw=Data{9};']);
      eval([ evstr 'P=[Data{10:11}];']);
      eval([ evstr 'T=[Data{12:13}];']);
      eval([ evstr 'B=[Data{14:15}];']);
      eval([ evstr 'F1=[Data{16:18}];']);
      eval([ evstr 'F2=[Data{19:21}];']);
      eval([ evstr 'Fault=Data{22}{:};']);
      eval([ evstr 'MXX=[Data{23:28}];']);
      eval([ evstr 'CXX=[Data{29:34}];']);
      eval([ evstr 'RMSERROR=[Data{35}];']);
      k = k + 1;
    end
  end
  j = j + 1;
end
fclose(fid);

%---- Read displacement information.
[fid,errmsg] = fopen(filename2,'r');

if fid == -1
  fclose(fid);
  error(errmsg);
end

j = 1;
while 1
  event_id = fscanf(fid,'%s',1); % event id
  n = fscanf(fid,'%d',1); % number of realizations
  if(isempty(event_id) || isempty(n))
    break;
  end
  for i=1:n
    sol_type = fscanf(fid,'%s',1); %#ok<NASGU>
    rejected_station_no = fscanf(fid,'%d',1); %#ok<NASGU>
    n_u = fscanf(fid,'%d',1);
    Data = textscan(fid,'%s %f %f', n_u);
    
    if matrixmode
      evstr = ['Solution{j}.' solution '.'];
      if i == 1
        % first is always "N" solution, prepare empty matrices.
        eval([ evstr 'Station=Data{1}'';']);
        eval([ evstr 'UMEASURED=nan(n,length(Data{1}));']);
        eval([ evstr 'UTH=nan(n,length(Data{1}));']);
        StationList = Data{1};
      end
      for m=1:length(Data{1})
        k = find(strcmp(Data{1}{m},StationList));
        if ~isempty(k)
          eval([ evstr 'UMEASURED(i,k)=Data{2}(m);']);
          eval([ evstr sprintf('UTH(%d,%d)=Data{3}(%d);',i,k,m)]);
        else
          error('!!!');
        end
      end
    else
      evstr = ['Solution{j}.' solution '(i).'];
      
      eval([ evstr 'Station=Data{1}'';']);
      eval([ evstr 'Umeasured=Data{2}'';']);
      eval([ evstr 'Uth=Data{3}'';']);
    end
  end
  j = j + 1;
end
fclose(fid);


