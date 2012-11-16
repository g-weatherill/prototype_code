function [vStationList] = import_stationlist(sFilename, fStartdate, fEnddate);
% function [vStationList] = import_stationlist(sFilename, fStartdate, fEnddate)
% -----------------------------------------------------------------------------
% Import list of stations from the SCSN/CISN network for the time period
% between fStartdate and fEnddate
%
% Input parameters:
%   sFilename       Filename of the stationlist
%   fStartdate      Startdate in Matlab form, e.g. datenum('19-May-2000')
%   fEnddate        Enddate in Matlab form, e.g. datenum('19-May-2000')
%
% Output parameters:
%   vStationList    Vector of active stations
%
% Copyright (C) 2005 by Danijel Schorlemmer & Jochen Woessner
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the
% Free Software Foundation, Inc.,
% 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

% Container
vStationList = [];

% Read in station list
hFile = fopen(sFilename, 'r');
while ~feof(hFile)
    sLine = fgetl(hFile);
    [rStation] = read_station(sLine);
    nLastEntry = length(vStationList);
    if station_be_used(rStation)
        if isempty(vStationList)
            vStationList = [vStationList; rStation];
        elseif strcmp(vStationList(nLastEntry).sName, rStation.sName)
            vStationList(nLastEntry).mOnTime = [vStationList(nLastEntry).mOnTime; rStation.mOnTime];
        else
            vStationList = [vStationList, rStation];
        end;
    end;
end;
fclose(hFile);

% Remove stations which have not been operating between fStartdate and fEnddate
for nStation = length(vStationList):-1:1
  rStation = vStationList(nStation);
  if ~operating_during_period(rStation, fStartdate, fEnddate);
    vStationList(nStation) = [];
  end;
end;

% --- Helper ----
function [rStation] = read_station(sLine);

rStation.sName = strtrim(sLine(5:9));
rStation.fLon = str2num(sLine(61:70));
rStation.fLat = str2num(sLine(51:59));
rStation.fElevation = str2num(sLine(72:76));
rStation.sNetwork = sLine(1:2);
% Seconds are set to zero as datenum format wants it!!
fOnDate = datenum(str2num(sLine(78:81)),str2num(sLine(83:84)),str2num(sLine(86:87)),0,0,0);
fOffDate = datenum(str2num(sLine(89:92)),str2num(sLine(94:95)),str2num(sLine(97:98)),0,0,0);
rStation.mOnTime = [fOnDate fOffDate];

% ---
function [bUsed] = station_be_used(rStation)

bUsed = 1;
% Stations from a mexican network
if strcmp(rStation.sNetwork, 'MX')
  bUsed = 0;
end;
% Temporary stations in California (for aftershock sequences) 
if strcmp(rStation.sNetwork, 'ZY')
  bUsed = 0;
end;
if strcmp(rStation.sNetwork, 'NR')
  bUsed = 0;
end;
if strcmp(rStation.sNetwork, 'RB')
  bUsed = 0;
end;
if strcmp(rStation.sNetwork, 'LI')
  bUsed = 0;
end;
% Stations from the Berkeley network
if strcmp(rStation.sNetwork, 'BK')
  bUsed = 0;
end;

% ---
function [bOperating] = operating_during_period(rStation, fStartdate, fEnddate);

bOperating = 0;
if ~isempty(rStation.mOnTime)
  for nCnt = 1:length(rStation.mOnTime(:,1))
    if (rStation.mOnTime(nCnt,1) <= fStartdate) & (rStation.mOnTime(nCnt,2) >= fEnddate)
      bOperating = 1;
    elseif (rStation.mOnTime(nCnt,1) >= fStartdate) & (rStation.mOnTime(nCnt,2) <= fEnddate)
      bOperating = 1;
    elseif (rStation.mOnTime(nCnt,1) <= fStartdate) & (rStation.mOnTime(nCnt,2) <= fEnddate) & (rStation.mOnTime(nCnt,2) >= fStartdate)
      bOperating = 1;
    elseif (rStation.mOnTime(nCnt,1) >= fStartdate) & (rStation.mOnTime(nCnt,2) >= fEnddate) & (rStation.mOnTime(nCnt,1) <= fEnddate)
      bOperating = 1;
    end;
  end;
end;
