function [vList] = import_assignpicks(vStationList, vEventList);
% function [vList] = import_assignpicks(vStationList, vEventList)
% ---------------------------------------------------------------
% Assigns picks from the eventlist to the stations. 
%
% Input parameters:
%   vStationList    Stationlist (imported with import_stationlist)
%   vEventList      Eventlist (imported with import_scsnphase)
%
% Output parameters:
%   vList           Vector of stations with pickinformation. Same structure 
%                   as vStationlist with additional field
%                   vList.mPick(1)   Distance of earthquake to station
%                   vList.mPick(2)   Magnitude of event
%                   vList.mPick(3)   Event picked at station (boolean)
%                   vList.mPick(4)   Event index (for lookup in vEventList)
%                   vList.mPick(5)   Distance converted to magnitude units
%                   vList.nPick      Number of pick-information of station
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

% Init
vList = [];

% Loop thru stations
nNumEvents = length(vEventList);
for nStation = 1:length(vStationList)
  rStation = vStationList(nStation);
  mPick = zeros(nNumEvents, 5);
  nPick = 0;
  disp(['Station: ' num2str(nStation)]);
  for nEvent = 1:length(vEventList);
    if rem(nEvent,10000) == 0
      disp(['Event: ' num2str(nEvent)]);
    end;
    rEvent = vEventList(nEvent);
    if event_during_ontime(rEvent, rStation)
      fDistanceXY = sqrt(((rEvent.lon-rStation.fLon)*cos(pi/180*rStation.fLat)*111).^2 + ((rEvent.lat-rStation.fLat)*111).^2);
      fDistanceXYZ = sqrt(fDistanceXY^2 + (-(rStation.fElevation/1000) - rEvent.depth)^2);
      fMagnitude = rEvent.mag;
      if any_picks(rEvent)
        mPick(nPick+1,:) = [fDistanceXYZ fMagnitude event_picked(rEvent, rStation) nEvent calc_magnitudefromdistance(fDistanceXYZ)];
        nPick = nPick + 1;
      end;
    end;
  end;
  rStation.mPick = mPick(1:nPick,:);
  rStation.nPick = nPick;
  vList = [vList; rStation];
end;

% ---

function [bInOnTime] = event_during_ontime(rEvent, rStation)

bInOnTime = 0;
for nCnt = 1:length(rStation.mOnTime(:,1))
  vOnTime = rStation.mOnTime(nCnt,:);
  if (rEvent.datenum >= vOnTime(1)) & (rEvent.datenum <= vOnTime(2))
    bInOnTime = 1;
  end;
end;  

% ---

function [bPicked] = event_picked(rEvent, rStation)

bPicked = 0;

sName = rStation.sName;
for nCnt = 1:length(rEvent.mPicks)
  if strcmp(strtrim(sName), strtrim(rEvent.mPicks(nCnt).sStation))
    bPicked = 1;
  end;
end;

% ---

function [bAnyPicks] = any_picks(rEvent)

bAnyPicks = ~isempty(rEvent.mPicks);