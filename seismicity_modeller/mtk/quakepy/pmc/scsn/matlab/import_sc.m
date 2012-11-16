function [vEventList, vStationList, vList] = import_sc;
% function [vEventList, vStationList, vList] = import_sc
% ------------------------------------------------------
% Wrapper-function to import stations, remove alias stations,
% import phase data, and assign picks to the stations
%
% Copyright (C) 2006 by Jochen Woessner & Danijel Schorlemmer
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

% Import events
[vEventList] = import_scsnphase(['../data/phase2001-200706.dat']);
save('../data/vEventList.mat','vEventList');
disp(['Number of events:' num2str(length(vEventList))]);

% Remove M=0 events
[vEventList] = import_removemagnitudezero(vEventList);
save('../data/vEventList.mat','vEventList');
disp(['Number of events:' num2str(length(vEventList))]);

% Import list of station aliases
[rAlias] = import_alias('../data/sc_alias.txt');
save('../data/tmp_Alias.mat','rAlias');
disp(['Number of aliases:' num2str(length(rAlias))]);

% Import all stations
[vStationList_all] = import_stationlist('../data/stationlist.dat',datenum('01-Jan-2001'),datenum('01-Jul-2007'));
save('../data/tmp_StationList_all.mat','vStationList_all');
disp(['Number of stations:' num2str(length(vStationList_all))]);

% Assign picks
[vList_all] = import_assignpicks(vStationList_all, vEventList);
save('../data/tmp_vList_all.mat','vList_all');
disp(['Number of stations:' num2str(length(vList_all))]);

% Remove aliases
[vList] = import_checkalias(vList_all, rAlias);
save('../data/vList.mat','vList');
disp(['Size of station list without alias stations:' num2str(length(vList))]);

% delete temporary files
!rm ../data/tmp_Alias.mat
!rm ../data/tmp_vList_all.mat
!rm ../data/tmp_StationList_all.mat