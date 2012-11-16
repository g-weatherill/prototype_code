function [vStationList] = import_checkalias(vStationList, rAlias);
% function [vStationList] = import_checkalias(vStationList, rAlias)
% -----------------------------------------------------------------
% Remove alias stations from list, add on-/off-times to station name and
% add pick information
%
% Input parameters:
%   vStationList    Record of stations with picks
%   rAlias          Record Alias stations (see import_alias for details)
%
% Output parameters:
%   vStationList    Vector of stations with alias stations removed
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

% Check for Alias stations
vDel = [];
for  nCnt = 1:length(rAlias)
    for nC = 1:length(vStationList)
        bMatch = strcmp(vStationList(nC).sName,rAlias(nCnt).sAlias);
        if bMatch
            for nX = 1:length(vStationList)
                bNoma = strcmp(vStationList(nX).sName,rAlias(nCnt).sName);
                if bNoma
                    % Add on-time
                    vStationList(nX).mOnTime = [vStationList(nX).mOnTime; vStationList(nC).mOnTime];
                    % Add picks
                    vStationList(nX).mPick = [vStationList(nX).mPick; vStationList(nC).mPick];
                end;
            end;
            vDel = [vDel; nC];
        end;
    end;
end;
% Remove alias stations from list
vStationList(vDel) = [];
