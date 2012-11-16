function [vUsedList, vDistances] = pmc_SortStationsForDistances(vList, fLon, fLat, fDepth);
% function [vUsedList, vDistances] = pmc_SortStationsForDistances(vList, fLon, fLat, fDepth)
% ------------------------------------------------------------------------------------------
% Computes distances of all stations in list to a given point and sorts
% station accordingly
%
% Input parameters:
%   vList           List of stations (from import_sc)
%   fLon            Longitude of point
%   fLat            Latitude of point
%   fDepth          Depth of point
%
% Output parameters:
%   vUsedList       Sorted list of stations
%   vDistances      Vector of distances of stations to the given point
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
vDistances = [];

% Compute distances
for nSta = 1:length(vList)
  fDistXY = sqrt(((fLon-vList(nSta).fLon)*cos(pi/180*vList(nSta).fLat)*111).^2 + ((fLat-vList(nSta).fLat)*111).^2);
  vDistances(nSta) = sqrt((vList(nSta).fElevation/1000 + fDepth)^2 + fDistXY^2);
end;

% Sort for distance
[vDummy, vIndices] = sort(vDistances);
vUsedList = vList(vIndices);
vDistances = vDistances(vIndices);
