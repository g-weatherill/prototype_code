function run_sc(fDepth);
% Computes the probability of detection for events of given magnitude
% for the area of Italy.
%
% Input parameter:
%   fMagnitude       Magnitude of events for which detection probabilities
%                    are computed
%
% Copyright (C) 2007 by Danijel Schorlemmer
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

addpath('../../code');

% Load station detection probabilities and events
load('../../data/vListOnlyDistros.mat');
load('../../data/vEventList.mat');
% Set the date/time of probability computation
fDatenum = datenum(2007,07,01,00,00,00);
% Select stations in operation
[vUsedList] = pmc_SelectStationsInOperation(vList, fDatenum);
% Set computational parameters
fMinMag = 0;
fMaxMag = 4;
fMinLon = -122;
fMaxLon = -113.5;
fMinLat = 31.5;
fMaxLat = 38;
fDepth  = 7.5;
fStep   = 0.05;
fMagStep = 0.1;

mMPValues = ones(22401,1) .* NaN;
mOldValues = zeros(22401,5);
for fMagnitude = fMinMag:fMagStep:fMaxMag
  mValues = zeros(22401,5);
  nCnt = 0;
  % Iterate through the grid
  for fLon = fMinLon:fStep:fMaxLon;
    for fLat = fMinLat:fStep:fMaxLat;
      nCnt = nCnt + 1;
      if mOldValues(nCnt,4) > 0.99999
        fProbability = 1;
      else
        % Compute detection probabilities for each station
        [vTmpList, vDistances] = pmc_SortStationsForDistances(vUsedList, fLon, fLat, fDepth);
        [vProbabilities] = pmc_ComputeProbabilitiesPerStation(vTmpList, vDistances, fMagnitude, 1);
        % Compute the combined probability of detecting the event
        [fProbability, nNumUsedStations] = pmc_ComputeCombinedProbability(vProbabilities);
        if fProbability > 0.99999
          mMPValues(nCnt) = fMagnitude;
          disp(['MP-value: ' num2str(fMagnitude)]);
        end;
      end;
      % Store the value
      mValues(nCnt,:) = [fLon fLat fDepth fProbability fMagnitude];
      disp(['Node: ' num2str(fLon) '/' num2str(fLat) ' | ' num2str(fProbability)]);
    end;
  end;
  % Generate string from magnitude value for filenames
  sMagnitude = sprintf('%3.1f', fMagnitude);
  sDepth     = sprintf('%03g', fDepth);
  % Save results and pretty print
  save(['data/tmp.dat'], 'mValues', '-ascii');
  sFilename = ['data/sc.p' sMagnitude '-' sDepth '-ALL.dat'];
  unix(['more data/tmp.dat | ./pretty.print > ' sFilename]);
  unix('rm data/tmp.dat');
  % Export stations and events for GMT
  gmt_exportstations(vUsedList, 'data/stations.20070701.dat');
  gmt_exportevents(vEventList, ['data/events.m' sMagnitude '.dat'], fMagnitude);
  mOldValues = mValues;
end;
mValues(:,4) = mMPValues;
mValues(:,5) = [];
save(['data/sc.mp0.99999-' sDepth '.dat'], 'mValues', '-ascii');