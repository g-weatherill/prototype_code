function [vProbabilities] = pmc_ComputeProbabilitiesPerStation(vList, vDistances, fMagnitude, bLookup)
% function [vProbabilities] = pmc_ComputeProbabilitiesPerStation(vList, vDistances, fMagnitude, bLookup)
% ---------------------------------------------------------------------------------------------
% Computes probabilities of detection for an event of given magnitude
% and given distance for each station
%
% Input parameters:
%   vList           List of stations (from import_sc)
%   vDistances      Vector of distances of stations to a point of interest
%   fMagnitude      Magnitude of an event for which the probabilities are computed
%   bLookup         0: Sample probabilities; 1: Lookup probabilities (vList.mPlusProbability)
%
% Output parameters:
%   vProbabilities  Vector of probabilities (same order as vList)
%                   vProbabilities(:,1) : Probability of detection
%                   vProbabilities(:,2) : Probability of non-detection
%                   vProbabilities(:,3) : Distance from point to station
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
nNumStations = length(vList);
vProbabilities = ones(nNumStations,3)*nan;

% Compute probabilites of detection of an event with fMagnitude at each station
for nSta = 1:nNumStations
  rStation = vList(nSta);
  if bLookup
    vProbabilities(nSta,1) = lookup_probability(rStation, vDistances(nSta), fMagnitude);
  else
    vProbabilities(nSta,1) = calc_probability(rStation, vDistances(nSta), fMagnitude);
  end;
end;

% Compute probabilities of non-detection
vProbabilities(:,2) = 1 - vProbabilities(:,1);

% Add distances to the result vector
vProbabilities(:,3) = vDistances(1:nNumStations)';
