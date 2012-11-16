function [fProbability, nNumUsedStations] = pmc_ComputeCombinedProbability(vProbabilities);
% function [fProbability, nNumUsedStations] = pmc_ComputeCombinedProbability(vProbabilities)
% ------------------------------------------------------------------------------------------
% Computes the combined probability that 4 or more stations detected a event.
% The probabilities of each single station are passed in the vProbabilities parameter.
% Use pmc_ComputeProbabilitiesPerStation for obtaining vProbabilities
%
% Input parameters:
%   vProbabilities  Vector of probabilities (from pmc_ComputeProbabilitiesPerStation)
%
% Output parameters:
%   fProbability    Combined probability that an event has been detected at
%                   4 or more stations, given the probabilities of each single station
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

% Remove stations with a probability of 0
vSel = vProbabilities(:,1) > 0;
vProbabilities = vProbabilities(vSel,:);
nNumUsedStations = length(vProbabilities(:,1));

if nNumUsedStations < 4
  % 4 stations needed for triggering
  fProbability = 0;
else
  V_STATION_INDICES = 1:nNumUsedStations;

  %Compute probability for picks at 0 stations
  vP = vProbabilities(:,2);
  fProd_0 = prod(vP);

  %Compute probabilities for picks at 1,2,3 stations
  for nNum = 1:3
    vProd(nNum) = 0;
    vCombinations = nchoosek(V_STATION_INDICES, nNum);
    for nCnt = 1:length(vCombinations(:,1))
      vP = vProbabilities(:,2);
      vComb = vCombinations(nCnt,:);
      vP(vComb) = vProbabilities(vComb,1);
      vProd(nNum) = vProd(nNum) + prod(vP);
    end;
  end;
  % Compute overall probability
  fProbability = 1 - fProd_0 - sum(vProd);
end;


