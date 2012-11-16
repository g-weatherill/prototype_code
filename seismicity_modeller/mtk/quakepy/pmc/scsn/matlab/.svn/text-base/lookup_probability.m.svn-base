function [fProbability, nNumSample] = lookup_probability(rStation, fDistance, fMagnitude);
% function [fProbability, nNumSample] = lookup_probability(rStation, fDistance, fMagnitude)
% ---------------------------------------------------------------------------------------
% Computes detection probability of a station given a distance/magnitude combination.
%
% Input parameters:
%   rStation        Station record (with pick-information) from vStationList
%   fDistance       Distance from hypothetical event to station
%   fMagnitude      Magnitude of hypothetical event
%
% Output parameters:
%   fProbability    Probability of detection
%   nNumSample      Number of pick-information sampled for calculation
%
% Copyright (C) 2005-2007 by Danijel Schorlemmer & Jochen Woessner
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

fMaxMagnitude = rStation.vBinX(length(rStation.vBinX));
if fMagnitude > fMaxMagnitude
  fMagnitude = fMaxMagnitude;
end;

fMaxDistance = rStation.vBinY(length(rStation.vBinY));
if (fMagnitude < 0) | (fDistance > fMaxDistance)
  fProbability = 0;
  nNumSample   = 0;
else
  if rStation.nStyle == 0
    nMagnitudeIndex = round(fMagnitude * 5) + 1;
    nDistanceIndex  = round(fDistance/5);
  elseif rStation.nStyle == 1
    nMagnitudeIndex = round(fMagnitude * 20) + 1;
    nDistanceIndex  = round(fDistance);
  else % rStation.nStyle == 2
    nMagnitudeIndex = round(fMagnitude * 10) + 1;
    nDistanceIndex  = round(fDistance/2);
  end;
  if nDistanceIndex == 0
    nDistanceIndex = 1;
  end;

  fProbability = rStation.mProbability(nDistanceIndex, nMagnitudeIndex);
  nNumSamples  = rStation.mNumSample(nDistanceIndex, nMagnitudeIndex);
end;

