function [fProbability, nNumSample] = calc_probability(rStation, fDistance, fMagnitude);
% function [fProbability, nNumSample] = calc_probability(rStation, fDistance, fMagnitude)
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

mValues = zeros(length(rStation.mPick(:,1)), 4);
mValues(:,1) = rStation.mPick(:,3);

% Magnitude differences
mValues(:,2) = rStation.mPick(:,2) - fMagnitude;

% Translate difference in distances into a magnitude difference
fLogDistance = calc_magnitudefromdistance(fDistance);
mValues(:,3) = rStation.mPick(:,5) - fLogDistance;

% Calculate total difference in "magnitude units"
mValues(:,4) = sqrt(mValues(:,2).^2 + mValues(:,3).^2);

% Select only picks with a "magnitude difference" of less than 0.1
vSel = (mValues(:,4) < 0.1);
mSample = mValues(vSel,:);
nNumSample = length(mSample(:,1));
if nNumSample < 10
  mNonSample = mValues(~vSel,:);
  vSel = (mNonSample(:,2) <= 0) & (mNonSample(:,3) >= 0);
  mSelection = mNonSample(vSel,:);
  [vVals, vIndices] = sort(mSelection(:,4));
  mSortedSelection = mSelection(vIndices(:,1),:);
  if length(mSortedSelection(:,1)) < (10 - nNumSample)
    fProbability = 0;
  else
    mSample = [mSample; mSortedSelection(1:(10-nNumSample),:)];
    fProbability = sum(mSample(:,1)==1)/10;
    nNumSample = 10;
  end;
else
  fProbability = sum(mSample(:,1)==1)/nNumSample;
end;
