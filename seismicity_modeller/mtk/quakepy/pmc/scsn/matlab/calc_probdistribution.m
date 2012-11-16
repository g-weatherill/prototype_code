function [rStation] = calc_probdistribution(rStation, bShow, nStyle, bSmooth);
% function [rStation] = calc_probdistribution(rStation, bShow, nStyle, bSmooth)
% -----------------------------------------------------------------------------
% Computes probability distribution for a given station
%
% Input parameters:
%   rStation          Station for which the probability distribution
%                     should be computed
%   bShow             Display the probability distribution
%   nStyle            Binning of probability distribution
%                     0: M=0:0.2:4; D=5:5:200; #=840;
%                        Low resolution for quick inspection
%                     1: M=0:0.05:4; D=1:1:200; #=16200;
%                        Classical resolution
%                     2: M=0:0.1:5; D=2:2:600; #=15300;
%                        Medium resolution with large range (Italy)
%   bSmooth           Smooth probability matrix according to
%                     Bachmann, Schorlemmer, and Kissling [2007]
%
% Output parameters:
%   rStation          Stationdata with rStation.mProbability
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

% Set defaults
if ~exist('bShow', 'var')
  bShow = 0;
end;
if ~exist('nStyle', 'var')
  nStyle = 1;
end;
if ~exist('bSmooth', 'var')
  bSmooth = 0;
end;

% Set size
if nStyle == 0
  vValues = ones(840,2)*nan;
  vMagnitudes = 0:0.2:4;
  vDistances = 5:5:200;
elseif nStyle == 1
  vValues = ones(16200,2)*nan;
  vMagnitudes = 0:0.05:4;
  vDistances = 1:1:200;
else
  vValues = ones(15300,2)*nan;
  vMagnitudes = 0:0.1:5;
  vDistances = 2:2:600;
end;

% Compute probabilities
nCnt = 1;
for fMagnitude = vMagnitudes
  for fDistance = vDistances
    [fProbability, nNumSample] = calc_probability(rStation, fDistance, fMagnitude);
    vValues(nCnt,:) = [fProbability nNumSample];
    nCnt = nCnt + 1;
  end;
end;

% Store parameters and probabilities
rStation.vBinX = vMagnitudes;
rStation.vBinY = vDistances;
rStation.mProbability = reshape(vValues(:,1), length(rStation.vBinY), length(rStation.vBinX));
rStation.nStyle = nStyle;
rStation.mNumSample = reshape(vValues(:,2), length(rStation.vBinY), length(rStation.vBinX));

% Smooth probability distribution
if bSmooth
  rStation.mProbability = calc_smooth_probdistribution(rStation.mProbability);
end;

% Show probability plot
if bShow
  show_probdistribution(rStation);
end;
