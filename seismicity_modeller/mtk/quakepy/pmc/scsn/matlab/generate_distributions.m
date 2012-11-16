function [vNewList] = generate_distributions(vList, sDirectory, nStyle);
% function [vNewList] = generate_distributions(vList, sDirectory, nStyle)
% -----------------------------------------------------------------------
% Computes probability distributions for all stations in vList and
% generates figures of these distributions in a directory
%
% Input parameters:
%   vList             List of station for which the probability distribution
%                     should be computed
%   sDirectory        Directory for saving probability distribution figures
%   nStyle            Binning of probability distribution
%                     See calc_probdistribution.m for details
%
% Output parameters:
%   vNewList          List of stations with additional mProbability
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

% Initialize output list
vNewList = [];

for i = 1:length(vList)
  rStation = vList(i);
  disp(['Generate Matrix for station #' num2str(i)]);
  [rStation] = calc_probdistribution(rStation, 0, nStyle, 1);

  sNumP = num2str(sum(rStation.mPick(:,3)==1), '%5.0f');
  for n = 1:(5 - length(sNumP))
    sNumP = ['0' sNumP];
  end;
  sNumT = num2str(sum(rStation.mPick(:,3)==1) + sum(rStation.mPick(:,3)==0), '%5.0f');
  for n = 1:(5 - length(sNumT))
    sNumT = ['0' sNumT];
  end;
  sIdx = num2str(i);
  sName = rStation.sName;
  vNewList = [vNewList; rStation];

  sFilenameEPS = [sDirectory '/' sNumT '_' sNumP '_' sIdx '_' sName '.eps'];
  [hFigure] = show_probdistribution(rStation, sFilenameEPS, 1);
  close(hFigure);

  %unix(['unset LD_LIBRARY_PATH;/usr/bin/convert ' sFilenameEPS ' ' sFilenamePNG]);
  %unix(['unset LD_LIBRARY_PATH;/bin/rm ' sFilenameEPS]);
end;
