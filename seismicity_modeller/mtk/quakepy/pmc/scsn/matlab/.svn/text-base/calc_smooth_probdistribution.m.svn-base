function [mSmoothProbability] = calc_smooth_probdistribution(mProbability)
% function [mSmoothProbability] = calc_smooth_probdistribution(mProbability)
% --------------------------------------------------------------------------
% Smoothes probability distribution according to Bachmann, Schorlemmer,
% and Kissling [2007]
%
% Input parameters:
%   mProbability         Probability disgtribution
%
% Output parameters:
%   rSmoothProbability   Smoothed probability distribution
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

% Initialize output
[nRows, nCols] = size(mProbability);
mSmoothProbability = ones(nRows, nCols) * nan;

% Smoth along distance
for nRow = 1:nRows
  mSmoothProbability(nRow,1) = mProbability(nRow,1);
  for nCol = 2:nCols
    if mProbability(nRow,nCol) < mSmoothProbability(nRow,nCol-1)
      mSmoothProbability(nRow,nCol) = mSmoothProbability(nRow,nCol-1);
    else
      mSmoothProbability(nRow,nCol) = mProbability(nRow,nCol);
    end;
  end;
end;

% Smooth along magnitudes
mProbability = mSmoothProbability;
for nCol=1:nCols
  mSmoothProbability(1,nCol) = mProbability(1,nCol);
  for nRow=(nRows-1):-1:1
    if mProbability(nRow,nCol) < mSmoothProbability(nRow+1,nCol)
      mSmoothProbability(nRow,nCol) = mSmoothProbability(nRow+1,nCol);
    else
      mSmoothProbability(nRow,nCol) = mProbability(nRow,nCol);
    end;
  end;
end;
