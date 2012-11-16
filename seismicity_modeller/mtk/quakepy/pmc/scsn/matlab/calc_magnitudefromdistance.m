function [fMagnitude] = calc_magnitudefromdistance(fDistance);
% function [fMagnitude] = calc_magnitudefromdistance(fDistance)
% -------------------------------------------------------------
% Converts distance to magnitude for pick sampling. Magnitude
% relation used at the SCSN implemented.
%
% Input parameters:
%   fDistance       Distance between event and station
%
% Output parameters:
%   fMagnitude      Corresponding magnitude
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

fMagnitude = (-log10(0.3173*exp(-0.00505*fDistance) * power(fDistance, -1.14)));