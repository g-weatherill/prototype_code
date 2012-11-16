function [vAlias] = import_alias(sFilename)
% function [vAlias] = import_alias(sFilename)
% -------------------------------------------
% Import station alias list of the SCSN/CISN
% Aliaslist from http://www.data.scec.org/stations/stamapping.html
%
% Input parameters:
%   sFilename       Filename of the aliaslist
%
% Output parameters:
%   vAlias          Vector of aliases
%
% Copyright (C) 2005 by Jochen Woessner & Danijel Schorlemmer
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
vAlias = [];

% Read in alias list
hFile = fopen(sFilename, 'r');
while ~feof(hFile)
  sLine = fgetl(hFile);
  [rAlias] = read_alias(sLine);
  nLastEntry = length(vAlias);
  vAlias = [vAlias, rAlias];
end;
fclose(hFile);

% ------------------------------------
% Helper functions
function [rAlias] = read_alias(sLine);

rAlias.sNetwork = strtrim(sLine(1:2));
rAlias.sAlias = strtrim(sLine(8:11));
rAlias.sName = strtrim(sLine(16:20));

