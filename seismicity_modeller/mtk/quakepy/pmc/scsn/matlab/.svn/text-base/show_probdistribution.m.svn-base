function [hFigure] = show_probdistribution(rStation, sFilename, bHidden);
% function [hFigure] = show_probdistribution(rStation, sFilename, bHidden)
% ------------------------------------------------------------------------
% Displays probability distribution for a given station
%
% Input parameters:
%   rStation          Station
%   sFilename         If given, image of probability distribution will
%                     be saved into a file with this name
%   bHidden           Hide display of image
%                     0: Display image
%                     1: Hide display (for automated image generation)
%
% Output parameters:
%   hFigure           Handle to figure
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

if ~exist('bHidden', 'var')
  bHidden = 0;
end;

if bHidden
  hFigure = figure('visible', 'off');
else
  hFigure = figure;
end;

surface(rStation.vBinX, rStation.vBinY, zeros(size(rStation.mProbability)), rStation.mProbability);
box on;
shading flat;

set(gca, 'Layer', 'top');
set(gca, 'xlim', [0.2 4], 'xtick', [0.2 0.5 1 1.5 2 2.5 3 3.5 4]);
if (rStation.nStyle == 0) | (rStation.nStyle == 1)
  set(gca, 'xlim', [0.2 4], 'xtick', [0.2 0.5 1 1.5 2 2.5 3 3.5 4]);
  set(gca, 'ylim', [0 200], 'ytick', [0 50 100 150 200]);
else % rStation == 2
  set(gca, 'xlim', [0 5], 'xtick', [0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5]);
  set(gca, 'ylim', [0 600], 'ytick', [0 100 200 300 400 500 600]);
end;
set(gca, 'linewidth', 2);
set(gca, 'tickdir', 'out', 'ticklength', [0.02 0.05]);
set(gca, 'fontname', 'Helvetica', 'fontsize', 16);
set(gca, 'clim', [0 1]);

% Y-label
hYLabel = text('String', 'Distance [km]');
set(hYLabel, 'FontName', 'Helvetica');
set(hYLabel, 'FontSize', 16);
set(hYLabel, 'FontWeight', 'normal');
set(hYLabel, 'FontAngle', 'normal');
set(gca, 'YLabel', hYLabel);

% X-label
hXLabel = text('String', 'Magnitude');
set(hXLabel, 'FontName', 'Helvetica');
set(hXLabel, 'FontSize', 16);
set(hXLabel, 'FontWeight', 'normal');
set(hXLabel, 'FontAngle', 'normal');
set(gca, 'XLabel', hXLabel);

% Colorbar
colormap(abs(flipud(gui_Colormap_ReadPovRay('rastafari_edited.pov', 128))));
hBar = colorbar;
set(hBar, 'position', [0.17 0.485 0.04 0.4]);
set(hBar, 'linewidth', 2);
set(hBar, 'fontname', 'Helvetica', 'fontsize', 14);
set(hBar, 'ylim', [0 1]);
set(hBar, 'tickdir', 'out', 'ticklength', [0.03 0.06], 'ytick', [0 0.2 0.4 0.6 0.8 1]);

% Colorbar-label
if (rStation.nStyle == 0) | (rStation.nStyle == 1)
  hText = text(0.9, 188, 'Probability of recording');
else
  hText = text(1.1, 564, 'Probability of recording');
end;
set(hText, 'fontname', 'Helvetica', 'fontsize', 16);

% Export figure
if exist('sFilename', 'var')
  if strcmp('SVG', upper(sFilename(length(sFilename)-2:length(sFilename))))
    plot2svg_2d(sFilename, hFigure);
  else
    exportfig(hFigure, sFilename, 'color', 'cmyk');
  end;
end;

