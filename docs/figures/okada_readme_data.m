% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% Data of the figure of the README: the surface displacements that OkadaSol.m computes for a
% sample thrust fault, a rectangle 60 km long and 40 km wide, dipping at 13 degrees, whose
% deeper edge lies at a depth of 20 km along the x-axis, with a slip of 1 m and the rake
% 90 degrees (pure dip slip). All lengths are in metres and all angles in radians.
%
% Run it from this folder with GNU Octave or MATLAB:
%   octave --no-gui okada_readme_data.m        or, in MATLAB,   run okada_readme_data
% It writes x.csv and y.csv (the grid, in metres) and u.csv, v.csv, w.csv (the displacements
% along x, along y and upwards, in metres) beside itself; plot_readme_figure.py then draws
% okada-displacement.png from them.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'MatlabScript'));

global delta phi L W
delta = 13*pi/180;          % dip angle
phi   = 0;                  % strike angle: the fault runs along the x-axis
L     = 60e3;               % length of the fault along strike
W     = 40e3;               % width of the fault along dip

x0 = 0; y0 = 0;             % corner of the fault (on its deeper edge)
d  = 20e3;                  % depth of the fault at that corner
U  = 1.0;                   % slip
theta = pi/2;               % rake

x = linspace(-60e3, 120e3, 361);
y = linspace(-80e3, 120e3, 401);
[X, Y] = meshgrid(x, y);
[u, v, w] = OkadaSol(X, Y, x0, y0, d, U, theta);

dlmwrite(fullfile(here, 'x.csv'), x(:), 'precision', '%.16e');
dlmwrite(fullfile(here, 'y.csv'), y(:), 'precision', '%.16e');
dlmwrite(fullfile(here, 'u.csv'), u, 'precision', '%.16e');
dlmwrite(fullfile(here, 'v.csv'), v, 'precision', '%.16e');
dlmwrite(fullfile(here, 'w.csv'), w, 'precision', '%.16e');
[top, i] = max(w(:)); [bottom, j] = min(w(:));
fprintf('uplift %.4f m at (%.1f, %.1f) km, subsidence %.4f m at (%.1f, %.1f) km\n', ...
        top, X(i)/1e3, Y(i)/1e3, bottom, X(j)/1e3, Y(j)/1e3);
