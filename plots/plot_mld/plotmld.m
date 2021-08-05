%%  plotMld - plotting kriged mld model inputs
%%-------------------------------------------------------------------------
%%  configure 
cd('/Users/perrindavidson/Research/whoi/current/globalThPOCModels/code/globalThPOCModels/plots/');
close('all'); 
clear

%%  add paths
addpath /Users/perrindavidson/Documents/MATLAB/toolboxes/plotting/contourOcean/

%%  set outputs
outputpath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/outputs/';
plotpath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/plots/';

%%  load data
%   set filename ::
filename = [outputpath 'makemld/mldok.asc']; 
filenamevar = [outputpath 'makemld/mldokvar.asc']; 

% Read your data
[A, R]= readgeoraster(filename, 'coordinateSystemType', 'geographic');
[Avar, ~]= readgeoraster(filenamevar, 'coordinateSystemType', 'geographic');

%%  make grids
%   get limits ::
lonlims = R.LongitudeLimits; 
latlims = R.LatitudeLimits; 

%   get spacing ::
dlon = R.SampleSpacingInLongitude;
dlat = R.SampleSpacingInLatitude; 

%   make vectors ::
lon = lonlims(1) : dlon : lonlims(2); 
lat = latlims(1) : dlat : latlims(2);

%   make grids ::
[lon, lat] = meshgrid(lon, lat);

%   get data ::
mld = rot90(A, 2); 
mldvar = rot90(Avar, 2); 

%%  plot mld
%   set limits
%%% latitude ::
latlims = [-90, 90];

%%% plotting ::
plotlims = 0 : 10 : 300;

%   set colormap ::
colormaptype = 'viridis';

%   set colorbar position ::
colorbarpos = 'east';

%   set title ::
titlename = '';

%   set colorbar title ::
colorbartitle = 'Kriged MLD (m)';

%   contour plot ::
contourOcean(lon, lat, mld, ...
             latlims, plotlims, ...
             colormaptype, colorbarpos, ...
             titlename, colorbartitle);
         
%%  plot mld variance
%   set limits
%%% latitude ::
latlims = [-90, 90];

%%% plotting ::
plotlims = 400 : 10 : 700;

%   set colormap ::
colormaptype = 'viridis';

%   set colorbar position ::
colorbarpos = 'east';

%   set title ::
titlename = '';

%   set colorbar title ::
colorbartitle = 'Kriged MLD Variance (m)';

%   contour plot ::
contourOcean(lon, lat, mldvar, ...
             latlims, plotlims, ...
             colormaptype, colorbarpos, ...
             titlename, colorbartitle);
         
%%  end plotting