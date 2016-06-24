%%--------------------------Example 2 :SPLM-------------------------------%
% Author:      Benjamin Campforts Katholieke Universiteit Leuven 
%              <benjamin.campforts@kuleuven.be>
% 
%References
%
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo%
%References
%
% Campforts, B., and G. Govers (2015),Keeping the edge: A numerical method
% that avoids knickpoint smearing when solving the stream power law, J.
% Geophys. Res. Earth Surf., 120, doi:10.1002/2014JF003376.
%
% Campforts, B., Schwanghart W, and G. Govers (2015),TTLEM 1.0: A numerical
% package for accurate simulation of transient landscape evolution in
% MATLAB. GMD
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo%
%
%-------------------------------------------------------------------------%


% clearvars
clc
close all force

%Numerical method
numM=4;
%Parameters
K=5e-6; m=.42; n=1;kappa=0;
parameters=[K m n kappa];
%Spatial
dx=100; x=1:dx:15E3;
hackFactor=2;
DA=x.^hackFactor;
spatial={dx x DA} ;
%Timing
t_end=1E6;
timing=[t_end nan];
%Uplift Scen
uScen=3; maxElevation=1000;
upliftData={uScen maxElevation};
%Original bed and baselevel evolution
nb =length(x);
iniSurf = x.*0;
baseLevelDescent=0;
oriBed={iniSurf baseLevelDescent};
visibleFlag=0;
plotOut=2;
%Run the model
[z, dt]=SPLM(numM,parameters,spatial,timing,upliftData,oriBed,visibleFlag,plotOut);

