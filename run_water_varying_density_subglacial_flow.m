 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%                                                       %%
 %%                 This script runs                      %%
 %%       water_varying_density_subglacial_flow.m         %%
 %%                                                       %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This model code is from Kingslake, J. and has been edited extensivesly by Jenson, A.

% This script is desinged to run the function density_FullNyeFowler.m
% 
%
%Uncommenting the save call at the bottom of the script enables you to save the results of each simulation.  
%
% See Appendix of Jenson et al. (20XX) for a decription of the numerical
% methods. 
%
% The parameter RunInfo.plots allows the user to choose which plots they
% want to see as the simulation runs. Define it as an array with either
% none, one or multiple entries:
%    
%      Entry     Plot
%      1         Along-channel profile of discharge, Q(s) 
%      2         Along-channel profile of effective pressure, N(s) 
%      3         Time series of discharge at the lake, Q(0,t)
%      4         Along-channel profile of channel cross-sectional area, S(s)
%      5         Time series of effective at the terminus, N(s0,t) 
%      6         Time series of the minumum area of the channel,    min(S(s))(t)
%      7         Along-channel profile of disharge and effective pressure together
%      8         3D phase-space plot of lake depth, discharge at the lake and channel crossectional are at the lake (useful!) 
%      9         Rate of change of discharge at the lake, dQR/dt
%      10        Along-channel profile of disharge, effective pressure and discharge together
%      11        Rate of change of discharge with lake hieght, dQR/dh 
%      12        Along-channel profile of channel radius at the lake, hr(s) 
%      13        Time series of lake input, Qin(t) 
%      14        Along-channel profile of brine concentration, beta(s)
%      15        Along-channel profile of salinity-dependent melting point, theta_hat(s)
%      16        Velocity at the lake

% e.g. if you want to plot Q(s,t) and N(s0,t) use:
% RunInfo.plots    = [1 5];

%
% Written by J. Kingslake, 2011, Department of Geography, University of
% Sheffield, UK.
% Edited by A. Jenson, 2022, Geophysical Institute, University of Alaska
% Fairbanks, US. 


% varying salinity
for i = 0:50:200   %salinity in ppt or psu 

RunInfo.DateAndTime            = datestr(now);
RunInfo.InitialLakeDepthDim    = 10;            % Initial lake level
RunInfo.plots                  = [];      % choose what plots to display (see table above)
RunInfo.PlotPeriod             = 100;             % The interval between plots in time steps - might need to change this
RunInfo.InitialrGuess          = 0.25;          % this is dimensional radius in meters
RunInfo.Initialbeta_psu        = i;              % this is the intitial condition on salinity in the channel [psu = ppt]
RunInfo.VLi                    = 1e+07;         % volume coefficient in lake shape parameterisation
RunInfo.s0                     = 1000;         % channel length scale, need to make the time step smaller if you make the channel length smaller
RunInfo.channel_geometry       = 1;             % 1 means circular, 1/2 means semi-circular
RunInfo.slope                  = 3;               % slope of channel and bed in degrees  
RunInfo.ice_thickness          = 100;             % ice thickness above channel

output = water_varying_density_subglacial_flow(RunInfo);


mkdir(['DIC_model/simulations/water_varying_density/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope)  '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/']);           
save(['DIC_model/simulations/water_varying_density/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope) '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/salinity=' num2str(RunInfo.Initialbeta_psu) '___.mat'],'output'); 
end