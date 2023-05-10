
% this runs all simulations needed for manuscript

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% baseline simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% varying salinity
for j= 0:50:200   %salinity in ppt or psu 


RunInfo.DateAndTime            = datestr(now);
RunInfo.InitialLakeDepthDim    = 10;            % Initial lake level in meters
RunInfo.plots                  = [];        % choose what plots to display (see table above)
RunInfo.PlotPeriod             = 100;           % The interval between plots in time steps - might need to change this
RunInfo.InitialrGuess          = 0.25;          % this is dimensional initial radius of channel in meters
RunInfo.Initialbeta_psu        = j;             % this is the intitial condition on salinity in the channel [psu = ppt]
RunInfo.VLi                    = 1e+06;         % volume coefficient in lake shape parameterisation
RunInfo.s0                     = 1000;          % channel length scale in meters
RunInfo.channel_geometry       = 1;             % 1 means circular, 1/2 means semi-circular
RunInfo.slope                  = 3;             % slope of channel and bed in degrees  
RunInfo.ice_thickness          = 100;           % ice thickness above channel

output = subglacial_brine_flow(RunInfo);

mkdir(['DIC_model/simulations/brine/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope)  '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/']);           
save(['DIC_model/simulations/brine/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope) '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/salinity=' num2str(RunInfo.Initialbeta_psu) '___.mat'],'output'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j= [0.1, 0.2, 0.3, 0.4, 0.5]   %salinity in ppt or psu 


RunInfo.DateAndTime            = datestr(now);
RunInfo.InitialLakeDepthDim    = 10;            % Initial lake level in meters
RunInfo.plots                  = [];        % choose what plots to display (see table above)
RunInfo.PlotPeriod             = 100;           % The interval between plots in time steps - might need to change this
RunInfo.InitialrGuess          = j;          % this is dimensional initial radius of channel in meters
RunInfo.Initialbeta_psu        = 100;             % this is the intitial condition on salinity in the channel [psu = ppt]
RunInfo.VLi                    = 1e+06;         % volume coefficient in lake shape parameterisation
RunInfo.s0                     = 1000;          % channel length scale in meters
RunInfo.channel_geometry       = 1;             % 1 means circular, 1/2 means semi-circular
RunInfo.slope                  = 3;             % slope of channel and bed in degrees  
RunInfo.ice_thickness          = 100;           % ice thickness above channel

output = subglacial_brine_flow(RunInfo);

mkdir(['DIC_model/simulations/brine/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope)  '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/']);           
save(['DIC_model/simulations/brine/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope) '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/salinity=' num2str(RunInfo.Initialbeta_psu) '___.mat'],'output'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=[0, 100]   %salinity in ppt or psu 


RunInfo.DateAndTime            = datestr(now);
RunInfo.InitialLakeDepthDim    = 10;            % Initial lake level in meters
RunInfo.plots                  = [];        % choose what plots to display (see table above)
RunInfo.PlotPeriod             = 100;           % The interval between plots in time steps - might need to change this
RunInfo.InitialrGuess          = 0.25;          % this is dimensional initial radius of channel in meters
RunInfo.Initialbeta_psu        = j;             % this is the intitial condition on salinity in the channel [psu = ppt]
RunInfo.VLi                    = 1e+06;         % volume coefficient in lake shape parameterisation
RunInfo.s0                     = 1000;          % channel length scale in meters
RunInfo.channel_geometry       = 1/2;             % 1 means circular, 1/2 means semi-circular
RunInfo.slope                  = 3;             % slope of channel and bed in degrees  
RunInfo.ice_thickness          = 100;           % ice thickness above channel

output = subglacial_brine_flow(RunInfo);

mkdir(['DIC_model/simulations/brine/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope)  '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/']);           
save(['DIC_model/simulations/brine/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope) '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/salinity=' num2str(RunInfo.Initialbeta_psu) '___.mat'],'output'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=[1e+05, 5e+06, 1e+07]   %salinity in ppt or psu 

RunInfo.DateAndTime            = datestr(now);
RunInfo.InitialLakeDepthDim    = 10;            % Initial lake level in meters
RunInfo.plots                  = [];        % choose what plots to display (see table above)
RunInfo.PlotPeriod             = 100;           % The interval between plots in time steps - might need to change this
RunInfo.InitialrGuess          = 0.25;          % this is dimensional initial radius of channel in meters
RunInfo.Initialbeta_psu        = 100;             % this is the intitial condition on salinity in the channel [psu = ppt]
RunInfo.VLi                    = j;         % volume coefficient in lake shape parameterisation
RunInfo.s0                     = 1000;          % channel length scale in meters
RunInfo.channel_geometry       = 1;             % 1 means circular, 1/2 means semi-circular
RunInfo.slope                  = 3;             % slope of channel and bed in degrees  
RunInfo.ice_thickness          = 100;           % ice thickness above channel

output = subglacial_brine_flow(RunInfo);

mkdir(['DIC_model/simulations/brine/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope)  '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/']);           
save(['DIC_model/simulations/brine/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope) '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/salinity=' num2str(RunInfo.Initialbeta_psu) '___.mat'],'output'); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=[2, 4]   %salinity in ppt or psu 
    for  i=[0, 100]

RunInfo.DateAndTime            = datestr(now);
RunInfo.InitialLakeDepthDim    = 10;            % Initial lake level in meters
RunInfo.plots                  = [];        % choose what plots to display (see table above)
RunInfo.PlotPeriod             = 100;           % The interval between plots in time steps - might need to change this
RunInfo.InitialrGuess          = 0.25;          % this is dimensional initial radius of channel in meters
RunInfo.Initialbeta_psu        = i;             % this is the intitial condition on salinity in the channel [psu = ppt]
RunInfo.VLi                    = 1e+06;         % volume coefficient in lake shape parameterisation
RunInfo.s0                     = 1000;          % channel length scale in meters
RunInfo.channel_geometry       = 1;             % 1 means circular, 1/2 means semi-circular
RunInfo.slope                  = j;             % slope of channel and bed in degrees  
RunInfo.ice_thickness          = 100;           % ice thickness above channel

output = subglacial_brine_flow(RunInfo);

mkdir(['DIC_model/simulations/brine/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope)  '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/']);           
save(['DIC_model/simulations/brine/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope) '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/salinity=' num2str(RunInfo.Initialbeta_psu) '___.mat'],'output'); 
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% water varying density 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 0:50:200   %salinity in ppt or psu 

RunInfo.DateAndTime            = datestr(now);
RunInfo.InitialLakeDepthDim    = 10;            % Initial lake level
RunInfo.plots                  = [];      % choose what plots to display (see table above)
RunInfo.PlotPeriod             = 100;             % The interval between plots in time steps - might need to change this
RunInfo.InitialrGuess          = 0.25;          % this is dimensional radius in meters
RunInfo.Initialbeta_psu        = i;              % this is the intitial condition on salinity in the channel [psu = ppt]
RunInfo.VLi                    = 1e+06;         % volume coefficient in lake shape parameterisation
RunInfo.s0                     = 1000;         % channel length scale, need to make the time step smaller if you make the channel length smaller
RunInfo.channel_geometry       = 1;             % 1 means circular, 1/2 means semi-circular
RunInfo.slope                  = 3;               % slope of channel and bed in degrees  
RunInfo.ice_thickness          = 100;             % ice thickness above channel

output = water_varying_density_subglacial_flow(RunInfo);

mkdir(['DIC_model/simulations/water_varying_density/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope)  '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/']);           
save(['DIC_model/simulations/water_varying_density/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope) '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/salinity=' num2str(RunInfo.Initialbeta_psu) '___.mat'],'output'); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% brine constant density 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j= 0:50:200   %salinity in ppt or psu 

RunInfo.DateAndTime            = datestr(now);
RunInfo.InitialLakeDepthDim    = 10;            % Initial lake level in meters
RunInfo.plots                  = [];        % choose what plots to display (see table above)
RunInfo.PlotPeriod             = 100;           % The interval between plots in time steps - might need to change this
RunInfo.InitialrGuess          = 0.25;          % this is dimensional initial radius of channel in meters
RunInfo.Initialbeta_psu        = j;             % this is the intitial condition on salinity in the channel [psu = ppt]
RunInfo.VLi                    = 1e+06;         % volume coefficient in lake shape parameterisation
RunInfo.s0                     = 1000;          % channel length scale in meters
RunInfo.channel_geometry       = 1;             % 1 means circular, 1/2 means semi-circular
RunInfo.slope                  = 3;             % slope of channel and bed in degrees  
RunInfo.ice_thickness          = 100;           % ice thickness above channel

output = density_subglacial_brine_flow(RunInfo);

mkdir(['DIC_model/simulations/brine_constant_density/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope)  '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/']);           
save(['DIC_model/simulations/brine_constant_density/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope) '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/salinity=' num2str(RunInfo.Initialbeta_psu) '___.mat'],'output'); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  outburst flood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 0:50:200   %salinity in ppt or psu 

RunInfo.DateAndTime            = datestr(now);
RunInfo.InitialLakeDepthDim    = 10;            % Initial lake level
RunInfo.plots                  = [];        % choose what plots to display (see table above)
RunInfo.PlotPeriod             = 100;             % The interval between plots in time steps - might need to change this
RunInfo.InitialrGuess          = 0.25;          % this is dimensional radius in meters
RunInfo.Initialbeta_psu        = i;              % this is the intitial condition on salinity in the channel [psu = ppt]
RunInfo.VLi                    = 1e+07;         % volume coefficient in lake shape parameterisation
RunInfo.s0                     = 1000;          % channel length scale, need to make the time step smaller if you make the channel length smaller
RunInfo.channel_geometry       = 1;             % 1 means circular, 1/2 means semi-circular
RunInfo.slope                  = 3;               % slope of channel and bed in degrees  
RunInfo.ice_thickness          = 100;             % ice thickness above channel

output = outburstflood_subglacial_flow(RunInfo);

mkdir(['DIC_model/simulations/outburst_flood/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope)  '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/']);           
save(['DIC_model/simulations/outburst_flood/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope) '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/salinity=' num2str(RunInfo.Initialbeta_psu) '___.mat'],'output'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 100 %0:50:200   %salinity in ppt or psu 

RunInfo.DateAndTime            = datestr(now);
RunInfo.InitialLakeDepthDim    = 10;            % Initial lake level
RunInfo.plots                  = [];      % choose what plots to display (see table above)
RunInfo.PlotPeriod             = 100;           % The interval between plots in time steps - might need to change this
RunInfo.InitialrGuess          = 0.25;          % this is dimensional initial radius of channel
RunInfo.Initialbeta_psu        = i;              % this is the intitial condition on salinity in the channel [psu = ppt]
RunInfo.VLi                    = 1e+06;         % volume coefficient in lake shape parameterisation
RunInfo.s0                     = 1000;         % channel length scale
RunInfo.channel_geometry       = 1;             % 1 means circular, 1/2 means semi-circular
RunInfo.slope                  = 3;               % slope of channel and bed in degrees  
RunInfo.ice_thickness          = 100;             % ice thickness above channel

output = temp_subglacial_brine_flow(RunInfo);

mkdir(['DIC_model/simulations/temp/H=' num2str(RunInfo.ice_thickness)  '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope)  '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/']);           
save(['DIC_model/simulations/temp/H=' num2str(RunInfo.ice_thickness) '/channel_geometry=' num2str(RunInfo.channel_geometry) '/slope=' num2str(RunInfo.slope) '/VLi=' num2str(RunInfo.VLi) '/InitialLakeDepth=' num2str(RunInfo.InitialLakeDepthDim) '/s0=' num2str(RunInfo.s0) '/radius=' num2str(RunInfo.InitialrGuess) '/salinity=' num2str(RunInfo.Initialbeta_psu) '___.mat'],'output', '-v7.3'); 
end