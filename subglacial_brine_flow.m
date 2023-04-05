function output = subglacial_brine_flow(RunInfo, varargin)

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                                          %%
%%         This function solves equations describing channelised subglacial brine flow      %%
%%  through cold ice coupled to a subglacial lake in terms of both pressure and discharge.  %%
%%                                                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% Function can be called by run_subglacial_brine_flow.m  
%
%
% Written by J. Kingslake, 2011, Department of Geography, University of
% Sheffield, UK.
% Edited by A. Jenson, 2022, Geophysical Insitute, University of Alaska
% Fairbanks, US.


% parse model input
InitialLakeDepthDim = RunInfo.InitialLakeDepthDim; 
plots = RunInfo.plots; % just specify which plots to display while this is running
PlotFreq = RunInfo.PlotPeriod; % how frequently to plot model output while the simulation is running
InitialrGuess = RunInfo.InitialrGuess; 
Initialbeta_psu = RunInfo.Initialbeta_psu;                                                                             
VLi = RunInfo.VLi;
s0 = RunInfo.s0;
channel_geometry = RunInfo.channel_geometry;
slope_degrees = RunInfo.slope;              
IceThickness = RunInfo.ice_thickness;

hLi = InitialLakeDepthDim; 
Slope = slope_degrees*pi/180;                  % slope of conduit bottom in radians 
InitialSGuess = pi*(InitialrGuess)^2;          % this is dimensional cross-sectional area
pL = 1;                                        % 1 means box-shaped lake, 2 means wedge, 3 means cone;   

%%%%%%%%%%%%%%%%%%%%%%%
Q_old = NaN;  %defining discharge for previous time step
hL_old = NaN; % defining lake depth for previous time step
P = 1;
Llo = 1;  
Lhi = 1;
%if ishandle(999) ==0
%    figure(999);
%    text(0.3,0.5,'Control-click to end simulation','FontSize',18)
%    set(999,'WindowButtonDownFcn',@EndRun);
%end

global UserReturn
UserReturn = 0;

LakeEmptied = 0;
ChannelClosed = 0;
Peak = NaN;
PeakTimeSecs = NaN;
PeakYearFrac = NaN;
Highstand = NaN;
Lowstand = NaN;

global paused
paused=0;
exitflag = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% MODEL PHYSICAL PARAMETERS %%%%%%%%%%%%
g=9.81;                                             % acceleration due to gravity
rho_i = 917;                                        % density of ice [kg m^-3]
L = 3.34e5;                                          % latent heat of fusion [J kg^-1]
ni = 0.0600;                                        % roughness of ice wall   [m^-1/3 s]
nb = 0.1629;                                        % roughness of bed material   [m^-1/3 s]

if channel_geometry == 1/2
    
%semicircular channel
n_full = ni * 0.611 + (1 - 0.611) * nb;             % roughness of channel when full  [m^-1/3 s]
f = 6.6 * n_full^2;                                 % friction factor for a semi-circular channel     

elseif channel_geometry == 1
%circular channel
n_full = ni * 0.611;                               % roughness of channel when full  [m^-1/3 s]
f = 5.4 * n_full^2;                                % friction factor for a circular channel    
end


Pi = g*rho_i*IceThickness;                           % function for ice-overburden pressure in N/m^2 or Pa 
p = Pi/100000;                                       % pressure in bars
sigma_i = 2093;                                     % specific heat of ice J/kg C


%Initial properties of brine

% initial density of brine as function of salinity
Initialrho_b  = (9.9780*10^(-10)*Initialbeta_psu^3 + 5.5328*10^(-08)*Initialbeta_psu^2 + 7.6346*10^(-04)*Initialbeta_psu + 9.9984*10^(-01))*(1000); % density of brine kg/m^3

% Convert psu to kg/m^3
Initialbeta = Initialbeta_psu*Initialrho_b/1000; % beta in kg/m^3 

%%%%%%%%%%%% this is for testing brine density!!!
%Initialrho_b  = 1000;
%%%%%%%%%%%%

% melting point of ice as function of salinity
theta_hat_salinity = -5.8202*10^(-07)*(Initialbeta_psu)^3 + 1.8653*10^(-06)*Initialbeta_psu^2 + -6.0536*10^(-02)*Initialbeta_psu + 2.5195*10^(-03);

% melting point of ice as function of pressure 
theta_hat_pressure = - 1.7628*10^(-26)*p^3 - 1.5226*10^(-16)*p^2 - 7.4477*10^(-8)*p + 1.1645*10^(-2); 

% melting point of ice as function of salinity and pressure 
Initialtheta_hat = theta_hat_salinity + theta_hat_pressure;

theta_b = Initialtheta_hat; 
theta_i = Initialtheta_hat;
theta0 = Initialtheta_hat;

%calculate flow law parameters/factors
n = 3;                                              % glens flow law exponent
A0 = 3.5*10^(-25);
R = 8.314;                                  % parameter needed for calculating A
T = theta_b + 273.15 + 7*10^(-8)*Pi;
T_star = 263 + 7*10^(-8);

% use constant A if near melting point
if theta_b >= 0
    A = 24*10^(-25);
elseif ((-10 + Pi*7*10^(-8))<theta_b) && (theta_b < 0)
    Q = 115;
    A = A0*exp((-Q/R)*(1/T - 1/T_star));
elseif theta_b < (-10 + Pi*7*10^(-8))
    Q = 6*10^4;
    A = A0*exp((-Q/R)*(1/T - 1/T_star));
end

K = 2* A* n^(-n);                                   % ice flow constant from Evatt (2006)

% Dimensional model scaling parameters
hL0  = hLi;                                         % lake depth scale [m] also flotation level
VL0  = (hL0/hLi)^pL * VLi;                          % lake volume scale [m^3]
QR0   = VLi/(10^5);                                        % discharge scale [m^3 s^-1]
psi0 = Initialrho_b*g*sin(Slope);                   % potential gradient scale   
SR0  = (f*Initialrho_b*g*QR0^2/psi0)^(3/8);        % area of channel scale [m^2]
m0   = psi0*QR0/L;                                  % melting of walls scale
t0   = rho_i*SR0/m0;                                % time scale [s]
N0   = (K * t0) ^(-1/3);                            % effective pressure scale [NR m^-2]

if channel_geometry==1
    hr0  = 2*sqrt((InitialSGuess/SR0)/pi);                              % channel roof height scale [m]
elseif channel_geometry ==1/2
    hr0  = sqrt((InitialSGuess/SR0)/(2*pi));                              % channel roof height scale [m]
end 

if Initialbeta_psu==0
    beta0=1;
else
    beta0 = Initialbeta;   % brine concentration scale [kg m^(-3)]
end

% dimensionless parameters of model              
zeta = t0 * (hLi^pL) * QR0/(pL * VLi * hL0^pL);
delta = N0/(s0*psi0);
gamma = theta0 * sigma_i /L; 
lambda = SR0*s0/(t0*QR0);

%%%%%%%%%%%%%%% SET UP TIME %%%%%%%%%%%%

% set up space grid
S_end = 1;
ds= 1/50;                                           % space step    usually 0.01
s = 0:ds:S_end;                                     % space vector
Ls = length(s);                                     % number of space steps

dt= ds/(1000);                                      % time step default 0.01 (0.01 in Kingslake)
T = 100;                                            % dimensionless simulation time default 500
t=0:dt:T;                                           % time vector
Lt = length(t);                                     % number of time steps
TSamp = t(1:10:Lt);                                 % set up a sampling time vector
TDays =TSamp*t0/3600/24;
tDays = t*t0/3600/24;


%%%%%%%%%%% PREALLOCATE ARRAYS %%%%%%%%
disp('Pre-allocating Sampling Arrays....')

hLSamp = zeros(length(TSamp),1);                    % array for sampling values of lake depth [m]
hrLakeSamp = zeros(length(TSamp),1);                % array for sampling values of channel roof height [m]
QRLakeSamp = zeros(length(TSamp),1);                % array for sampling values of the discharge at the lake [QR0]
QRendSamp  = zeros(length(TSamp),1);                % array for sampling values of the discharge at the end of the tunnel [QR0]
QRSamp = NaN(length(TSamp),Ls);                     % array for sampling values of the channel discharge profile
SRSamp = NaN(length(TSamp),Ls);                     % array for sampling values of the channel area profile
NRSamp = NaN(length(TSamp),Ls);  
betaSamp = NaN(length(TSamp),Ls);                 % array for sampling values of the brine concentration
beta_psuSamp = NaN(length(TSamp),Ls);              % array for sampling values
theta_hatSamp = NaN(length(TSamp),Ls);             % array for sampling values of the salinity-dependent melting point of ice
rho_bSamp = NaN(length(TSamp),Ls);  
psiSamp = NaN(length(TSamp),Ls);

disp('Done')

disp('Pre-allocating Variable Arrays....')
% Lake variables
hL  = zeros(2,1);                                   % lake level  [m]

% R Channel variables
hr = zeros(2,Ls);                                   % channel roof height [m]
SR = zeros(2,Ls);                                   % channel cross-section
NR = zeros(2,Ls);                                   % channel effective pressure
QR  = zeros(2,Ls);                                  % discharge in channel  [m^3 s^-1]
beta = zeros(2, Ls);                                % brine concentration
beta_psu = zeros(2,Ls); 
theta_hat = zeros(2,Ls);                            % salinity-dependent melting point  
psi = zeros(2,Ls);   

disp('Done')

% Sampling Counter
SampNumber = 2;
SpinUpFin = 1;

%%%%%%% initial conditions %%%%%%%%%

disp('Starting Initial Conditions...')


%%% define hydraulic potential
Initialpsi= Initialrho_b*g*sin(Slope)/psi0;         % idealised dimensionless hydraulic gradient

SR_temp(1,:) = ones(1,Ls).*(InitialSGuess/SR0);
hL(1,1) = InitialLakeDepthDim/hL0;

% Define N at top end
NR(1,1) = 0; % BC for effective pressure at lake
%NR(1,Ls) = 0; % BC for effective pressure at channel end

beta_temp(1,:) = ones(1,Ls).*(Initialbeta/beta0);                  % this defines the IC for beta
beta_psu_temp(1,1:Ls) = ones(1,Ls).*(Initialbeta_psu);
theta_hat_temp(1,:) = ones(1,Ls).*(Initialtheta_hat/theta0);     
rho_b_temp(1,:)= ones(1,Ls).*Initialrho_b;
psi_temp = ones(1,Ls).*Initialpsi;

 % Boundary Layer Method in section 2.2.5 of Kingslake (2013)
        % initial conditions

        % S at far end of channel fixes Q there
        QR_temp(1,Ls) = sqrt((SR_temp(1,Ls)^(8/3))*psi_temp(1,Ls));
        % fixes Q for the whole channel
        QR_temp(1,:) = QR_temp(1,Ls);
        
        %initial N Guess is important 
        %NR_temp = NR(1,1)*(1-s./1); %linear gradient from BC at the lake 
        %NR_temp = NR(1,1);
        %NR_temp(1,1) = NR(1,1);
        NR_temp = -s;
 
           %for k = 1:Ls-1
               %NR_temp(1,k+1) = NR_temp(1,k) + ds/delta * ( QR_temp(1,k)*abs(QR_temp(1,k))/(SR_temp(1,k).^(8/3)) - psi_temp(1,k));
           %end
   
disp('...')

QR(1,:) = QR_temp;
NR(1,:) = NR_temp;
SR(1,:) = SR_temp;
beta(1,:) = beta_temp;
beta_psu(1,:) = beta_psu_temp;
theta_hat(1,:) = theta_hat_temp;    
rho_b(1,:) = rho_b_temp;
psi(1,:) = psi_temp;

disp('Done')
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% MAIN LOOP %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Starting Main Loop...')

for i = 2:Lt

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Step Channel Crossection Forward %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update the R channel cross-sectional area

    SR(2,:) = SR(1,:)+ dt*(abs((QR(1,:).^3)) ./ ((SR(1,:).^(8/3)).*(1+gamma.*(theta_hat(1,:) - theta_i/theta0))) - SR(1,:).*(NR(1,:)).^3);

    if any(SR(2,:)<=0)
        disp 'SR has gone to zero - reducing the time step can help with this.'
        
        ChannelClosed = 1;
        return
    end 


% calculate non-dimensional channel height from this
    if channel_geometry ==1
        hr(2,:) = 2*sqrt(SR(2,:)/pi);
    elseif channel_geometry==1/2
        hr(2,:) = sqrt(SR(2,:)/(2*pi));
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Step Lake and lake effective pressure forward %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    hL(2,1) = max(0,(hL(1,1) + (dt * zeta/(hL(1,1)^(pL-1)) * (-QR(1,1)))));

    if hL(2,1)*hL0 <= hr(2,:)*hr0
        
        disp 'Lake level dropped below the height of the channel roof! - This version of the model cannot deal with this. A later version described in Ch. 6 of Kingslake (2013) [thesis] implements an open-channel model to simulate what happens when the lake level drops below the level of the channel roof. '
        LakeEmptied = 1;
        return
        %             error 'Lake Emptied!!!!'
    end
    
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Find NR and QR which fit this channel shape and BC's %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     % Boundary Layer Method

            QR(2,:) = sqrt(psi(1,Ls)*(SR(2,Ls)^(8/3)));

      % for dirichlet BC at start of channel 
               NR(2,1) = 0;
      % for dirichlet BC at end of channel
               %NR(2,Ls) = 0;
      %for neumann BC at start of channel 
              % NR(2,1) = NR(2,2);
      %for neumann BC at end of channel 
               NR(2,Ls) = NR(2,Ls-1);

        for k = 1:Ls-1
             NR(2,k+1) = NR(2,k) + (ds/delta)*(QR(2,k)*abs(QR(2,k))/(SR(2,k)^(8/3)) - psi(1,k));
        end     

        %NR(2,:) = 0;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%% Find beta which fit BC's %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Upwind difference scheme to solve brine equation 
    
       % define boundary condition at the lake 
       beta(2,1)= Initialbeta/beta0;
    
       beta(2,2:Ls) = beta(1,2:Ls) - dt*((beta(1,2:Ls) ./ SR(2,2:Ls)).* ((SR(2,2:Ls) - SR(1,2:Ls))/dt) ...
              + (QR(2,2:Ls)./(lambda * SR(2,2:Ls))).*((beta(1,2:Ls)-beta(1,1:Ls-1))/ds));
    
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % used for testing influence of brine but can be left
    %rho_b(1,:)  = (9.9780*10^(-10)*beta_psu(1,:).^3 + 5.5328*10^(-08)*beta_psu(1,:).^2 + 7.6346*10^(-04)*beta_psu(1,:) + 9.9984*10^(-01))*(1000); % density of brine kg/m^3

    
    %convert brine from kg/m^3 to psu 
    beta_psu(2,:) = ((beta(2,:)*beta0)*1000./rho_b(1,:));      % psu  - Rutishauser (140-160 psu)    dimensional!                                                     
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%% Find new salinity-dependent melting point %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % melting point of ice as function of salinity and pressure 
    theta_hat(2,:) =  (theta_hat_pressure + (-5.8202*10^(-07)*(beta_psu(2,:)).^3 + 1.8653*10^(-06)*(beta_psu(2,:)).^2 - 6.0536*10^(-2)*(beta_psu(2,:)) + 2.5195*10^(-3)))/theta0;


    % density of brine as function of salinity and pressure 
    rho_b(2,:)  = (9.9780*10^(-10)*beta_psu(2,:).^3 + 5.5328*10^(-08)*beta_psu(2,:).^2 + 7.6346*10^(-04)*beta_psu(2,:) + 9.9984*10^(-01))*(1000); % density of brine kg/m^3

    %%%%%%%%%%%%%%%%%this is for testing denisty!!!
    %rho_b(2,:) = 1000;
     %%%%%%%%%%%%%%%%

    psi(2,:) = rho_b(2,:)*g*sin(Slope)/psi0;         % idealised dimensionless hydraulic gradient


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% Plotting %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if rem(i,PlotFreq)==0
        

  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%% Calculate day, hours, min, secs. and display %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %sprintf('%3.3f %% Complete  ',t(i)/T * 100)
        %TotalSeconds = t(i)*t0;
        %TotalDays = TotalSeconds/(24*3600);
        %WholeDays = fix(TotalDays)
        %         disp(WholeDays)
        %disp(datestr(datenum(0,0,0,0,0,TotalSeconds),'HH:MM:SS'))
        %         disp(t(i))
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%           Plotting           %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       
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
%      14        Along-channel profile of brine concentration, beta(s)
%      15        Along-channel profile of salinity-dependent melting point, theta_hat(s)
% 
        
        
        
        % QR profile   
        if nnz(plots==1)~=0
            if ishandle(1) ==0; figure(1);set(1,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',1)
            plot(s*s0,QR(2,:)*QR0,'b')
            title 'Discharge Along-Channel Profiles'
            pause(0.000001)
        end
        
        % effective pressure profile
        if nnz(plots==2)~=0
            if ishandle(2) ==0; figure(2);set(2,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',2)
            plot(s*s0,NR(2,:)*N0,'b')
            %       axis([0 1 -3 1])
            title 'Effective Pressure Along-Channel Profiles'
            pause(0.000001)
        end

        % QR time series
        if nnz(plots==3)~=0
            if ishandle(3) ==0; figure(3);set(3,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',3)
            hold on
            plot(tDays(i),QR(2,1)*QR0,'.r','MarkerSize',5)
            title 'QR time series (at lake)'
            pause(0.000001)
        end

        % SR along channel profile 
        if nnz(plots==4)~=0
            if ishandle(4) ==0; figure(4);set(4,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',4)
            hold off
            plot(s*s0,SR(2,:)*SR0,'b')
            %             hold on
            %             plot(s,SC(2,:)*SC0,'r')
            title 'SR Along-Channel Profiles '
            pause(0.000001)
        end
        
        % NR time series at terminus
        if nnz(plots==5)~=0
            if ishandle(5) ==0; figure(5);set(5,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',5)
            hold on
            if i>1
            plot(tDays(i),NR(2,end)*N0,'r.')
            end
            title 'N at the terminus'
            pause(0.000001)
        end

        %  Time series of the minumum area of the channel
        if nnz(plots==6)~=0
            if ishandle(6) ==0; figure(6);set(6,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',6)
            hold on
            if i>1
            plot(tDays(i),min(SR(2,:)*SR0),'r.')
            end
            title ' Time series of the minumum area of the channel, min(SR)'
            pause(0.000001)
        end

        % Along-channel profile of disharge and effective pressure together
        % (QR, NR and QC and NC profiles)
        if nnz(plots==7)~=0
            if ishandle(7) ==0; figure(7);set(7,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',7)
            hold off
            plot(s*s0,QR(2,:)*QR0,'b',s*s0,NR(2,:)*N0,'--b')
            hold on
%             Line2 = line([XR XR],[-1 1]);
%             set(Line2,'Color','b')
            axis([0 s0 -0.1 1])
            title 'Profiles of Q and N'
            pause(0.000001)
        end

        % 3D phase-space plot of lake depth, discharge at the lake and channel crossectional are at the lake (useful!)  
        % QR-hL-SR 3d phase space
        if nnz(plots==8)~=0
            if ishandle(8) ==0; figure(8);set(8,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',8)
            hold on
            plot3(hL(2,1)*hL0,QR(2,1)*QR0,SR(2,1)*SR0,'.r','MarkerSize',5)
            xlabel hL
            ylabel QR
            zlabel SR
            grid on
            pause(0.000001)
        end

        % Rate of change of discharge at the lake, dQR/dt
        if nnz(plots==9)~=0
            if ishandle(9) ==0; figure(9);set(9,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',9)
            hold on
            if i>1
            plot(t(i),((QR(2,1)-QR(1,1))*QR0)/(tDays(i)-tDays(i-1)),'r.')
            end
            title 'dQR/dt'
            pause(0.000001)
        end
        
% Along-channel profile of disharge, effective pressure and discharge together
        if nnz(plots==10)~=0
            if ishandle(10) ==0; figure(10);set(10,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',10)          
            plot(s*s0,NR(2,:)*N0,'r',s*s0,SR(2,:)*SR0,s*s0,QR(2,:)*QR0)          
            title 'NR, SR and QR profiles'
            hold off
           pause(0.000001)
        end
     %   Rate of change of discharge with lake hieght, dQR/dh 
        if nnz(plots==11)~=0
            if ishandle(11) ==0; figure(11);set(11,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',11)
            hold on
            plot(t(i),((QR(2,1)-QR(1,1))*QR0)/((hL(2,1)-hL(1,1))*hL0),'r.')
            title 'dQR/dh'
           pause(0.000001)
        end
 
        % Channel radius profile
        if nnz(plots==12)~=0
            if ishandle(12) ==0; figure(12);set(12,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',12)
            hold off
            plot(s,sqrt(SR(2,:)*SR0/pi)) 
            title 'Channel radius profile'
            pause(0.000001)
        end

        
        % brine profile
        if nnz(plots==14)~=0
            if ishandle(14) ==0; figure(14);set(14,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',14)
            plot(s*s0,beta_psu(2,1:Ls),'o')
            %       axis([0 1 -3 1])
            title 'Brine concentration along-channel profile'
            pause(10)
        end
 
         % melting point profile
        if nnz(plots==15)~=0
            if ishandle(15) ==0; figure(15);set(15,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',15)
            plot(s*s0,theta_hat(2,1:Ls)*theta0,'b')
            %       axis([0 1 -3 1])
            title 'Ice-brine interface melting point along channel'
           pause(0.000001)
        end
        
             % Average Velocity along channel profile 
        if nnz(plots==16)~=0
            if ishandle(16) ==0; figure(16);set(16,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',16)
            hold off
            plot(s*s0,(QR(2,:)*QR0)./(SR(2,:)*SR0),'b')
            %             hold on
            %             plot(s,SC(2,:)*SC0,'r')
            title 'Average Velocity Along-Channel Profiles '
            pause(0.000001)
        end


    end

    %%%%%%%%%%%%%%%%%%
    %%%  SAMPLING  %%%
    %%%%%%%%%%%%%%%%%%

    if TSamp(SampNumber)==t(i)
        hrLakeSamp(SampNumber,1) = hr(2,1);
        QRLakeSamp(SampNumber,1) = QR(2,1);
        QRendSamp(SampNumber,1) = QR(2,end);
        hLSamp(SampNumber,1) = hL(2,1);
        QRSamp(SampNumber,:) = QR(2,:);
        SRSamp(SampNumber,:) = SR(2,:);
        NRSamp(SampNumber,:) = NR(2,:);
        betaSamp(SampNumber,:) = beta(2,:);
        theta_hatSamp(SampNumber,:) = theta_hat(2,:);
        beta_psuSamp(SampNumber,:) = beta_psu(2,:);
        rho_bSamp(SampNumber,:) = rho_b(2,:);
        psiSamp(SampNumber,:) = psi(2,:);

        SampNumber = SampNumber + 1;

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Detect Peak of floods %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if QR(2,1) <QR(1,1) && QR(1,1)>Q_old
    Peak(P) = QR(1,1);
    PeakTimeSecs(P) = t(i-1)*t0;          % the seconds since the beginning
    PeakTimeDays(P) = PeakTimeSecs(P)/3600/24;
    PeakYearFrac(P) = PeakTimeDays(P)/365-floor(PeakTimeDays(P)/365);
    P=P+1;
    
    
    
end
Q_old = QR(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Detect Lake Highsatnd %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hL(2,1) <hL(1,1) && hL(1,1)>hL_old
    Highstand(Lhi) = hL(1,1);
    HighstandTimeSecs(Lhi) = t(i-1)*t0;          % the seconds since the beginning
    HighstandTimeDays(Lhi) = HighstandTimeSecs(Lhi)/3600/24;
    HighstandYearFrac(Lhi) = HighstandTimeDays(Lhi)/365-floor(HighstandTimeDays(Lhi)/365);
    Lhi=Lhi+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Detect Lake Lowsatnd %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hL(2,1) >hL(1,1) && hL(1,1)<hL_old
    Lowstand(Llo) = hL(1,1);
    LowstandTimeSecs(Llo) = t(i-1)*t0;          % the seconds since the beginning
    LowstandTimeDays(Llo) = LowstandTimeSecs(Llo)/3600/24;
    LowstandYearFrac(Llo) = LowstandTimeDays(Llo)/365-floor(LowstandTimeDays(Llo)/365);
    Llo=Llo+1;
end
hL_old = hL(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Detect Convergence to a limit cycle %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Converged = P>=4 && Lhi>=4 && Llo>=4 ...
    && abs((1-Peak(P-1)/Peak(P-2))) <0.0005 && abs((1-Peak(P-1)/Peak(P-3))) <0.0005...
    && abs((1-Highstand(Lhi-1)/Highstand(Lhi-2))) <0.0005 && abs((1-Highstand(Lhi-1)/Highstand(Lhi-3))) <0.0005...
    && abs((1-Lowstand(Llo-1)/Lowstand(Llo-2))) <0.0005 && abs((1-Lowstand(Llo-1)/Lowstand(Llo-3))) <0.0005;
if Converged
    disp(['Limit Cycle Reached - change line code around line 640 in ' mfilename  '.m to change this behaivour'])

    return

end


%%% create output variable %%%

output.QR0 = QR0;
output.hL0 = hL0;
output.VL0 = VL0;
output.VLi = VLi;
output.SR0 = SR0;
output.t = tDays;
output.QR = QRSamp;
output.QRLake = QRLakeSamp;
output.QRend = QRendSamp;
output.hL = hLSamp;
output.SR = SRSamp;
output.hrLake = hrLakeSamp;
output.hr0 = hr0;
output.Ls = Ls;
output.beta0 = beta0;
output.beta = betaSamp;
output.beta_psu = beta_psuSamp;
output.theta_hat = theta_hatSamp;
output.theta0 = theta0;
output.rho_b = rho_bSamp;
output.T = TSamp;
output.TDays = TDays;
output.t0 = t0;
output.s0 = s0;
output.dt = dt;
output.ds = ds;
output.s = s;
output.psi = psiSamp;
output.NR = NRSamp;
output.N0 = N0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Loop values round %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % copy values into position 1

    hL(1,1) = hL(2,1);
    hr(1,:) = hr(2,:);
    SR(1,:)  = SR(2,:);
    QR(1,:) = QR(2,:);
    NR(1,:) = NR(2,:);
    beta(1,:) = beta(2,:);
    theta_hat(1,:) = theta_hat(2,:);
    beta_psu(1,:) = beta_psu(2,:);
    rho_b(1,:) = rho_b(2,:);
    psi(1,:) = psi(2,:); 

    % wipe old values
    hL(2,1) = 0;
    hr(2,:) = 0;
    SR(2,:)  = 0;
    QR(2,:) = 0;
    NR(2,:) = 0;
    beta(2,:) = 0;
    theta_hat(2,:) = 0;
    beta_psu(2,:) = 0;
    rho_b(2,:) = 0;
    psi(2,:) = 0;

     %%%%%%%%%%%%%%%%%%%%
     %%%  PAUSE CODE  %%%
     %%%%%%%%%%%%%%%%%%%%
 
     
     if UserReturn == 1
         button = questdlg('Return?');
         if strcmp(button,'Yes')
             return
         else
             UserReturn = 0;
         end
     end

    
end
disp('Model reached t=T')
exitflag = 1
end
 %%% Pause function %%%
 
 function [paused] = PauseSim(src,evt);
 
 global paused
 if paused == false;
     paused =true;
     disp('Paused')
     return
 end
 if paused == true;
     paused =false;
     disp('UnPaused')
     return
 end
 end
 
 function [UserReturn] = EndRun(src,evnt)
    global UserReturn
 if strcmp(get(src,'SelectionType'),'alt')
       disp Return?
       UserReturn = 1;
    else
       disp('Use control-click return')
    end
 end