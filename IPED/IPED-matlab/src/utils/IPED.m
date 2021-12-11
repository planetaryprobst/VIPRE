function IPED(spicepath,body)

%--------------------------------------------------------------------------
%   Main script for computations called by setup_IPED.m. Computes entry 
%   conditions for all safe entry opportunities for one trajectory to the 
%   input 'body'. 
%
%   Returns plots for visualization.
% 
%--------------------------------------------------------------------------
%   Form:
%   IPED(spicepath,body)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   spicepath   str         path to the generic SPICE kernels folder
%   body        str         body name
%
%   ------
%   Output
%   ------
%   
%   plots
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 10.11.2021 |  A. Probst    | First release
%*************************************************************************%

disp('#*************************************************************************#')
disp('#')
disp('#    Welcome to IPED')
disp('#')
disp('#*************************************************************************#')
disp('#')


%% LOAD STANDARD SPICE KERNELS

success = loadKernels(spicepath);

if success
    disp('#    The SPICE kernels have been successfully loaded ...')
    disp('#')
    disp('#*************************************************************************#')
    disp('#')
end

%% Configurations

% Planet Data

    % reference frame definitions
    from = 'ECLIPJ2000';
    to = ['IAU_',body];

    % Creating key-value set for planet's entry altitude, in km
    keySet = {'Mars','Saturn','Uranus','Neptune'};
    valueSet = [700 700 700 700];
    M = containers.Map(keySet,valueSet);

    % planet data retrieval
    [mu,radii,~,rings] = planetData(body);
    % trajectory data retrieval
    [trajectories,vInf,epoch] = trajData(body);

% Loop Data
    % number of computed states between initial and periapsis
    n = 1e3;

    % clock angle B-plane interval
    theta_min = -cspice_pi;
    theta_delta = cspice_pi/(2 * 90 - 1);
    theta_max = cspice_pi;
    
    % B-vector length interval step
    Babs_delta = radii(1)/25; % km
    
% Graphics Data
    % figure counter
    fig = 1;

    % number of patches of ellipsoid
    nrFacets = 150;

    % subplot definition
    row = 4;
    column = 2;
    plotNr = 1;

%% Entry Definition

% entry altitude definition, in km
hEntry = M(body);

% % % % % % % 
trajectory = trajectories{1,:};
v_inf = vInf(:,1);
epoch = epoch(1);

% Body-fixed coordinate system
xPlanet = [1;0;0];
yPlanet = [0;1;0];
zPlanet = [0;0;1];

disp('#    You are computing the following case:')
disp('#')
disp(['#    body:           ',body])
disp(['#    entry altitude: ',num2str(hEntry),' km'])
disp('#')
disp('#*************************************************************************#')
disp('#')
disp('#    Status Output:')
disp('#')



%% Variable Definition

% ephemeris time, in s past J2000 epoch
et = epoch * cspice_spd;

% rotation matrix for vInfinity transformation             
rotate = cspice_pxform( from, to, et );
v_inf =  rotate * v_inf;

%v_inf2 = bsxfun(@(x,y) x * y,rotate,reshape(vInf,3,1,length(vInf)));

% B-vector length interval step
Babs_min = 0; % 10; 
Babs_max = 8 * (radii(1)+radii(2))/2;

% Color Maps
%CMap_generic = NaN(nrFacets + 1,nrFacets + 1,1);
CMap_phiMin = NaN(nrFacets + 1,nrFacets + 1,1);
CMap_phiMax = NaN(nrFacets + 1,nrFacets + 1,1);
CMap_vEntryMin = NaN(nrFacets + 1,nrFacets + 1,1);
CMap_vEntryMax = NaN(nrFacets + 1,nrFacets + 1,1);
CMap_vRelEntryMin = NaN(nrFacets + 1,nrFacets + 1,1);
CMap_vRelEntryMax = NaN(nrFacets + 1,nrFacets + 1,1);
CMap_vrot = NaN(nrFacets + 1,nrFacets + 1,1);


%% Create the Planet's Sphere 
% as a preparation for the plots

% ellipsoidal references frame axes, body fixed
u = [1,0,0;0,1,0;0,0,1];

% ellipsoid point matrices
[x,y,z]=Ellipsoid(radii,u,nrFacets);
% % atmosphere point matrices
% [aX,aY,aZ]=Ellipsoid(radii + hEntry,u,nrFacets/4);

% ellipsoid point vector
[a,b] = size(x);
j = 0;
P = zeros(a * b,5);

for k    = 1:b 
    for i = 1:a
        j = j + 1; 
        P(j,1) = x(i,k); 
        P(j,2) = y(i,k);
        P(j,3) = z(i,k);
        P(j,4) = i;
        P(j,5) = k;
    end
end
  
    
%% Planet's Rotation omega
% in rad/s
omega = planetRotation(body);
        
%% Atmospheric Speed depending on Latitude 
% rotational speed for latitude range [-90 deg:90 deg] in m/s

lat = [-cspice_pi/2:cspice_pi/(2 * 89):cspice_pi/2]';
[vRot,~] = rotationalSpeed(lat,omega,radii+hEntry);


%% Varying the Incoming Trajectory
%[~,~,B_abs1,theta] = RV2B_plane(rInit,vInit,mu);
i = 1;
trajStates = [];
v_rot_entry = [];
state_entry = [];
states_hazard = [];
reclat = [];
tic

% status print: number of current iteration
disp(['#    Step ',num2str(i),'/360']);

% loop over B_theta
for B_theta = theta_min:theta_delta:theta_max
    
    % Creating B-plane
    [~,~,~,~,h_uv,~,a] = VinfThe2B_plane(v_inf,B_theta,mu);
    
    % loop over length of B_vec
    for Babs = Babs_min:Babs_delta:Babs_max
        
        if Babs == 0
            [radius, lon0, lat0] = cspice_reclat(-v_inf);
            reclat = [reclat,[radius;rad2deg(lon0);rad2deg(lat0)]];
            Babs_min = Babs_min + Babs_delta;
            continue
        end
        
        % Calculating Orbit Parameters:
        
%         % B-vector, in km
%         B_vec = Babs * B_uv;
        % eccentricity, -
        e_abs = sqrt((Babs/a)^2 + 1);
        % angular momentum, in km^2/s^2
        h = sqrt(mu * a * (1 - e_abs^2)) * h_uv;
        % eccentricity vector, -
        e = cross(v_inf,h)/mu + v_inf/vecnorm(v_inf,2,1);
        % Radius @ periapsis in magnitude and vector, km
        rP_abs = a * (1 - e_abs);
        rP = rP_abs * e/vecnorm(e,2,1);
        % Velocity @ periapsis in magnitude and vector, km/2
        vP_abs = sqrt(mu * (2/rP_abs - 1/a));
        vP = vP_abs * cross(h_uv,rP/vecnorm(rP,2,1),1);
        
        % Orbit state in XYZ-coordinates
        state = [rP;vP];
        % State in conical elements
        elts0 = cspice_oscelt(state,et,mu);        
        
        % Computing the orbital states until periapsis
        [states,elts] = computingOrbStates(elts0,n,true,et);
        
        % Locating point of closest approach, or, if given, closest point
        % to entry height
        [altitude,~,idx,flag] = ...
            closestApproach(states(1:3,:),radii,hEntry);
            
        
        % If trajectory is an entry trajectory    
        if flag
            
            % Check if it is a hazardous trajectory or not
            [state_hazard,safe] = hazardCrssng(states,elts,et,radii,rings);
            states_hazard = [states_hazard,[state_hazard;safe]];
            
            if safe
                % Calculating the Flight Path Angle at entry / closest approach

                % determining entry point of trajectory
                diff = altitude(1,idx) - hEntry;
                if diff > 0.1
                    % call function to determine entry point
                    [altitude,~,~,states,idx,diff] = ...
                                  entryPoint(elts,idx,radii,hEntry,10,et,diff);

                    if altitude(1,idx) > 700.1
                        disp(['err in ',num2str(Babs),', diff = ',...
                            num2str(diff),', idx = ',num2str(idx)])
                    end
    %                 else
    %                     % adapt the states vector in length
    %                     states = states(:,idx-4:end);
    %                     idx = 5;
    % 
    %                     if size(states,2) > 10
    %                         states = states(:,1:10);
    %                     else
    %                         l = 10 - size(states,2);
    %                         states = [states,NaN(6,l)];
    %                     end

                    % collecting trajectory states before entry for test plotting 
                    trajStates = [trajStates,states(1:3,idx)];

                end

                % Calculating the lon/lat of the trajectory radius vector
                [radius,lon_idx,lat_idx] = cspice_reclat(states(1:3,idx));
                
                reclat = [reclat,[radius;rad2deg(lon_idx);rad2deg(lat_idx)]];
                
                % Rotational velocity of planet at point idx
                v_rot_idx = rotationalSpeed(lat_idx,omega,radii+hEntry);
                v_rot_x = - v_rot_idx * sin(lon_idx);
                v_rot_y = v_rot_idx * cos(lon_idx); 
                % vectorized
                v_rot = [v_rot_x;v_rot_y;0];

                phi_idx = RV2FPA(states(1:3,idx),states(4:6,idx),mu);
                % entry velocity including atmospheric rotation
                v_entry_idx = states(4:6,idx) - v_rot;

                % testing variables
                v_rot_entry = [v_rot_entry,[v_rot;B_theta]]; 
                state_entry = [state_entry,states(:,idx)];

                % identify point on sphere
                k = dsearchn(P(:,1:3),states(1:3,idx)');

                CMap_temp = NaN(nrFacets + 1,nrFacets + 1,1);

                % allocate values to identified point(s)
                % flight path angle
                CMap_temp(k) = rad2deg(phi_idx);
                CMap_phiMin = min(CMap_temp,CMap_phiMin);
                CMap_phiMax = max(CMap_temp,CMap_phiMax);
                % orbital velicity @ entry
                CMap_temp(k) = vecnorm(states(4:6,idx),2);
                CMap_vEntryMin = min(CMap_temp,CMap_vEntryMin);
                CMap_vEntryMax = max(CMap_temp,CMap_vEntryMax);
                % relative entry velocity 
                CMap_temp(k) = vecnorm(v_entry_idx,2);
                CMap_vRelEntryMin = min(CMap_temp,CMap_vRelEntryMin);
                CMap_vRelEntryMax = max(CMap_temp,CMap_vRelEntryMax);
                % rotational velocity
                CMap_temp(k) = vecnorm(v_rot,2) * 1e3;
                CMap_vrot = max(CMap_temp,CMap_vrot);
                
            end

        else
            % FPA for closest approach (if no entry) is periapsis, 
            % therefore FPA ~ 0
            continue
%        end 
        end
        
    end
    
    i = i + 1;
    disp(['#    Step ',num2str(i),'/360']);

end


%% Create Plots

h1 = figure(fig);
set(h1,'Position',[150,150,1500,1500])

sgtitle({[body,'\newline']
    
    ['v_{\infty} = [',num2str(v_inf(1)),',',num2str(v_inf(2)),...
    ',',num2str(v_inf(3)),'] km/s']
    
    ['h_{entry} = ',num2str(hEntry),' km']
    
    [' \Delta \beta = [',num2str(theta_min),':',...
    num2str(theta_delta),':',num2str(theta_max),']']
    });

% sgtitle(planet)

% Surf Plot, Min Flight Path Angle Plot, Subplot 1
% ------------------------------------------------

ax1 = subplot(row,column,plotNr);

% surf plot of ellipsoid
s1 = surf(x,y,z);

% colormap and color data matrix definition
colormap(ax1,flipud(jet))
grid on
% CMap1 = NaN(nrFacets + 1,nrFacets + 1,1);
% CMap2 = NaN(nrFacets + 1,nrFacets + 1,1);

% colormap properties
s1.CDataMapping = 'scaled';
s1.CData = CMap_phiMin;
s1.EdgeColor = 'none';
caxis([-90 0])
axis_limit = 1.5 * radii(1);
colorbar
legend off
title('Min. flight path angle, deg')
xlabel('X, km')
ylabel('Y, km')
zlabel('Z, km')
xlim([-axis_limit,axis_limit])
ylim([-axis_limit,axis_limit])
zlim([-axis_limit,axis_limit])

hold on 
daspect([1 1 1])
arrow3(-2 * axis_limit * v_inf'/norm(v_inf),[0,0,0], 'b-')


% Surf Plot, Max Flight Path Angle Plot, Subplot 2
% ------------------------------------------------
plotNr = plotNr + 1;

ax2 = subplot(row,column,plotNr);

% surf plot of ellipsoid
s2 = surf(x,y,z);

% colormap and color data matrix definition
colormap(ax2,flipud(jet))
grid on
CMap1 = NaN(nrFacets + 1,nrFacets + 1,1);
CMap2 = NaN(nrFacets + 1,nrFacets + 1,1);

% colormap properties
s2.CDataMapping = 'scaled';
s2.CData = CMap_phiMax;
s2.EdgeColor = 'none';
caxis([-90 0])
axis_limit = 1.5 * radii(1);
colorbar
legend off
title('Max. flight path angle, deg')
xlabel('X, km')
ylabel('Y, km')
zlabel('Z, km')
xlim([-axis_limit,axis_limit])
ylim([-axis_limit,axis_limit])
zlim([-axis_limit,axis_limit])

hold on
daspect([1 1 1])
arrow3(-2 * axis_limit * v_inf'/norm(v_inf),[0,0,0], 'b-')


% Surf Plot, Min. Orbital Entry Velocity Plot, Subplot 3
% ------------------------------------------------------
plotNr = plotNr + 1;

ax3 = subplot(row,column,plotNr);

% surf plot of ellipsoid
s3 = surf(x,y,z);

% colormap and color data matrix definition
colormap(ax3,jet)
grid on

% plot formatting
hold on

% colormap properties
s3.CDataMapping = 'scaled';
s3.CData = CMap_vEntryMin;
s3.EdgeColor = 'none';
%caxis([-90 0])
axis_limit = 1.5 * radii(1);
colorbar
legend off
title('Min. orbital velocity @ entry, km/s')
xlabel('X, km')
ylabel('Y, km')
zlabel('Z, km')
xlim([-axis_limit,axis_limit])
ylim([-axis_limit,axis_limit])
zlim([-axis_limit,axis_limit])
hold on 
daspect([1 1 1])
arrow3(-2 * axis_limit * v_inf'/norm(v_inf),[0,0,0], 'b-')

% Surf Plot, Max Orbital Entry Velocity Plot, Subplot 4
% -----------------------------------------------------
plotNr = plotNr + 1;
ax4 = subplot(row,column,plotNr);

% surf plot of ellipsoid
s4 = surf(x,y,z);

% colormap and color data matrix definition
colormap(ax4,jet)
grid on

% plot formatting
hold on

% colormap properties
s4.CDataMapping = 'scaled';
s4.CData = CMap_vEntryMax;
s4.EdgeColor = 'none';
%caxis([-90 0])
axis_limit = 1.5 * radii(1);
colorbar
legend off
title('Max. orbital velocity @ entry, km/s')
xlabel('X, km')
ylabel('Y, km')
zlabel('Z, km')
xlim([-axis_limit,axis_limit])
ylim([-axis_limit,axis_limit])
zlim([-axis_limit,axis_limit])
hold on 
daspect([1 1 1])
arrow3(-2 * axis_limit * v_inf'/norm(v_inf),[0,0,0], 'b-')



% Surf Plot, Min Rel. Entry Velocity Plot, Subplot 5
% --------------------------------------------------
plotNr = plotNr + 1;

ax5 = subplot(row,column,plotNr);

% surf plot of ellipsoid
s5 = surf(x,y,z);

% colormap and color data matrix definition
colormap(ax5,jet)
grid on

% plot formatting
hold on

% colormap properties
s5.CDataMapping = 'scaled';
s5.CData = CMap_vRelEntryMin;
s5.EdgeColor = 'none';
%caxis([-90 0])
axis_limit = 1.5 * radii(1);
colorbar
legend off
title('Min. equal atmospheric-relative entry velocity, km/s')
xlabel('X, km')
ylabel('Y, km')
zlabel('Z, km')
xlim([-axis_limit,axis_limit])
ylim([-axis_limit,axis_limit])
zlim([-axis_limit,axis_limit])
hold on 
daspect([1 1 1])
arrow3(-2 * axis_limit * v_inf'/norm(v_inf),[0,0,0], 'b-')


% Surf Plot, Max Rel. Entry Velocity Plot, Subplot 6
% --------------------------------------------------
plotNr = plotNr + 1;
ax6 = subplot(row,column,plotNr);

% surf plot of ellipsoid
s6 = surf(x,y,z);

% colormap and color data matrix definition
colormap(ax6,jet)
grid on

% plot formatting
hold on

% colormap properties
s6.CDataMapping = 'scaled';
s6.CData = CMap_vRelEntryMax;
s6.EdgeColor = 'none';
%caxis([-90 0])
axis_limit = 1.5 * radii(1);
colorbar
legend off
title('Max. equal atmospheric-relative entry velocity, km/s')
xlabel('X, km')
ylabel('Y, km')
zlabel('Z, km')
xlim([-axis_limit,axis_limit])
ylim([-axis_limit,axis_limit])
zlim([-axis_limit,axis_limit])
hold on 
daspect([1 1 1])
arrow3(-2 * axis_limit * v_inf'/norm(v_inf),[0,0,0], 'b-')


% Surf Plot, Rotational Speed, Subplot 7
% --------------------------------------
plotNr = plotNr + 1;
ax7 = subplot(row,column,plotNr);

% surf plot of ellipsoid
s7 = surf(x,y,z);

% colormap and color data matrix definition
colormap(ax7,jet)
grid on

% colormap properties
s7.CDataMapping = 'scaled';
s7.CData = sign(omega) * CMap_vrot;
s7.EdgeColor = 'none';
if sign(omega) == 1
    caxis([0 ceil(max(vRot)*1e3/50) * 50])
else 
    caxis([floor(min(vRot)*1e3/50) * 50 0])
end
cbh = colorbar;
legend off
title(['Rotational Speed, m/s'])
xlabel('X, km')
ylabel('Y, km')
zlabel('Z, km')
xlim([-axis_limit,axis_limit])
ylim([-axis_limit,axis_limit])
zlim([-axis_limit,axis_limit])
hold on 
daspect([1 1 1])
arrow3(-2 * axis_limit * v_inf'/norm(v_inf),[0,0,0], 'b-')




% Plot, Rotational Speed, Subplot 8
% ---------------------------------
plotNr = plotNr + 1;
 
ax8 = subplot(row,column,plotNr);

s8 = plot(ax8,vRot.*1e3,rad2deg(lat));
title('Rotational Speed, m/s')
xlabel('vRot, m/s')
ylabel('latitude, deg')
if 0 < omega
    caxis([0 ceil(max(vRot)*1e3/50) * 50])
else 
    caxis([floor(min(vRot)*1e3/50) * 50 0])
end
ylim([-100 100])

% rotational period of planet, in d

hours = (2 * cspice_pi / abs(omega))/cspice_spd * 24;
minutes = rem(hours,1) * 60; 
seconds = rem(minutes,1) * 60;

dim = [.75 .225 .2 .2];
 
str = ['rotational period: ', num2str(floor(hours)),' h ',num2str(floor(minutes)),...
    ' min ',num2str(floor(seconds)),'s'];
annotation('textbox',dim,'String',str,'FitBoxToText','on');



% % Plot, inertial interface velocity over incoming vInf, Subplot 5
% % ---------------------------------------------------------------
%  
% ax5 = subplot(row,column,5);
% 
% % plot: inertial vInterface over vInf
% vInf = [1:1:12]'; % km/s
% % Hyperbolic Velocity v_0 at interface altitude, [km/s]
% for i = 1:length(rInterface)
%     [vHyp(i,:),vEsc(i,:)] = hyperbolicVelocity(mu,rInterface(1,i),vInf);
% end
% rInterface = zeros(1,101);
% rInterface(1,i+1) = norm(r_idx);
% 
% s5 = surf(vInf,rInterface,vHyp);
% s5.CDataMapping = 'scaled';
% title(['Hyperbolic velocity @ entry over excess velocity and entry radius'])
% xlabel('incoming excess velocity, km/s')
% ylabel('radius @ entry, km')
% zlabel('hyperbolic velocity @ entry, km/s')

toc

%% Validation plot with rings
fig = fig + 1;

figure(fig)

% planet ellipsoid

    % surface ellipsoid point matrices
    [aX,aY,aZ]=Ellipsoid(radii,u,ceil(nrFacets/4));

    % surf plot of ellipsoid
    s11 = surf(aX,aY,aZ,'FaceColor',[0 0.45 0.75]);
    s11.EdgeColor = 'black';
    % s11.FaceColor = 'none';

for i = 1:3:size(trajStates,1)-2
    
    hold on
    scatter3(trajStates(i,:),trajStates(i+1,:),trajStates(i+2,:),...
        '.','MarkerEdgeColor',[0.3 0.6 0],'MarkerFaceColor',[0.3 0.6 0])
    
end

hold on
states_temp = states_hazard(:,states_hazard(7,:)==1);
scatter3(states_temp(1,:),states_temp(2,:),states_temp(3,:),...
        '.','MarkerEdgeColor','green','MarkerFaceColor','green')

hold on
states_temp = states_hazard(:,states_hazard(7,:)==0);
scatter3(states_temp(1,:),states_temp(2,:),states_temp(3,:),...
        '.','MarkerEdgeColor','red','MarkerFaceColor','red')

% max hazardous state component
axMax = max(max(states_hazard(1,:)),max(states_hazard(2,:)));    

keySet = keys(rings);
for i = 1:length(keySet)
    hold on
    borders = rings(keySet{i});
    
    [bx,by,bz] = circle([0;0;0],borders(1),[0;0;1]);
    [bx1,by1,bz1] = circle([0;0;0],borders(2),[0;0;1]);
    
    fill3([bx;bx1],[by;by1],[bz;bz1],[0.95 0.95 0.95],'EdgeColor','none')

    axMax = max(axMax,borders(2));
    
end



xlim([-axMax axMax])
ylim([-axMax axMax])
zlim([-axMax/2 axMax/2])
title(body)
xlabel('X, km')
ylabel('Y, km')
zlabel('Z, km')

% body inertial reference frame
hold on
daspect([1 1 1])
arrow3([0 0 0],[1 0 0] * 0.5 * axMax, 'b-')
arrow3([0 0 0],[0 1 0] * 0.5 * axMax, 'b-')
arrow3([0 0 0],[0 0 1] * 0.5 * axMax/2, 'b-')

% hyperbolic excess velocity
hold on 
daspect([1 1 1])
arrow3(-1 * axMax * v_inf'/vecnorm(v_inf,2),...
    -0.5 * axMax * v_inf'/vecnorm(v_inf,2), 'k-',0.5,1)


view(-45,-45)

%% Rotational Velocity Plot

fig = fig + 1;

figure(fig)

% atmosphere ellipsoid point matrices
[aX,aY,aZ]=Ellipsoid(radii + hEntry,u,ceil(nrFacets/4));

% ellipsoid point vector
[aa,bb] = size(aX);
jj = 0;
aP = zeros(aa * bb,5);
for kk    = 1:bb 
    for ii = 1:aa
        jj = jj + 1; 
        aP(jj,1) = aX(ii,kk); 
        aP(jj,2) = aY(ii,kk);
        aP(jj,3) = aZ(ii,kk);
        aP(jj,4) = ii;
        aP(jj,5) = kk;
    end
end

% surf plot of atmosphere
scatter3(aP(:,1),aP(:,2),aP(:,3),...
    '.','MarkerEdgeColor','green','MarkerFaceColor','green')

% limit axes
axMax = max(max(radii(1),radii(2)),radii(3));

% rotational velocity vectors
daspect([1 1 1])
for i = 1:40:length(v_rot_entry)
    hold on
    arrow3(state_entry(1:3,i)',...
        state_entry(1:3,i)' + [v_rot_entry(1:3,i) * 1e3]', 'b-')
end


% % all hyperbolic veclocities
% daspect([1 1 1])
% for i = 2:size(vInf,2)
%     hold on
%     arrow3(-1.5 * axMax * [rotate * vInf(:,i)/vecnorm(vInf(:,i),2)]',[0,0,0], 'r-')
% end

% hyperbolic excess velocity
hold on 
daspect([1 1 1])
arrow3(- 1.5 * axMax * v_inf'/vecnorm(v_inf,2),...
    -2 * axMax * v_inf'/vecnorm(v_inf,2), 'y-',0.5,1)


% body inertial reference frame
hold on 
daspect([1 1 1])
arrow3([0 0 0],[1 0 0] * 1.5 * axMax, 'b-')
arrow3([0 0 0],[0 1 0] * 1.5 * axMax, 'b-')
arrow3([0 0 0],[0 0 1] * 1.5 * axMax, 'b-')

% hyperbolic excess velocity
hold on 
daspect([1 1 1])
arrow3(-1.75 * axMax * v_inf'/vecnorm(v_inf,2),...
    -1.25 * axMax * v_inf'/vecnorm(v_inf,2), 'k-',0.5,1)


xlim([-2 * axMax 2 * axMax])
ylim([-2 * axMax 2 * axMax])
zlim([-2 * axMax 2 * axMax])
title(body)
xlabel('X, km')
ylabel('Y, km')
zlabel('Z, km')

view(0,90)
% 
% 
% fig = fig + 1;
% figure(fig)
% plot(v_rot_entry(1,:),v_rot_entry(2,:))
% 
% fig = fig + 1;
% figure(fig)
% plot(v_rot_entry(1,:),v_rot_entry(2,:))
% plot(reclat(2,2:end),reclat(3,2:end))

%% Interplanetary Trajectory Plots

% one planet
if strcmp(body,'Mars')
    % do nothing
else
    fig = fig + 1;
    h2 = figure(fig);
    set(h2,'Position',[130,175,1500,400])
    sgtitle(body)

    % subplot definition
    row = 1;
    column = 3;
    plotNr = 1;

    v_inf_mag = vecnorm([trajectories{:,3} trajectories{:,4} trajectories{:,5}] ...
        - [trajectories{:,9} trajectories{:,10} trajectories{:,11}],2,2);
    string = cspice_etcal(trajectories{:,1}'*cspice_spd);
    ToF = (trajectories{:,2} - trajectories{:,1})/365;
    
    ax9 = subplot(row,column,plotNr);
    scatter(trajectories{:,1}/365,ToF,1,v_inf_mag)
    colormap(ax9,jet)
    colorbar
    ylim([floor(min(ToF)) ceil(max(ToF))+1])
    xlim([floor(min(trajectories{:,1}/365))-1 ceil(max(trajectories{:,1}/365))+1])
    title(['v_{\infty}, km/s'])
    xlabel('Launch Date, yr past 2000')
    ylabel('ToF, yr')
    
    plotNr = plotNr + 1;
    ax10 = subplot(row,column,plotNr);
    scatter(trajectories{:,1}/365,ToF,1,trajectories{:,14})
    colormap(ax10,jet)
    colorbar
    xlim([floor(min(trajectories{:,1}/365))-1 ceil(max(trajectories{:,1}/365))+1])
    ylim([floor(min(ToF)) ceil(max(ToF))+1])
    title(['S/C mass @ arrival, kg'])
    xlabel('Launch Date, yr past 2000')
    ylabel('ToF, yr')
    
    plotNr = plotNr + 1;
    ax11 = subplot(row,column,plotNr);
    scatter(trajectories{:,1}/365,ToF,1,trajectories{:,13})
    colormap(ax11,jet)
    colorbar
    xlim([floor(min(trajectories{:,1}/365))-1 ceil(max(trajectories{:,1}/365))+1])
    ylim([floor(min(ToF)) ceil(max(ToF))+1])
    title(['Deltav, km/s'])
    xlabel('Launch Date, yr past 2000')
    ylabel('ToF, yr')
    
    
end

%% EPSC Abstract Plots

fig = fig + 1;
h3 = figure(fig);
set(h3,'Position',[150,150,1500,1500])

% sgtitle({[planet,'\newline']
%     
%     ['v_{\infty} = [',num2str(v_inf(1)),',',num2str(v_inf(2)),...
%     ',',num2str(v_inf(3)),'] km/s']
%     
%     ['h_{entry} = ',num2str(hEntry),' km']
%     
%     [' \Delta \beta = [',num2str(theta_min),':',...
%     num2str(theta_delta),':',num2str(theta_max),']']
%     });

sgtitle(body)

% Surf Plot, Min Flight Path Angle Plot, Subplot 1
% ------------------------------------------------

ax31 = subplot(2,2,1);

% surf plot of ellipsoid
s31 = surf(x,y,z);

% colormap and color data matrix definition
colormap(ax31,flipud(jet))
grid on
% CMap1 = NaN(nrFacets + 1,nrFacets + 1,1);
% CMap2 = NaN(nrFacets + 1,nrFacets + 1,1);

% colormap properties
s31.CDataMapping = 'scaled';
s31.CData = CMap_phiMin;
s31.EdgeColor = 'none';
caxis([-90 0])
axis_limit = 1.5 * radii(1);
colorbar
legend off
title('Entry flight path angle, deg')
xlabel('X, km')
ylabel('Y, km')
zlabel('Z, km')
xlim([-axis_limit,axis_limit])
ylim([-axis_limit,axis_limit])
zlim([-axis_limit,axis_limit])

hold on 
daspect([1 1 1])
arrow3(-2 * axis_limit * v_inf'/norm(v_inf),[0,0,0], 'b-')


% Surf Plot, Orbital Entry Velocity Plot, Subplot 2
% -------------------------------------------------

ax32 = subplot(2,2,2);

% surf plot of ellipsoid
s32 = surf(x,y,z);

% colormap and color data matrix definition
colormap(ax32,jet)
grid on

% plot formatting
hold on

% colormap properties
s32.CDataMapping = 'scaled';
s32.CData = CMap_vEntryMin;
s32.EdgeColor = 'none';
%caxis([-90 0])
axis_limit = 1.5 * radii(1);
colorbar
legend off
title('Orbital velocity @ entry, km/s')
xlabel('X, km')
ylabel('Y, km')
zlabel('Z, km')
xlim([-axis_limit,axis_limit])
ylim([-axis_limit,axis_limit])
zlim([-axis_limit,axis_limit])
hold on 
daspect([1 1 1])
arrow3(-2 * axis_limit * v_inf'/norm(v_inf),[0,0,0], 'b-')

% Surf Plot, Rel. Entry Velocity Plot, Subplot 3
% ----------------------------------------------

ax33 = subplot(2,2,3);

% surf plot of ellipsoid
s33 = surf(x,y,z);

% colormap and color data matrix definition
colormap(ax33,jet)
grid on

% plot formatting
hold on

% colormap properties
s33.CDataMapping = 'scaled';
s33.CData = CMap_vRelEntryMin;
s33.EdgeColor = 'none';
%caxis([-90 0])
axis_limit = 1.5 * radii(1);
colorbar
legend off
title('Equal atmospheric-relative entry velocity, km/s')
xlabel('X, km')
ylabel('Y, km')
zlabel('Z, km')
xlim([-axis_limit,axis_limit])
ylim([-axis_limit,axis_limit])
zlim([-axis_limit,axis_limit])
hold on
daspect([1 1 1])
arrow3(-2 * axis_limit * v_inf'/norm(v_inf),[0,0,0], 'b-')


% Surf Plot, Rotational Speed, Subplot 4
% --------------------------------------

ax34 = subplot(2,2,4);

% surf plot of ellipsoid
s34 = surf(x,y,z);

% colormap and color data matrix definition
colormap(ax34,jet)
grid on

% colormap properties
s34.CDataMapping = 'scaled';
s34.CData = sign(omega) * CMap_vrot;
s34.EdgeColor = 'none';
if sign(omega) == 1
    caxis([0 ceil(max(vRot)*1e3/50) * 50])
else 
    caxis([floor(min(vRot)*1e3/50) * 50 0])
end
cbh = colorbar;
legend off
title(['Rotational Speed, m/s'])
xlabel('X, km')
ylabel('Y, km')
zlabel('Z, km')
xlim([-axis_limit,axis_limit])
ylim([-axis_limit,axis_limit])
zlim([-axis_limit,axis_limit])
hold on 
daspect([1 1 1])
arrow3(-2 * axis_limit * v_inf'/norm(v_inf),[0,0,0], 'b-')




% Plot, Rotational Speed
% ----------------------
fig = fig + 1;

figure(fig)

plot(vRot.*1e3,rad2deg(lat));
title('Rotational Speed, m/s')
xlabel('vRot, m/s')
ylabel('latitude, deg')
if 0 < omega
    caxis([0 ceil(max(vRot)*1e3/50) * 50])
else 
    caxis([floor(min(vRot)*1e3/50) * 50 0])
end
ylim([-100 100])

% rotational period of planet, in d

hours = (2 * cspice_pi / abs(omega))/cspice_spd * 24;
minutes = rem(hours,1) * 60; 
seconds = rem(minutes,1) * 60;

dim = [.75 .225 .2 .2];
%  
% str = ['rotational period: ', num2str(floor(hours)),' h ',num2str(floor(minutes)),...
%     ' min ',num2str(floor(seconds)),'s'];
% annotation('textbox',dim,'String',str,'FitBoxToText','on');


% %% Create Launch Opportunity Plot for all Planets
% 
% % all planets
% 
% fig = fig + 1;
% h3 = figure(fig);
% set(h3,'Position',[50,50,1500,1000])
% 
% % subplot definition
% keySet = keys(M);
% row = length(keySet);
% column = 3;
% plotNr = 1;
% 
% titletext = 'Launch Opportunities for ';
% 
% for i = 2:row
%     
%     
%     planet = keySet{i};
%     disp(planet)
%     
%     if i == row
%         titletext = [titletext,planet];
%     elseif i == row - 1
%         titletext = [titletext,planet,' and '];
%     else
%         titletext = [titletext,planet,', '];
%     end
%     
%     [trajectories,~,~,~,~,~,~,~] = planetData(planet);
%     trajectory = trajectories{1,:};
%     
%     v_inf_mag = vecnorm([trajectories{:,3} trajectories{:,4} trajectories{:,5}] ...
%     - [trajectories{:,9} trajectories{:,10} trajectories{:,11}],2,2);
%     string = cspice_etcal(trajectories{:,1}'*cspice_spd);
%     ToF = (trajectories{:,2} - trajectories{:,1})/365;
%     
% %     switch i
% %         case 2
% %             h1 = text(-15, 0.35,planet);
% %             set(h1, 'rotation', 90)
% %         case 3
% %             h1 = text(-10, 0.25,planet);
% %             set(h1, 'rotation', 90)
% %         case 4
% %             h1 = text(-15, 0.15,planet);
% %             set(h1, 'rotation', 90)
% %     end
%     
%     for k = 1:column
%         ax = subplot(row,column,plotNr);
%         switch k
%             case 1
%                 scatter(trajectories{:,1}/365,ToF,1,v_inf_mag)
%                 title(['v_{\infty}, km/s'])
%             case 2
%                 scatter(trajectories{:,1}/365,ToF,1,trajectories{:,14})
%                 title(['S/C mass @ arrival, kg'])
%             case 3
%                 scatter(trajectories{:,1}/365,ToF,1,trajectories{:,13})
%                 title(['Deltav, km/s'])
%         end
%         colormap(ax,jet)
%         colorbar
%         ylim([floor(min(ToF)) ceil(max(ToF))+1])
%         xlim([floor(min(trajectories{:,1}/365))-1 ceil(max(trajectories{:,1}/365))+1])
%         xlabel('Launch Date, yr past 2000')
%         ylabel('ToF, yr')
% 
%         plotNr = plotNr + 1;
%     end
% 
% end
%     
% sgtitle(titletext);



%% clear spice kernels
cspice_kclear
