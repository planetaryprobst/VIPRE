function out = createPlots(planet,vInf,hEntry,theta,CMap,radii,omega,...
    trajStates,states_hazard,rings,state_entry,v_rot_entry,nrFacets)
%--------------------------------------------------------------------------
%   Summarizes the code for creating plots with IPED. Same plots as in
%   IPED_matlab. This function was initially called by IPED_python, however
%   hasn't been used since after the initial testing phase.
%
%   Returns a boolean and opens plots for graphics.
%--------------------------------------------------------------------------
%   Form:
%   out = createPlots(planet,vInf,hEntry,theta,CMap,radii,omega,...
%          trajStates,states_hazard,rings,state_entry,v_rot_entry,nrFacets)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%
%   ------
%   Output
%   ------  
% 
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 25.05.2020 |  A. Probst    | First revision
%*************************************************************************%

%% Atmospheric Speed depending on Latitude FOR GRAPHICS

% rotational speed for latitude range [-90 deg:90 deg] in m/s
lat = [-cspice_pi/2:cspice_pi/(2 * 89):cspice_pi/2]';

[vRot,~] = rotationalSpeed(lat,omega,radii+hEntry);


%% Graphics Data

    % figure counter
    fig = 1;

    % subplot definition
    row = 4;
    column = 2;
    plotNr = 1;
    
%% Create the Planet's Sphere 

% as a preparation for the plots
[x,y,z] = plot_planetSphere(radii,nrFacets);
    
%% Create Plots

h1 = figure(fig);
set(h1,'Position',[150,150,1500,1500])

sgtitle({[planet,'\newline']

    ['v_{\infty} = [',num2str(vInf(1)),',',num2str(vInf(2)),...
    ',',num2str(vInf(3)),'] km/s']

    ['h_{entry} = ',num2str(hEntry),' km']

    [' \Delta \beta = [',num2str(theta.min),':',...
    num2str(theta.delta),':',num2str(theta.max),']']
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
s1.CData = CMap.phiMin;
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
arrow3(-2 * axis_limit * vInf'/norm(vInf),[0,0,0], 'b-')


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
s2.CData = CMap.phiMax;
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
arrow3(-2 * axis_limit * vInf'/norm(vInf),[0,0,0], 'b-')


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
s3.CData = CMap.vEntryMin;
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
arrow3(-2 * axis_limit * vInf'/norm(vInf),[0,0,0], 'b-')

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
s4.CData = CMap.vEntryMax;
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
arrow3(-2 * axis_limit * vInf'/norm(vInf),[0,0,0], 'b-')



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
s5.CData = CMap.vRelEntryMin;
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
arrow3(-2 * axis_limit * vInf'/norm(vInf),[0,0,0], 'b-')


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
s6.CData = CMap.vRelEntryMax;
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
arrow3(-2 * axis_limit * vInf'/norm(vInf),[0,0,0], 'b-')


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
s7.CData = sign(omega) * CMap.vrot;
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
arrow3(-2 * axis_limit * vInf'/norm(vInf),[0,0,0], 'b-')




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


%% Validation plot with rings
fig = fig + 1;

figure(fig)

% planet ellipsoid
[aX,aY,aZ] = plot_planetSphere(radii,ceil(nrFacets/4));

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
title(planet)
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
arrow3(-1 * axMax * vInf'/vecnorm(vInf,2),...
    -0.5 * axMax * vInf'/vecnorm(vInf,2), 'k-',0.5,1)


view(-45,-45)

%% Rotational Velocity Plot

fig = fig + 1;

figure(fig)

% atmosphere ellipsoid point matrices
[aX,aY,aZ] = plot_planetSphere(radii + hEntry,ceil(nrFacets/4));

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
arrow3(- 1.5 * axMax * vInf'/vecnorm(vInf,2),...
    -2 * axMax * vInf'/vecnorm(vInf,2), 'y-',0.5,1)


% body inertial reference frame
hold on 
daspect([1 1 1])
arrow3([0 0 0],[1 0 0] * 1.5 * axMax, 'b-')
arrow3([0 0 0],[0 1 0] * 1.5 * axMax, 'b-')
arrow3([0 0 0],[0 0 1] * 1.5 * axMax, 'b-')

% hyperbolic excess velocity
hold on
daspect([1 1 1])
arrow3(-1.75 * axMax * vInf'/vecnorm(vInf,2),...
    -1.25 * axMax * vInf'/vecnorm(vInf,2), 'k-',0.5,1)


xlim([-2 * axMax 2 * axMax])
ylim([-2 * axMax 2 * axMax])
zlim([-2 * axMax 2 * axMax])
title(planet)
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

% %% Interplanetary Trajectory Plots
% 
% % one planet
% if strcmp(planet,'Mars')
%     % do nothing
% else
%     fig = fig + 1;
%     h2 = figure(fig);
%     set(h2,'Position',[130,175,1500,400])
%     sgtitle(planet)
% 
%     % subplot definition
%     row = 1;
%     column = 3;
%     plotNr = 1;
% 
%     v_inf_mag = vecnorm([trajectories{:,3} trajectories{:,4} trajectories{:,5}] ...
%         - [trajectories{:,9} trajectories{:,10} trajectories{:,11}],2,2);
%     string = cspice_etcal(trajectories{:,1}'*cspice_spd);
%     ToF = (trajectories{:,2} - trajectories{:,1})/365;
% 
%     ax9 = subplot(row,column,plotNr);
%     scatter(trajectories{:,1}/365,ToF,1,v_inf_mag)
%     colormap(ax9,jet)
%     colorbar
%     ylim([floor(min(ToF)) ceil(max(ToF))+1])
%     xlim([floor(min(trajectories{:,1}/365))-1 ceil(max(trajectories{:,1}/365))+1])
%     title(['v_{\infty}, km/s'])
%     xlabel('Launch Date, yr past 2000')
%     ylabel('ToF, yr')
% 
%     plotNr = plotNr + 1;
%     ax10 = subplot(row,column,plotNr);
%     scatter(trajectories{:,1}/365,ToF,1,trajectories{:,14})
%     colormap(ax10,jet)
%     colorbar
%     xlim([floor(min(trajectories{:,1}/365))-1 ceil(max(trajectories{:,1}/365))+1])
%     ylim([floor(min(ToF)) ceil(max(ToF))+1])
%     title(['S/C mass @ arrival, kg'])
%     xlabel('Launch Date, yr past 2000')
%     ylabel('ToF, yr')
% 
%     plotNr = plotNr + 1;
%     ax11 = subplot(row,column,plotNr);
%     scatter(trajectories{:,1}/365,ToF,1,trajectories{:,13})
%     colormap(ax11,jet)
%     colorbar
%     xlim([floor(min(trajectories{:,1}/365))-1 ceil(max(trajectories{:,1}/365))+1])
%     ylim([floor(min(ToF)) ceil(max(ToF))+1])
%     title(['Deltav, km/s'])
%     xlabel('Launch Date, yr past 2000')
%     ylabel('ToF, yr')
% 
% 
% end

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
%     trajecories = trajectories{1,:};
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

   
out = true;
end
