function [data,out] = IPED_python(body,ID,vInf,epoch,hEntry,swPath,prints)
%--------------------------------------------------------------------------
%   Main script for computations called by IPED_pythonInterface. Computes
%   the entry conditions for all safe entry opportunities for one
%   trajectory  to the input 'body' with the identifier 'ID' and the
%   hyperbolic velocity 'vInf', the arrival time 'epoch', at the entry
%   altitude 'hEntry'.
%
%   The computations are equivalent to the computations done with the
%   function IPED.m. 
%
%   Saves results in a text file.
%
%--------------------------------------------------------------------------
%   Form:
%   [data,out] = IPED_python(body,ID,vInf,epoch,hEntry,swPath,prints)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   planet      str      -          SPICE code or string for planet
%   ID          str      -          trajectory ID, assigned by python
%                                   script
%   vInf        list/3  km/s        Python object - hyperbolic entry 
%                                   veclocity per trajectory in IAU_ body  
%                                   inertial reference frame
%   epoch       (1,1)    days       arrival time, days after JD2000 
%   hEntry      (1,1)    km         probe entry altitude 
%   swPath    str     -             path to software directory 
%   prints      boolean             True/False to (de)activate screen
%                                   output
%
%   ------
%   Output
%   ------
%   data
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 25.05.2020 |  A. Probst    | First revision
%*************************************************************************%
   

%% Reference Frame Definitions

    from = 'ECLIPJ2000';
    to = ['IAU_',body];
    

%% Loop Data

    % number of computed states between initial and periapsis
    n = 1e3;
    
    % clock angle B-plane interval
    theta.min = -cspice_pi;
    theta.max = cspice_pi;
    
    switch body
        case 'Saturn'
            % clock angle B-plane interval step
            theta.delta = cspice_pi/89;
            % B-vector length interval step
            Babs_delta = 1750; % km (1500)
        otherwise
            theta.delta = cspice_pi/(2 * 89 + 1);
            % B-vector length interval step
            Babs_delta = 750; % km
    end

%% Planet Definition

% planet data retrieval
[mu,radii,~,rings] = planetData(body);


%% Input Data Definition

% ephemeris time, in s past J2000 epoch
et = epoch .* cspice_spd;

% rotation matrix for vInf vector transformation             
rotate = cspice_pxform( from, to, et );

% vInf transformation
v_inf = rotate * vInf;

% number of patches of ellipsoid (graphics)
nrFacets = 150;

%% Planet's Rotation omega
% in rad/s
omega = planetRotation(body);

%% Variable Definition

% B-vector length interval step
Babs_min = 0; % 10; 
Babs_max = 8 * (radii(1)+radii(2))/2;

% % Color Maps
% CMap_generic = NaN(nrFacets + 1,nrFacets + 1,1);
% % defined based on CMap_generic
% CMap.phiMin = CMap_generic;
% CMap.phiMax = CMap_generic;
% CMap.vEntryMin = CMap_generic;
% CMap.vEntryMax = CMap_generic;
% CMap.vRelEntryMin = CMap_generic;
% CMap.vRelEntryMax = CMap_generic;
% CMap.vrot = CMap_generic;


% %% Create the Planet's Sphere 
% % as a preparation for the plots
% 
% [x,y,z]=plot_planetSphere(radii,nrFacets);
% % % atmosphere point matrices
% % [aX,aY,aZ]=ellipsoid3(radii + hEntry,u,nrFacets/4);
% 
% % ellipsoid point vector
% [a,b] = size(x);
% j = 0;
% P = zeros(a * b,5);
% 
% for k    = 1:b 
%     for i = 1:a
%         j = j + 1; 
%         P(j,1) = x(i,k); 
%         P(j,2) = y(i,k);
%         P(j,3) = z(i,k);
%         P(j,4) = i;
%         P(j,5) = k;
%     end
% end

%% Prepare file to safe data

if hEntry < 10000
    prefix = '0';
    if hEntry < 1000
        prefix = [prefix,'0'];
    end
else
    prefix = [];
end
   

% Create and open file
savePath = [swPath,body,'/',prefix,num2str(hEntry)];

if ~exist(savePath, 'dir')
   mkdir(savePath)
end

filename = [savePath,'/',ID,'.txt'];
fileID  = fopen(filename,'w');

% header
fprintf(fileID, ['Bvec-theta(rad) \t Bvec-abs(km) \t entryTrajec(boolean) ', ...
    '\t stateEqu-rX(km) \t stateEqu-rY(km) \t stateEqu-rZ(km) ', ...
    '\t stateEqu-vX(km/s) \t stateEqu-vY(km/s) \t stateEqu-vZ(km/s)', ...
    '\t safe(boolean) \t entryState-rX(km) \t entryState-rY(km) \t entryState-rZ(km) ', ...
    '\t entryState-vX(km/s) \t entryState-vY(km/s) \t entryState-vZ(km/s) ', ...
    '\t lon_entry(rad) \t lat_entry(rad) ', ...
    '\t vRot-x(km/s) \t vRot-y(km/s) \t vRot-z(km/s) \t FPA(rad) ',...
    '\t vRel_entry-x(km/s) \t vRel_entry-y(km/s) \t vRel_entry-z(km/s)\n']);


%% Varying the Incoming Trajectory

i = 1;

reclat = [];

% % number of total iterations
% iter = ((theta.max - theta.min)/theta.delta + 1) * ...
%     ((Babs_max - Babs_min)/Babs_delta + 1);
% states_hazard = [];
% nrRngCrssng = [];

%% LOOP

if prints
    % status print: number of current iteration
    disp(['#    Step ',num2str(i),'/360']);
end

% loop over B_theta
for Btheta = theta.min:theta.delta:theta.max

    % Creating B-plane
    [~,~,~,B_uv,h_uv,~,a] = VinfThe2B_plane(v_inf,Btheta,mu);

    % loop over length of B_vec
    for Babs = Babs_min:Babs_delta:Babs_max
        
        % resetting of print variables
            state_hazard = NaN(6,1);
            safe = NaN;
            lon_idx = NaN;
            lat_idx = NaN ;
            v_rot = NaN(3,1);
            phi_idx = NaN; 
            v_entry_idx = NaN(3,1);

        if Babs == 0
            [radius, lon0, lat0] = cspice_reclat(-v_inf);
            reclat = [reclat,[radius;rad2deg(lon0);rad2deg(lat0)]];
            Babs_min = Babs_min + Babs_delta;
            continue
        end

        % Calculating Orbit Parameters:

        % B-vector, in km
        B_vec = Babs * B_uv;
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

        % Computing the orbit until periapsis
        [states,elts] = computingOrbStates(elts0,n,true,et);

        % Locating point of closest approach, or, if given, closest point
        % to entry height
        [altitude,~,idx,flag] = ...
            closestApproach(states(1:3,:),radii,hEntry);
        
        % If trajectory is an entry trajectory    
        if flag
            % Check if it is a hazardous trajectory or not
            [state_hazard,safe] = hazardCrssng(states,elts,et,radii,rings);
            
            %states_hazard = [states_hazard,[state_hazard;safe]];
            if safe
                % Calculating the Flight Path Angle at entry / closest approach
                % determining entry point of trajectory
                diff = altitude(1,idx) - hEntry;
                if diff > 0.1
                    % call function to determine entry point
                    [~,~,~,states,idx,~] = ...
                                  entryPoint(elts,idx,radii,hEntry,10,et,diff);

%                     % collecting trajectory state at entry
%                     trajStates = [trajStates,states(1:3,idx)];

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

%                 % identify point on sphere
%                 k = dsearchn(P(:,1:3),states(1:3,idx)');
% 
%                 CMap_temp = NaN(nrFacets + 1,nrFacets + 1,1);
% 
%                 % allocate values to identified point(s)
%                 % flight path angle
%                 CMap_temp(k) = rad2deg(phi_idx);
%                 CMap.phiMin = min(CMap_temp,CMap.phiMin);
%                 CMap.phiMax = max(CMap_temp,CMap.phiMax);
%                 % orbital velicity @ entry
%                 CMap_temp(k) = vecnorm(states(4:6,idx),2);
%                 CMap.vEntryMin = min(CMap_temp,CMap.vEntryMin);
%                 CMap.vEntryMax = max(CMap_temp,CMap.vEntryMax);
%                 % relative entry velocity 
%                 CMap_temp(k) = vecnorm(v_entry_idx,2);
%                 CMap.vRelEntryMin = min(CMap_temp,CMap.vRelEntryMin);
%                 CMap.vRelEntryMax = max(CMap_temp,CMap.vRelEntryMax);
%                 % rotational velocity
%                 CMap_temp(k) = vecnorm(v_rot,2) * 1e3;
%                 CMap.vrot = max(CMap_temp,CMap.vrot);            
                
            else
                % determining entry point of trajectory
                diff = altitude(1,idx) - hEntry;
                if diff > 0.1
                    % call function to determine entry point
                    [~,~,~,states,idx,~] = ...
                                  entryPoint(elts,idx,radii,hEntry,10,et,diff);
                end
            end
                   
            % data vector
            data = [Btheta Babs flag state_hazard' safe states(:,idx)' ...
                lon_idx lat_idx v_rot' phi_idx v_entry_idx']; 

            fprintf(fileID, '%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', data); 

        else
            % FPA for closest approach (if no entry) is periapsis, 
            % therefore FPA ~ 0
            continue
            
%        end 
        end

    end
    
    % print out on screen as running variable, to indicate status of
    % computation
    if prints
        i = i + 1;
        disp(['#    Step ',num2str(i),'/360']);
    end
end

if prints
    disp('#*************************************************************************#')
    disp('#')
end

%% Closing File

fclose(fileID);

%% Creating Plots
out = false;
% out = createPlots(planet,v_inf,hEntry,theta,CMap,radii,omega,...
%     trajStates,states_hazard,rings,state_entry,v_rot_entry,nrFacets);

%% Return success boolean
data = true;
end
