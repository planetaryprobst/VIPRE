function [data,out] = IPED_pythonInterface(body,ID,vInf,epoch,hEntry,swPath,prints)

%--------------------------------------------------------------------------
%   Interface script called by Python setup function. Computes entry 
%   conditions for all safe entry opportunities for one trajectory. 
%
%   Takes as input the name of the 'body', the trajectory 'ID' assigned by
%   the Python script, the hyperbolic excess velocity 'vInf' of the
%   trajectory, the time of arrival 'epoch' and the entry altitude of the
%   probe 'hEntry'.
%
%   Adapts the entry variables to Matlab format and passes them on to the
%   IPED_python main script.
%
%   IPED is part of the software package VIPRE, consisting of VAPRE and 
%   IPED,  a software to Visualize the Impact of the PRobe Entry location 
%   on the spacecraft and mission design.
%
%--------------------------------------------------------------------------
%   Form:
%   [data,out] = IPED_pythonInterface(body,ID,vInf,epoch,hEntry,swPath,prints)
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
%   swPath       str     -          optional: path to software directory 
%                                   If empty, local directory is assumed
%   prints      boolean             True/False to (de)activate screen
%                                   output
%
%   ------
%   Output
%   ------
%   data        bool    -           success boolean
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 25.05.2020 |  A. Probst    | First revision
%*************************************************************************%

% When swPath empty, local directory is assumed
if isempty(swPath)
    swPath = '../IPED-python/src/data/';
end

% Adapting variable format
vInf = vInf';

%% NAIF ID associated to body

if ischar(body)
    
    [ ~, found] = cspice_bodn2c(body);

    if found==0
        disp('The input planet is not a valid planet!')
        return
    end
    
else 
    
    [body, found] = cspice_bodc2n(body);
    
end

if found==0
    disp('The input planet is not a valid planet!')
    return
end

%% Call IPED 

% ID = [planet,num2str(hEntry),'-vID',num2str(vID)];

[data,out] = IPED_python(body,ID,vInf,epoch,hEntry,swPath,prints);
    

end
