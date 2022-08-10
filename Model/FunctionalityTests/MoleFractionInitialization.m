%% Testing mole fraction initialization
clear all;close all;clc;

%% Set Parameters
FLAG.KnownZ = 1;

load_filename = 'WileyData.mat';

AN.L         = 100e-6;
AN.A_c       = 0.1;
AN.eps_ed    = 0.320700992697794;
% AN.eps_ed    = 0.471;
AN.C_Li_max  = 28.4851292093828;
% AN.C_Li_max  = 31.542;
AN.EeqHandle = @E_eqGraphite;

CA.L         = 100e-6;
CA.A_c       = 0.1;
CA.eps_ed    = 0.442994582876546;
CA.C_Li_max  = 40.544149086751;
% CA.eps_ed    = 0.297;
% CA.C_Li_max  = 49.668;
CA.EeqHandle = @E_eqNMC;

x_lim = [0,1];
y_lim = [0,1];
% x_lim = [0 , 1];
% y_lim = [0.3, 0.7];


SIM.options_curvefit = optimoptions('lsqcurvefit');
SIM.options_curvefit.Display = 'off';
% SIM.options_curvefit.OptimalityTolerance = 1e-8; % Original is 1e-6
SIM.options_curvefit.FunctionTolerance   = 1e-8; % Original is 1e-6
SIM.options_curvefit.StepTolerance       = 1e-8; % Original is 1e-6
% SIM.options_curvefit.MaxIterations = 500; % Original is 400

CONS.F = 96485338.3; % [C kmol^-1], Faraday's Constant 

%% Pointers
i = 1;
P.y_int = i; i = i+1;
P.x_min = i; i = i+1;
P.x_max = i; i = i+1;
if ~FLAG.KnownZ
    P.z = i; i = i+1;
    num_par = P.z;
else
    num_par = P.x_max;
end

%% Load SOC vs OCP data
load(load_filename)

%% Geometry Calc (Capacity Ratio: z)
AN.Vol  = AN.A_c * AN.L;
AN.V_ed = AN.eps_ed * AN.Vol;
AN.Cap  = AN.C_Li_max * CONS.F *(1/3600)*AN.V_ed;

CA.Vol  = CA.A_c * CA.L;
CA.V_ed = CA.eps_ed * CA.Vol;
CA.Cap  = CA.C_Li_max * CONS.F *(1/3600)*CA.V_ed;

SIM.z = (CA.C_Li_max*CA.V_ed)/(AN.C_Li_max*AN.V_ed);

%% Initialize
p0 = zeros(num_par,1);
lb = zeros(num_par,1);
ub = zeros(num_par,1);

p0(P.y_int) = 1;
lb(P.y_int) = 0;
ub(P.y_int) = 2;

p0(P.x_min) = x_lim(1);
lb(P.x_min) = 0;
ub(P.x_min) = 0.5;

p0(P.x_max) = x_lim(2);
lb(P.x_max) = 0.5;
ub(P.x_max) = 1;

if ~FLAG.KnownZ
    p0(P.z) = 1;
    lb(P.z) = 0.5;
    ub(P.z) = 2;
end

%% Solution
params = lsqcurvefit(@(params,xdata)myFun(params,xdata,AN,CA,P,SIM,FLAG),p0,xdata,ydata,lb,ub,SIM.options_curvefit);

%% Post-Calcs
x_lim_soln = [params(P.x_min),params(P.x_max)];
y_int_soln = params(P.y_int);
if FLAG.KnownZ
    z_soln = SIM.z;
else
    z_soln = params(P.z);
end

% Same size as SOC
x_soln   = (xdata/100)*(params(P.x_max) - params(P.x_min)) + params(P.x_min);
y_soln   = YfromX(x_soln,z_soln,y_int_soln);
Eeq_an   = AN.EeqHandle(x_soln);
Eeq_ca   = CA.EeqHandle(y_soln);
OCV_soln = Eeq_ca - Eeq_an;



%% Plot
% x,y plot
    figure
    plot(x_soln , y_soln, 'Linewidth' , 2 )
    xlim([0 , 1])
    ylim([0 , 1])
    %Add stoich limits

% SOC vs x,y
    figure
    title('X_{surf} vs SOC')
    xlabel('SOC (%)')

    yyaxis left
    plot(xdata , x_soln , 'Linewidth' , 2 )
    ylabel('Anode')

    yyaxis right
    plot(xdata , y_soln, 'Linewidth' , 2 )
    ylabel('Cathode')

% SOC vs V and OCP
    figure
    hold on
    plot(xdata,ydata,'LineWidth',2,'DisplayName','Experimental Data')
    plot(xdata,OCV_soln,'ok','DisplayName','Match Data')
    xlabel('SOC (%)')
    ylabel('OCV (V)')
    lgn = legend;
    lgn.Location = 'southeast';
    xlim([0 , 100])
    ylim([3.2 , 4.4])


%% Functions
%% Fit Function
function OCV = myFun(params,xdata,AN,CA,P,SIM,FLAG)
% Inputs
    % params  : The parameters that I am trying to fit
    % xdata   : dep variable used to evaluate the parameter fit. This is SOC
    % Anode   : The OCV function for the anode
    % Cathode : The OCV function for the cathode
    % P       : pointers used for params

% Outputs
    % OCV     : The open circuit potential from the fitted parameters

%% Convert SOC to mole fraction (x and y)
x = (xdata/100)*(params(P.x_max) - params(P.x_min)) + params(P.x_min);

if ~FLAG.KnownZ
    y = YfromX(x , params(P.z) , params(P.y_int));
else
    y = YfromX(x , SIM.z       , params(P.y_int));
end

%% Solve for voltage
Eeq_an = AN.EeqHandle(x);
Eeq_ca = CA.EeqHandle(y);
OCV    = Eeq_ca - Eeq_an;

end

%% Get y from x
function y_out = YfromX(x,z,y_intcep)
    y_out = -1/z*x + y_intcep;
end