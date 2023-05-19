%% Compare DT and CT Nyquist (EIS)
clear all; close all; clc;

%% Parameters
exp_min  = -7; % -3
exp_max  =  2; %  5

PLOT.Nyquist    = 1;
PLOT.Bode_mag   = 0;
PLOT.Bode_phase = 0;


%% Load SS
load('ROM_SS_SOC50_Ts1.mat')


%% Separate systems
sys = sys_HK{1, end}(1, 1); % end is for the combined ROM, (1,1) is for the cell voltage output
sys_DT = sys;
sys_CT = d2c(sys_DT);


%% Desired Frequency
exp_diff = exp_max - exp_min;
des_freq  = logspace(exp_min ,exp_max ,exp_diff*10+1);


%% Get Nyquist
[re_DT,im_DT,wout_DT] = nyquist(sys_DT,des_freq);
[re_CT,im_CT,wout_CT] = nyquist(sys_CT,des_freq);

re_DT = reshape(re_DT,[],1);
im_DT = reshape(im_DT,[],1);
re_CT = reshape(re_CT,[],1);
im_CT = reshape(im_CT,[],1);

% Get the negative of imaginary impedance
im_DT = - im_DT;
im_DT = - im_DT;

% Sort Results (Sort By Re)
    [~ , newI] = sort(re_DT);%,'descend');
    re_DT_new = re_DT(newI);
    im_DT_new = im_DT(newI);
    
    [~ , newI] = sort(re_CT);%,'descend');
    re_CT_new = re_CT(newI);
    im_CT_new = im_CT(newI);


%% Bode
[mag_DT,phase_DT,wout_DT] = bode(sys_DT,des_freq);
[mag_CT,phase_CT,wout_CT] = bode(sys_CT,des_freq);

mag_DT = reshape(mag_DT,[],1);
phase_DT = reshape(phase_DT,[],1);
mag_CT = reshape(mag_CT,[],1);
phase_CT = reshape(phase_CT,[],1);


%% Plot Results
if PLOT.Nyquist
    % Nyquist
    figure
    hold on
    plot(re_CT_new, im_CT_new ,'ro', 'LineWidth',2,'DisplayName','CT')
    plot(re_DT_new, im_DT_new ,'k' , 'LineWidth',2,'DisplayName','DT')
    lgn = legend;
    axis equal
    % axis square
end

if PLOT.Bode_mag
    % Magnitude
    figure
    semilogx(wout_CT, mag_CT ,'ro', 'LineWidth',2,'DisplayName','CT')
    hold on
    semilogx(wout_DT, mag_DT ,'k' , 'LineWidth',2,'DisplayName','DT')
    lgn = legend;
end

if PLOT.Bode_phase
    % Phase
    figure
    semilogx(wout_CT, phase_CT ,'ro', 'LineWidth',2,'DisplayName','CT')
    hold on
    semilogx(wout_DT, phase_DT ,'k' , 'LineWidth',2,'DisplayName','DT')
    lgn = legend;
end