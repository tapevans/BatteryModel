%% Is their a relationship between the singular value and P_infty
% This is using Input Q step and sine functions
clear all; close all; clc;
%%

FLAG.PLOT = 1;
%%
% svd_n_data: normalized singular value for each simulation
% cpct_data:  the diagonal of the error covariance matrix


i = 0;
i = i+1;
svd_n_data(i,:) = [0.992810668	0.145976983	0.003088145	0.000304958		0.056557093	0.999093262	0.028251629	0.003992017		0.003670038	0.000303175	0.988813961	0.18515142		0.027986863	0.003103248	0.02802977	0.999494132];
cpct_data(i,:)  = [4.894000E-08	9.874234E-09	1.126464E-07	9.753653E-09		6.649257E-08	9.944203E-09	1.329934E-07	9.973881E-09		4.880596E-08	9.800842E-09	9.249775E-08	7.331697E-09		6.634788E-08	9.942359E-09	1.286243E-07	9.369165E-09];

i = i+1;
svd_n_data(i,:) = [0.992810668	0.145976983	0.003088145	0.000304958		0.056557093	0.999093262	0.028251629	0.003992017		0.003670038	0.000303175	0.988813961	0.185151420		0.027986863	0.003103248	0.028029770	0.999494132];
cpct_data(i,:)  = [4.89400E-05	9.87423E-06	1.12646E-04	9.75365E-06		6.64926E-05	9.94420E-06	1.32993E-04	9.97388E-06		4.88060E-05	9.80084E-06	9.24978E-05	7.33170E-06		6.63479E-05	9.94236E-06	1.28624E-04	9.36916E-06];

i = i+1;
svd_n_data(i,:) = [0.992810668	0.145976983	0.003088145	0.000304958		0.056557093	0.999093262	0.028251629	0.003992017		0.003670038	0.000303175	0.988813961	0.185151420		0.027986863	0.003103248	0.028029770	0.999494132];
cpct_data(i,:)  = [4.89400E-02	9.87423E-03	1.12646E-01	9.75365E-03		6.64926E-02	9.94420E-03	1.32993E-01	9.97388E-03		4.88060E-02	9.80084E-03	9.24978E-02	7.33170E-03		6.63479E-02	9.94236E-03	1.28624E-01	9.36916E-03];

i = i+1;
svd_n_data(i,:) = [0.956187115	0.342940611	0.042276736	0.023805843		0.184266903	0.999025630	0.089446021	0.077959659		0.071745162	0.019105460	0.767990749	0.732722618		0.134492433	0.043658532	0.137952717	0.988211155];
cpct_data(i,:)  = [2.46205E-07	4.54284E-08	5.64599E-07	4.86244E-08		3.31567E-07	4.58436E-08	6.63973E-07	4.97204E-08		2.45254E-07	4.50614E-08	4.70565E-07	3.65787E-08		3.30841E-07	4.57677E-08	6.42242E-07	4.68007E-08];

i = i+1;
svd_n_data(i,:) = [0.956187115	0.342940611	0.042276736	0.023805843		0.184266903	0.999025630	0.089446021	0.077959659		0.071745162	0.019105460	0.767990749	0.732722618		0.134492433	0.043658532	0.137952717	0.988211155];
cpct_data(i,:)  = [2.46205E-04	4.54284E-05	5.64599E-04	4.86244E-05		3.31567E-04	4.58436E-05	6.63973E-04	4.97204E-05		2.45254E-04	4.50614E-05	4.70565E-04	3.65787E-05		3.30841E-04	4.57677E-05	6.42242E-04	4.68007E-05];

i = i+1;
svd_n_data(i,:) = [0.956187115	0.342940611	0.042276736	0.023805843		0.184266903	0.999025630	0.089446021	0.077959659		0.071745162	0.019105460	0.767990749	0.732722618		0.134492433	0.043658532	0.137952717	0.988211155];
cpct_data(i,:)  = [2.46205E-01	4.54284E-02	5.64599E-01	4.86244E-02		3.31567E-01	4.58436E-02	6.63973E-01	4.97204E-02		2.45254E-01	4.50614E-02	4.70565E-01	3.65787E-02		3.30841E-01	4.57677E-02	6.42242E-01	4.68007E-02];

i = i+1;
svd_n_data(i,:) = [0.980889112	0.233334866	0.010412646	0.002133049		0.106876541	0.997965855	0.053220123	0.015962974		0.014194895	0.002069055	0.956080890	0.358966760		0.055598443	0.010542971	0.055901685	0.997987993];
cpct_data(i,:)  = [9.80712E-09	1.95373E-09	2.25476E-08	1.95002E-09		1.32936E-08	1.96857E-09	2.65937E-08	1.99403E-09		9.77741E-09	1.93904E-09	1.85824E-08	1.46603E-09		1.32647E-08	1.96734E-09	2.57205E-08	1.87405E-09];

i = i+1;
svd_n_data(i,:) = [0.980889112	0.233334866	0.010412646	0.002133049		0.106876541	0.997965855	0.053220123	0.015962974		0.014194895	0.002069055	0.956080890	0.358966760		0.055598443	0.010542971	0.055901685	0.997987993];
cpct_data(i,:)  = [9.80712E-06	1.95373E-06	2.25476E-05	1.95002E-06		1.32936E-05	1.96857E-06	2.65937E-05	1.99403E-06		9.77741E-06	1.93904E-06	1.85824E-05	1.46603E-06		1.32647E-05	1.96734E-06	2.57205E-05	1.87405E-06];

i = i+1;
svd_n_data(i,:) = [0.980889112	0.233334866	0.010412646	0.002133049		0.106876541	0.997965855	0.053220123	0.015962974		0.014194895	0.002069055	0.956080890	0.358966760		0.055598443	0.010542971	0.055901685	0.997987993];
cpct_data(i,:)  = [9.80712E-03	1.95373E-03	2.25476E-02	1.95002E-03		1.32936E-02	1.96857E-03	2.65937E-02	1.99403E-03		9.77741E-03	1.93904E-03	1.85824E-02	1.46603E-03		1.32647E-02	1.96734E-03	2.57205E-02	1.87405E-03];

i = i+1;
svd_n_data(i,:) = [0.962140218	0.320850828	0.030894149	0.013509206		0.169648850	0.998401019	0.083292470	0.055441132		0.050021536	0.011725246	0.839610686	0.636609757		0.108943179	0.031688386	0.110893615	0.992223874];
cpct_data(i,:)  = [1.96737E-08	3.74676E-09	4.51498E-08	3.89428E-09		2.65506E-08	3.77905E-09	5.31476E-08	3.98210E-09		1.96029E-08	3.71739E-09	3.74879E-08	2.92887E-09		2.64927E-08	3.77393E-09	5.14057E-08	3.74629E-09];

i = i+1;
svd_n_data(i,:) = [0.962140218	0.320850828	0.030894149	0.013509206		0.169648850	0.998401019	0.083292470	0.055441132		0.050021536	0.011725246	0.839610686	0.636609757		0.108943179	0.031688386	0.110893615	0.992223874];
cpct_data(i,:)  = [1.96737E-05	3.74676E-06	4.51498E-05	3.89428E-06		2.65506E-05	3.77905E-06	5.31476E-05	3.98210E-06		1.96029E-05	3.71739E-06	3.74879E-05	2.92887E-06		2.64927E-05	3.77393E-06	5.14057E-05	3.74629E-06];

i = i+1;
svd_n_data(i,:) = [0.962140218	0.320850828	0.030894149	0.013509206		0.169648850	0.998401019	0.083292470	0.055441132		0.050021536	0.011725246	0.839610686	0.636609757		0.108943179	0.031688386	0.110893615	0.992223874];
cpct_data(i,:)  = [1.96737E-02	3.74676E-03	4.51498E-02	3.89428E-03		2.65506E-02	3.77905E-03	5.31476E-02	3.98210E-03		1.96029E-02	3.71739E-03	3.74879E-02	2.92887E-03		2.64927E-02	3.77393E-03	5.14057E-02	3.74629E-03];

N_dataSets = i;
dataSet = 1:1:N_dataSets;

%% Total Data Set
% [r,c] = size(svd_n_data)
% N_dataPoints = r*c;

svd_n_data_total = [];
cpct_data_total  = [];
for j = dataSet
    svd_n_data_total = [svd_n_data_total, svd_n_data(j,:)];
    cpct_data_total  = [cpct_data_total,  cpct_data(j,:) ];
end

%% Covariance Calc
for j = dataSet
    [mu_x] = calcMean(svd_n_data(j,:));
    [error_x] = calcError(svd_n_data(j,:),mu_x);
%     [mu_y] = calcMean(cpct_data(j,:));
%     [error_y] = calcError(cpct_data(j,:),mu_y);
    [mu_y] = calcMean(log10(cpct_data(j,:)));
    [error_y] = calcError(log10(cpct_data(j,:)),mu_y);

    [Covar] = calcPopCovar(error_x, error_y);
end
    [mu_x] = calcMean(svd_n_data_total);
    [error_x] = calcError(svd_n_data_total,mu_x);
%     [mu_y] = calcMean(cpct_data_total);
%     [error_y] = calcError(cpct_data_total,mu_y);
    [mu_y] = calcMean(log10(cpct_data_total));
    [error_y] = calcError(log10(cpct_data_total),mu_y);
    [Covar_total] = calcPopCovar(error_x, error_y)


if FLAG.PLOT
    %% Plots
%     for j = dataSet
%         figure
%         scatter(svd_n_data(j,:),log10(cpct_data(j,:)))
% 
%         %     figure
%         %     semilogy(svd_n_data(j,:),cpct_data(j,:),'o')
%     end

    figure
    %     scatter(svd_n_data(j,:),cpct_data(j,:))
    scatter(svd_n_data_total , log10(cpct_data_total))
    title('All Data')
    xlabel('Normalized Singular Values (Of SV)')
    ylabel('log10(Error Covariance)')

%% Combine 2 sets
set1 = 1; set2 = 4;
% set1 = 2; set2 = 5;
% set1 = 3; set2 = 6;
% set1 = 7; set2 = 10;
% set1 = 8; set2 = 11;
% set1 = 9; set2 = 12;

combined_SVD = [svd_n_data(set1,:),svd_n_data(set2,:)];
combined_CPC = [cpct_data(set1,:) ,cpct_data(set2,:) ];

[mu_x] = calcMean(combined_SVD);
[error_x] = calcError(combined_SVD,mu_x);
[mu_y] = calcMean(log10(combined_CPC));
[error_y] = calcError(log10(combined_CPC),mu_y);
[Covar_combined] = calcPopCovar(error_x, error_y)

figure
scatter(combined_SVD , log10(combined_CPC))
title('Combined Data')
xlabel('Normalized Singular Values (Of SV)')
ylabel('log10(Error Covariance)')


    %% Arrange Figures
    FigArrange = 1;
    if FigArrange == 1
        fig = gcf;
        NumFig = fig.Number;

        Ncol = 3;

        for i = 1:NumFig
            f = figure(i);
            k = mod(i-1,Ncol);
            row = mod(fix((i-1)/Ncol),2);
            if row == 0
                r = 575;
                %             r = 540;
            elseif row == 1
                r = 62;
            end
            f.Position = [k*575+15 r 560 420];
        end
    end
end


%% Functions
function [mu] = calcMean(x)
mu = mean(x);
end
function [error] = calcError(x,mu)
% difference between x and the expected value (mu)
error = x-mu;
end
function [Covar] = calcPopCovar(error_x, error_y)
Covar = (error_x*error_y')/length(error_x);
end