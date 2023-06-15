%% Monte Carlo type Error Variance Initialization
clear all; close all; clc;

%% Inputs
% C = [1 2 ; 3 4];
% sigma = [1e-1 1e-2];
% N_states = length(sigma);

N_outputs = 10;
N_states  = 10; 
N_samples = 10000;

% C = 1: N_outputs*N_states;
C = randn(N_outputs*N_states,1);

C = reshape(C,N_outputs,[]);

sigma = logspace(-2,-1,N_states);


%% Test states distribution
% % randn(r,c) produces an r-by-c array of numbers from a standard normal
% % distribution, mean = 0, standard deviation = 1
%     r3 = randn(1e6,1);
% 
% % Plot of the histogram of the data
%     % figure
%     % h = histogram(r3);
%     % h.EdgeAlpha = 0;
% 
% % Change the number of bins to represent the data
%     % n = 7;
%     % figure
%     % histogram(r3,n) % ,n is the number of bins
% 
% % Change data to have a standard deviation of 2
%     % std = 2;
%     std = 1e-1;
%     r3 = std * randn(1,1e6);
% 
% % Plot of the histogram of the data
%     figure
%     h = histogram(r3);
%     h.EdgeAlpha = 0;
%     % h.Normalization = 'countdensity';
%     % h.Normalization = 'probability';
%     h.Normalization = 'pdf'; % Normalize the data using the probability density function estimate
% 
% % Add the plot of the normal distribution function on the histogram
%     f = @(x,mu,std) (2*pi*std^2)^(-1/2) * exp( - ((x-mu).^2) / (2*std^2) ); % Probability density function
%     x_vec = linspace(-3*std , 3*std , 100);
%     f_vec = f(x_vec,0,std);
% 
%     hold on
%     plot(x_vec,f_vec,'k')
% 
% % Calculate the variance of the data set
%     covar = data2Covar(r3)
%     std_calc = ( diag(covar) ).^(1/2)
%     std


%% Test Multiple States
% Initialize State vectors
    states_vec =   nan(N_states , N_samples);
    states_IC  = zeros(N_states , 1);
    for i = 1 : N_states
        states_vec(i,:) = sigma(i) * randn(1,N_samples);
    end

    states_vec = [states_IC states_vec];

% % Create a histogram plot for each
%     for i = 1 : N_states
%         figure
%         h = histogram(states_vec(i,:));
%         h.EdgeAlpha = 0;
%         h.Normalization = 'pdf';
%         title(['State ' num2str(i) ', \sigma = ' num2str(sigma(i))])
%     end

% Calculate the variance of the data set
    covar = data2Covar(states_vec);
    disp('covar')
    disp(num2str(covar))
    std_calc = ( ( diag(covar) ).^(1/2) )';
    disp(newline)
    disp('std_calc')
    disp(num2str(std_calc))
    disp('sigma actual')
    disp(num2str(sigma))


%% Get a predicted output error covariance
    CPCT_predicted = C * covar * C';
    disp(newline)
    disp('CPCT_predicted')
    disp(num2str(CPCT_predicted))

    CPCT_predicted_actual = C * diag(sigma.^2) * C';
    disp(newline)
    disp('CPCT_predicted_actual')
    disp(num2str(CPCT_predicted_actual))


%% Calculate Change in Output to each state individually
    y = nan(N_outputs , N_samples+1 , N_states);
    for i = 1: N_states
        x_vec      = zeros(size(states_vec));
        x_vec(i,:) = states_vec(i,:);
        y(:,:,i)   = C * x_vec;
    end

%% Calculate the output error covar
    y_reshape = [];
    for i = 1:N_states
        y_reshape = [y_reshape , y(:,:,i)];
    end
    CPCT_calc = data2Covar(y_reshape);

    % disp(newline)
    % disp('CPCT_calc (Uses all data)')
    % disp(num2str(CPCT_calc))

    % disp(newline)
    % disp('CPCT_calc * 2')
    % disp(num2str(CPCT_calc*2))

    disp(newline)
    disp('CPCT_calc (Uses all data)')
    disp('CPCT_calc * N_states')
    disp(num2str(CPCT_calc*N_states))


%% Calculate Individual Covariance
    CPCT_calc_ind = zeros(N_outputs,N_outputs,N_states+1);
    for i = 1:N_states
        CPCT_calc_ind(:,:,i)   = data2Covar(y(:,:,i));
        CPCT_calc_ind(:,:,end) = CPCT_calc_ind(:,:,end) + CPCT_calc_ind(:,:,i);
    end    

    % CPCT_calc_ind(:,:,1)
    % CPCT_calc_ind(:,:,2)
    disp(newline)
    disp('Sum of CPCT_calc_ind')
    disp(num2str(CPCT_calc_ind(:,:,end)))

    % CPCT_calc_1 = data2Covar(y(:,:,1))
    % CPCT_calc_2 = data2Covar(y(:,:,2))
    % CPCT_calc_sum = CPCT_calc_1  + CPCT_calc_2


%% Reverse Calc P_k
    CPCT = CPCT_calc_ind(:,:,end);
    C_inv  = getPinv(C);
    CT_inv = getPinv(C');

    P_reverse = C_inv * CPCT * CT_inv;
    disp(newline)
    disp('Reverse calc for P_k')
    disp(num2str(P_reverse))

    disp(newline)
    disp('sigma of Reverse calc P_k')
    disp(num2str( ( (diag(P_reverse)).^(1/2) )'  ))
    
    % disp(newline)
    disp('orig sigma')
    disp(num2str(sigma))

    disp(newline)
    disp('This works if C is square')


%% Test Generating vectors using a covariance matrix
    % mu    = [2 3];
    % Sigma = [1 1.5; 1.5 3];
    
    % mu = [0,0];
    % Sigma = [1,0
    %          0,1];
    
    mu    = [0,0];
    Sigma = [1   , 3/5
             3/5 , 2  ];
    
    n = 1000;
    rng('default')  % For reproducibility
    R = mvnrnd(mu,Sigma,n);
    R = R';

    X_min = -4.5; X_max = 4.5;
    Y_min = -4.5; Y_max = 4.5;

% Find the single pdf for each variable
    f = @(x,mu,std) (2*pi*std^2)^(-1/2) * exp( - ((x-mu).^2) / (2*std^2) ); % Probability density function

    X_vec = R(1,:); %X_min = min(X_vec); X_max = max(X_vec);
    Y_vec = R(2,:); %Y_min = min(Y_vec); Y_max = max(Y_vec);
    X_dis = linspace(X_min,X_max,50);
    Y_dis = linspace(Y_min,Y_max,50);

    f_vec_X = f(X_dis,mu(1),Sigma(1,1));
    f_vec_Y = f(Y_dis,mu(2),Sigma(2,2));

% PLOTTING OF THE DATA AND CONFIDENCE ELLIPSE
    % Determine largest and smallest eig val and vec
        [eigVec,eigVal] = eig(Sigma);
        
        [IDX_largeEigVal] = find(diag(eigVal) == max(diag(eigVal)));
        if length(IDX_largeEigVal)==1
            if IDX_largeEigVal == 1
                IDX_smallEigVal = 2;
            else
                IDX_smallEigVal = 1;
            end
        else
            IDX_largeEigVal = 1;
            IDX_smallEigVal = 2;
        end
        
        largest_eigenvec  = eigVec(:,IDX_largeEigVal);
        smallest_eigenvec = eigVec(:,IDX_smallEigVal);
        largest_eigenval  = eigVal(IDX_largeEigVal,IDX_largeEigVal);
        smallest_eigenval = eigVal(IDX_smallEigVal,IDX_smallEigVal);
    
    % Calculate the angle between the x-axis and the largest eigenvector
        angle = atan2(largest_eigenvec(2), largest_eigenvec(1));
    
    % This angle is between -pi and pi.
    % Let's shift it such that the angle is between 0 and 2pi
        if(angle < 0)
            angle = angle + 2*pi;
        end
    
    % Get the coordinates of the data mean
        avg = mean(R,2);
    
    % Get the 95% confidence interval error ellipse
        chisquare_val = 2.4477;
        % chisquare_val = 5.991;
        theta_grid = linspace(0,2*pi);
        phi = angle;
        X0=avg(1);
        Y0=avg(2);
        a=chisquare_val*sqrt(largest_eigenval);
        b=chisquare_val*sqrt(smallest_eigenval);
    
    % the ellipse in x and y coordinates 
        ellipse_x_r  = a*cos( theta_grid );
        ellipse_y_r  = b*sin( theta_grid );
    
    %Define a rotation matrix
        Rot = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
    
    %let's rotate the ellipse to some angle phi
        r_ellipse = [ellipse_x_r;ellipse_y_r]' * Rot;
        % if min(ellipse_x_r) < X_min
        %     X_min = min(ellipse_x_r);
        % end
        % if min(ellipse_y_r) < Y_min
        %     Y_min = min(ellipse_y_r);
        % end
        % if max(ellipse_x_r) > X_max
        %     X_max = max(ellipse_x_r);
        % end
        % if max(ellipse_y_r) > Y_max
        %     Y_max = max(ellipse_y_r);
        % end
    
        % eta = 0.1;
        % X_min = X_min - eta;
        % Y_min = Y_min - eta;
        % X_max = X_max + eta;
        % Y_max = Y_max + eta;
    
    % % Plot the data
    %     figure
    %     hold on
    %     plot(R(1,:) , R(2,:),'k+')
    % 
    % % Plot the eigenvectors
    %     quiver(X0, Y0,  largest_eigenvec(1)*sqrt(largest_eigenval) ,  largest_eigenvec(2)*sqrt(largest_eigenval) , '-m', 'LineWidth',2);
    %     quiver(X0, Y0, smallest_eigenvec(1)*sqrt(smallest_eigenval), smallest_eigenvec(2)*sqrt(smallest_eigenval), '-g', 'LineWidth',2);
    % 
    % % Draw the error ellipse
    %     plot(r_ellipse(:,1) + X0 , r_ellipse(:,2) + Y0,'g-','Linewidth',3)
    
    % Plot 3D Version
        figure
        hold on
        plot3(X_dis                    ,Y_max*ones(length(X_dis)),f_vec_X                   ,'b-','Linewidth',2)
        plot3(X_min*ones(length(X_dis)),Y_dis                    ,f_vec_Y                   ,'r-','Linewidth',2)
        plot3(R(1,:)                   , R(2,:)                  ,zeros(n,1)        ,'k+')
        plot3(r_ellipse(:,1) + X0      , r_ellipse(:,2) + Y0     ,zeros(length(theta_grid)) ,'g-','Linewidth',3)
        xlabel('X')
        ylabel('Y')
        xlim([X_min,X_max])
        ylim([Y_min,Y_max])
        grid on
        view([46 54])

% % Plot surf of the PDF
%     [X_mesh,Y_mesh] = meshgrid(X_dis,Y_dis);
%     for i = 1:length(X_dis)
%         for j = 1:length(Y_dis)
%             x = [X_dis(i) Y_dis(j)];
%             [pdf(i,j)] = getPDFMultiVar(x,mu,Sigma);
%         end
%     end
%     figure
%     s = surf(X_mesh,Y_mesh,pdf);
%     s.EdgeColor = 'none';
%     xlabel('X')
%     ylabel('Y')
%     zlabel('PDF')
%     xlim([X_min,X_max])
%     ylim([Y_min,Y_max])




%% Helper Function
function [X_inv] = getPinv(X)
    [U,S,V] = svd(X);
    S_inv = pinv(S);
    X_inv = V * S_inv * U';
end

% getPDFMultiVariateNormalDistribution
function [pdf] = getPDFMultiVar(x,mu,Sigma)
    d   = length(mu);
    pdf = ((2*pi)^(-d/2)) * ((det(Sigma))^(-1/2)) * exp(-0.5 * (x - mu) * inv(Sigma) * (x-mu)' );
end

