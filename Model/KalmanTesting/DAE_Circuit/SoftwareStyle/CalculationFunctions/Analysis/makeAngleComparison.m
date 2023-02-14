%% Angle Comparison
function makeAngleComparison(MyData,FLAG)
% Use Selected Data
    Q_0_vec_spec = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1]; % Process Noise
    R_0_vec_spec = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1]; % Measurement Noise
    Ts_vec_spec  = [1e-7 5e-7 1e-6 5e-6 1e-5 1e-4];
    C_vec_spec   = [1 2 3];

%     Q_0_vec_spec = [ 1e0 1e1]; % Process Noise
%     R_0_vec_spec = [1e-6 ]; % Measurement Noise
%     Ts_vec_spec  = [1e-4];
%     C_vec_spec   = [1];

N_sims = length(MyData);

Q_idx = zeros(N_sims,1);
R_idx = zeros(N_sims,1);
C_idx = zeros(N_sims,1);
T_idx = zeros(N_sims,1);
for QQ = Q_0_vec_spec
    for RR = R_0_vec_spec
        for  CC = C_vec_spec
            for TT = Ts_vec_spec
                for j = 1:N_sims
                    if MyData(j).Q0 == QQ
                        Q_idx(j) = 1;
                    end
                    if MyData(j).R0 == RR
                        R_idx(j) = 1;
                    end
                    if MyData(j).C == CC
                        C_idx(j) = 1;
                    end
                    if MyData(j).Ts == TT
                        T_idx(j) = 1;
                    end
                end
            end
        end
    end
end

combined = Q_idx & R_idx & C_idx & T_idx;
idx = find(combined == 1);

% Show how the angle has a relationship between S_Orm and sing_val
    angle_total = [];
    angle_min_norm_to_max_total = [];
    angle_shift_then_norm_total = [];
    PD_total    = [];
    cpct_data   = [];
    % for i = 1:length(MyData)
    for i = idx'
        for j = 1:3 %%%%Hardcoded
            % Find smallest angle
            [M,I] = min(MyData(i).deg(j,:));
            angle_min(j) = M;
            perc_diff(j) = abs(MyData(i).sing_val(j) - MyData(i).S_Orm(I,I))/MyData(i).S_Orm(I,I);
        end
        angle_min_norm_to_max = angle_min/max(angle_min);
        angle_shift = angle_min - min(angle_min);
        angle_shift_then_norm = angle_shift / max(angle_shift);

        MyData(i).angle_min = angle_min;
        MyData(i).percent_diff = perc_diff;
        angle_total = [angle_total angle_min];
        angle_min_norm_to_max_total = [angle_min_norm_to_max_total angle_min_norm_to_max];
        angle_shift_then_norm_total = [angle_shift_then_norm_total angle_shift_then_norm];
        PD_total    = [PD_total perc_diff];
        if FLAG.UseNormCPC
            cpct_data = [cpct_data ,  MyData(i).CPC_normalized ];
        elseif FLAG.UseMAXNormCPC
            cpct_data = [cpct_data ,  MyData(i).CPC_normMAX ];
        else
            cpct_data = [cpct_data ,  MyData(i).CPC ];
        end
    end
    
%     figure
%     semilogy(angle_total ,PD_total,'ko','MarkerFaceColor','k')

% Plot CPC values
% figure
% % hold on
% semilogy(angle_total(1:4),cpct_data(1:4),'ko','MarkerFaceColor','g')
% semilogy(angle_total(5:8),cpct_data(5:8),'ko','MarkerFaceColor','r')
% xlabel('Angle')
% ylabel('CPC')
% title("Angle CPC' Comparison")
% 
% figure
% % hold on
% semilogy(angle_min_norm_to_max_total(1:4),cpct_data(1:4),'ko','MarkerFaceColor','g')
% semilogy(angle_min_norm_to_max_total(5:8),cpct_data(5:8),'ko','MarkerFaceColor','r')
% xlabel('Angle')
% ylabel('CPC')
% title("Angle CPC' Comparison Norm")
% 
% figure
% % hold on
% semilogy(angle_shift_then_norm_total(1:4),cpct_data(1:4),'ko','MarkerFaceColor','g')
% semilogy(angle_shift_then_norm_total(5:8),cpct_data(5:8),'ko','MarkerFaceColor','r')
% xlabel('Angle')
% ylabel('CPC')
% title("Angle CPC' Comparison Shift Then Norm")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% hold on
semilogy(angle_total,cpct_data,'ko','MarkerFaceColor','k')
xlabel('Angle')
ylabel('CPC')
title("Angle CPC' Comparison")

figure
% hold on
semilogy(angle_min_norm_to_max_total,cpct_data,'ko','MarkerFaceColor','k')
xlabel('Angle')
ylabel('CPC')
title("Angle CPC' Comparison Norm")

figure
% hold on
semilogy(angle_shift_then_norm_total,cpct_data,'ko','MarkerFaceColor','k')
xlabel('Angle')
ylabel('CPC')
title("Angle CPC' Comparison Shift Then Norm")
end