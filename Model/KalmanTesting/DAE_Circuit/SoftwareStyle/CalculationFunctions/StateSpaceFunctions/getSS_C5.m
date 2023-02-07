function [C_des, C_m] = getSS_C5(SIM,N,P,FLAG)
    C_des = eye(5);
    switch FLAG.C_mode
        case 1
            C_m = C_des(1,:);
        case 2
            C_m = C_des(2,:);
        case 3
            C_m = C_des(3,:);
        case 4
            C_m = C_des;
    end
end
