%% Activity Function
% C_e is the concentration of Li^+ in the electrolyte in kmol m^-3
% Also depends on temperature but for now I'll keep that constant
% T = 303.15; %[K]

function dfdc = activity_HZ(Ce, T)
    dfdc = (0.54*exp(329.0 ./ T).*(Ce.^2.0)- 0.00225*exp(1360.0 ./ T).*Ce + 0.341*exp(261.0 ./ T));
end

% A = 0.54     * exp(329.0./T);
% B = -0.00225 * exp(1360 ./T);
% % C = 0.341    * exp(261.0/T) - 1 ; %-1 doesn't show up in the C++
% C = 0.341    * exp(261.0./T);
% 
% dfdc = A .* Ce.^2 + B .* Ce + C;