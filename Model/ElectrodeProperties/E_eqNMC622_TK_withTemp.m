% E^eq (NMC622)
% Only valid between x = [0.121427623511542 , 0.998958038269834]

function E_eq = E_eqNMC622_TK_withTemp(x,T)
    % x = reshape(x,1,[]);

    P.order = 13;
    P.A = P.order+2;
    P.B = P.A + 1;
    P.C = P.B + 1;

    coeff = [     5.83523164478016
                -39.4968532413686
                450.392568164418
              -2896.65449152598
              10948.1743192078
             -24093.2683867470
              26277.2065772594
                -76.2714790742666
             -27122.6324578352
                120.545955153735
              63900.3707340657
             -83428.5358467504
              45476.9384587984
              -9519.07288483844
                 -0.0166586312926941
                  2.81001403558238
                286.250353781549];

    E_eq = zeros(1,length(x));
    for i = 1:P.order+1
        E_eq = E_eq + coeff(i) * x .^ (i-1);
    end
    E_eq = E_eq + coeff(P.A) * exp(coeff(P.B) * x .^ coeff(P.C));

    idx_low  = find(x < 0.121427623511542);
    idx_high = find(x > 0.998958038269834);
    if ~isempty(idx_low)
        nan_vec = 3*ones(length(idx_low),1);
        E_eq(idx_low) = nan_vec;
        % nan_vec = nan(length(idx_low),1);
        % E_eq(idx_low) = nan_vec;
    end
    if ~isempty(idx_high)
        nan_vec = 5*ones(length(idx_high),1);
        E_eq(idx_high) = nan_vec;
        % nan_vec = nan(length(idx_high),1);
        % E_eq(idx_high) = nan_vec;
    end
    
    % Temperature Adjustments
    p = [-0.007136 , 0.01314 , -0.003029 , -0.007028 , 0.004858 , -0.00104 , 2.823e-05];
    dUdT = p(1)*x.^6 + p(2)*x.^5 + p(3)*x.^4 + p(4)*x.^3 + p(5)*x.^2 + p(6)*x + p(7);
    T_ref  = 25 + 273.15;
    E_eq   = E_eq + (T - T_ref) .* dUdT;
end