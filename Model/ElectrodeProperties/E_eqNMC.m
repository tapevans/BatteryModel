% E^eq (NMC532)
function E_eq = E_eqNMC(x)
E_eq = 5.314735633000300 ...
    - 3640.117692001490*x.^14 ...
    + 13176.57544484270*x.^13 ...
    - 14557.42062291360*x.^12 ...
    - 1571.094264365090*x.^11 ...
    + 12656.30978512400*x.^10 ...
    - 2057.808873526350*x.^9  ...
    - 10743.74333186190*x.^8  ...
    + 8698.112755348720*x.^7  ...
    - 829.7904604107030*x.^6  ...
    - 2073.765547574810*x.^5  ...
    + 1190.223421193310*x.^4  ...
    - 272.4851668445780*x.^3  ...
    + 27.23409218042130*x.^2  ...
    - 4.158276603609060*x     ...
    - (5.573191762723310e-4)*exp(6.560240842659690*x.^41.48209275061330);
end