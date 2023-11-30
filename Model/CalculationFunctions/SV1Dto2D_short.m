%% SV1Dto2D_short
function SV_out = SV1Dto2D_short(SV_in , SIM__SV_nan , N__IDX_1Dto2D)
    SV_out = SIM__SV_nan;
    SV_out(N__IDX_1Dto2D) = SV_in;
end