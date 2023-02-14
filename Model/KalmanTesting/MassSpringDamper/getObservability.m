%% Calculate Observability Matrix
function observability = getObservability(A,C,N_multiple)
% Nstates = length(A);
[Noutputs , Nstates] = size(C);
% CA = C;
% observability = CA;
% for i = 1:Nstates*N_multiple-1
%     CA = CA*A;
%     observability(end+1:end+Noutputs,:) = CA;
% end


CA = C*inv(A)^(Nstates*N_multiple-1);
observability = CA;
for i = 1:Nstates*N_multiple-1
    CA = CA*A;
    observability(end+1:end+Noutputs,:) = CA;
end


end