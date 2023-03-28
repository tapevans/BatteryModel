%% Hankel
% The first column in g_k is k = 0. Therefore, the first element (zero
% vector) will not be included
%
% This version handles multivariables
% g_k must be a row vector for each variable
%

function H = myHankel(g_k)
[Nvar, Nk] = size(g_k);

p = floor( (Nk+1)/2 );
H = zeros(Nvar*p,p);

for r = 1:p
    % offset = r - 1;
    offset = r;
    H((Nvar*r-(Nvar-1)) : (Nvar*r),:) = g_k(:,1+offset:p+offset);
    % I can do this with the rows or, just do end+1:end+Nvar
end

end