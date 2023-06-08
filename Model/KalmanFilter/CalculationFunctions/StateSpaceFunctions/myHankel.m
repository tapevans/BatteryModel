%% Hankel
% The first column in g_k is k = 0. Therefore, the first element (zero
% vector) will not be included
%
% This version handles multivariables
% g_k must be a row vector for each variable
%
% For a single set of data, a Hankel matrix is square (p x p)
%   * 
% For a multi  set of data, a Hankel matrix is tall (Nvar*p x p)

function H = myHankel(g_k)
[Nvar, Nk] = size(g_k);
% Since the first column in g_k is k = 0, the first element (zero vector) will not be included
Nk = Nk - 1;

p = floor( (Nk+1)/2 );
H = zeros(Nvar*p,p);

for r = 1:p
    offset = r;
    H((Nvar*r-(Nvar-1)) : (Nvar*r),:) = g_k(:,1+offset:p+offset);
    % I can do this with the rows or, just do end+1:end+Nvar
end

end