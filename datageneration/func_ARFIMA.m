function [ X ] = func_ARFIMA( n, d, eps, n_inf )
%UNTITLED31 Summary of this function goes here
%   Detailed explanation goes here

switch nargin
    case 1
        d = 0;
        eps = randn(n,1);
        n_inf = 100;
    case 2
        eps = randn(n,1);
        n_inf = 100;
    case 3
        n_inf = 100;
        if n ~= length(eps)
            error('length of eps has to be equal to n')
        end
    case 4
        if n ~= length(eps)
            error('length of eps has to be equal to n')
        end
    otherwise
        error('incorrect number of input arguments')
end

X = zeros(n,1);
for t = 1:n
    if t-n_inf <= 0
        an = gamma((0:1:t-1)+d)./(gamma((0:1:t-1)+1)*gamma(d));
        xt = an*flipud(eps(1:t));
    else
        an = gamma((0:1:n_inf-1)+d)./(gamma((0:1:n_inf-1)+1)*gamma(d));
        xt = an*flipud(eps(t-n_inf+1:t));
    end
    X(t) = xt;
end
end

