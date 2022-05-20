function [ ts_MC_ARFIMA ] = func_mcARFIMA( input_struct )
%UNTITLED33 Summary of this function goes here
%   Detailed explanation goes here

N = input_struct.N;
alpha = input_struct.alpha;
beta = input_struct.beta;
gamma = input_struct.gamma;
delta = input_struct.delta;
d1 = input_struct.d1;
d2 = input_struct.d2;
d3 = input_struct.d3;
d4 = input_struct.d4;
rho = input_struct.rho;

eps = func_gen_corr_noise(N, 0, 0, 1, 1, rho);

X1 = func_gen_ARFIMA_0d0(N, d1, randn(N,1));
X2 = func_gen_ARFIMA_0d0(N, d2, eps(:,1));
X3 = func_gen_ARFIMA_0d0(N, d3, eps(:,2));
X4 = func_gen_ARFIMA_0d0(N, d4, randn(N,1));

X = alpha*X1 + beta*X2;
Y = gamma*X3 + delta*X4;

ts_MC_ARFIMA = [X, Y];
end

function [ eps_corr ] = func_gen_corr_noise(n, mu1, mu2, sigma1, sigma2, rho )
% Inputs:
% - n: time series length (data points)
% - mu1, mu2: means of gaussian random variables 1 and 2
% - sigma1, sigma2: variances of gaussian random variables 1 and 2
% - rho: correlation between gaussian random variables 1 and 2

eps1 = randn(n,1);
eps2 = rho*eps1 + sqrt(1-rho^2)*randn(n,1);

eps_corr = [mu1 + sigma1*eps1, mu2 + sigma2*eps2];


end

function [ X ] = func_gen_ARFIMA_0d0( N, d, eps, n_inf )
%UNTITLED31 Summary of this function goes here
%   Detailed explanation goes here

switch nargin
    case 1
        d = 0;
        eps = randn(N,1);
        n_inf = 100;
    case 2
        eps = randn(N,1);
        n_inf = 100;
    case 3
        n_inf = 100;
        if N ~= length(eps)
            error('length of eps has to be equal to n')
        end
    case 4
        if N ~= length(eps)
            error('length of eps has to be equal to n')
        end
    otherwise
        error('incorrect number of input arguments')
end

X = zeros(N,1);
for t = 1:N
    if t == 1
        xt = eps(t);
    elseif t-n_inf <= 0
        an = d*gamma((1:1:t-1)-d)./(gamma(1-d)*gamma((1:1:t-1)+1));
        xt = an*flipud(X(1:t-1)) + eps(t);
    else
        an = d*gamma((1:1:n_inf)-d)./(gamma(1-d)*gamma((1:1:n_inf)+1));
        xt = an*flipud(X(t-n_inf:t-1)) + eps(t);
    end
    X(t) = xt;
end
end



