function S = lhsamp(m, n, a, b)
%LHSAMP  Latin hypercube distributed random numbers
%
% Call:    S = lhsamp
%          S = lhsamp(m)
%          S = lhsamp(m, n)
%
% m : number of sample points to generate, if unspecified m = 1
% n : number of dimensions, if unspecified n = m
%
% a : lower bounds
% b : upper bounds
%
% S : the generated n dimensional m sample points chosen from
%     uniform distributions on m subdivions of the interval ({a}, {b})

% hbn@imm.dtu.dk  
% Last update April 12, 2002

if nargin < 1, m = 1; end
if nargin < 2, n = m; end
if nargin > 2 & nargin < 3, a = 0*ones(n,1),b = 0*ones(n,1); end

S = zeros(n,m);
for i = 1 : n
  S0 = (rand(1, m) + (randperm(m) - 1)) / m;
  S(i, :) = a(i)*ones(1,m)+(b(i)-a(i))*S0;
end
