%***************************************************************************************
%Author: Amir Ehsan Ranginkaman
%Last Edited:  Sunday - 2014 12 May
%Email: ranginkaman@qiau.ac.ir; aeranginkaman@gmail.com
%Reference: Unpublished work
%DOI:
%*****************************************************************************************
function SOSCA(fun,foutput,nloop,nsol,nvar,nbit,a,b)

format long;
format compact;

% Set the number of variables
D = nvar;

% Set the population size
popsize = nsol; % 40


freq = 0.25;



lu(:,1) = a;
lu(:,2) = b;

% Record the best results
outcome = [];





rand('seed', sum(100 * clock));

% Initialize the population
p = repmat(lu(:, 1), 1, popsize) + rand(D, popsize) .* (repmat(lu(:, 2) - lu(:, 1), 1, popsize));

% Evaluate the population
for i=1:size(p,2)
    [fit(1,i),f(i),g(:,i)]=feval(fun,p(:,i));
end

fp=fit;
[fpmin,npmin]=min(fp);
fpminhist=min(fp);fpavghist=mean(fp);fpmaxhist=max(fp);
fhist=f(npmin);ghist=g(:,npmin);
Xbest=p(:,npmin);

% Compute the number of FES
FES = 0;
iter = 0;

MAX_ITER = nloop;

best_error = [];
while FES < MAX_ITER*popsize
    
    F = 0.5 * (sin(2*pi*freq*iter)*(iter/MAX_ITER)+1);
    CR = 0.5 * (sin(2*pi*freq*iter+pi)*(iter/MAX_ITER)+1);
    
    % V: the set of mutant vectors
    % U: the set of trial vectors
    % fit_U: the set of objective function values of trial vectors
    V = p;
    U = p;
    fit_U = fit;
    
    % Get indices for mutation
    [r1, r2, r3] = getindex_rand_1(popsize);
    
    
    
    % Implement DE/rand/1 mutation
    V = p(:, r1) + F * (p(:, r2) - p(:, r3));
    
    % Check whether the mutant vector violates the boundaries or not
    V = repair(V, lu);
    
    % Implement binomial crossover
    for i = 1:popsize
        
        j_rand = floor(rand * D) + 1;
        t = rand(D, 1) < CR;
        t(j_rand, 1) = 1;
        t_ = 1 - t;
        U(:, i) = t .* V(:, i) + t_ .* p(:, i);
        
        [fit_U(i),fu(i),gu(:,i)] = feval(fun,U(:,i));
        FES = FES +1;
         if fit_U(i) <= fit(i)
            p(:, i) = U(:, i);
            fit(i) = fit_U(i);
            f(i)=fu(i);
            g(:,i)=gu(:,i);
        end
    end
 
    
    fp=fit;
    [fpmin,npmin]=min(fp);
    fpminhist=[fpminhist fpmin];fpavghist=[fpavghist mean(fp)];
    fpmaxhist=[fpmaxhist max(fp)];
    fhist=[fhist f(npmin)];ghist=[ghist g(:,npmin)];
    Xbest=p(:,npmin);
    iter = iter + 1;
end

xmin=Xbest;
[fpmin,fmin,gmin]=feval(fun,xmin);
maxeval=FES;
save(foutput,'xmin','fpmin','fmin','gmin','maxeval',...
    'fpminhist','fpavghist','fpmaxhist','fhist','ghist')

function [r1, r2, r3]  = getindex_rand_1(popsize)

r1 = zeros(1, popsize);
r2 = zeros(1, popsize);
r3 = zeros(1, popsize);

for i = 1 : popsize
    
    sequence = 1 : popsize;
    sequence(i) = [];

    temp = floor(rand * (popsize - 1)) + 1;
    r1(i) = sequence(temp);
    sequence(temp) = [];

    temp = floor(rand * (popsize - 2)) + 1;
    r2(i) = sequence(temp);
    sequence(temp) = [];

    temp = floor(rand * (popsize - 3)) + 1;
    r3(i) = sequence(temp);

end

function V = repair(V, lu)

[D, NP]  = size(V);

xl = repmat(lu(:, 1), 1, NP);
xu = repmat(lu(:, 2), 1, NP);

% if any variable of the mutant vector violates the lower bound
pos = V < xl;
V(pos) = 2 .* xl(pos) - V(pos);
pos_ = V(pos) > xu(pos); 
V(pos(pos_)) = xu(pos(pos_));

% if any variable of the mutant vector violates the upper bound
pos = V > xu;
V(pos) = 2 .* xu(pos) - V(pos);
pos_ = V(pos) < xl(pos); 
V(pos(pos_)) = xl(pos(pos_));