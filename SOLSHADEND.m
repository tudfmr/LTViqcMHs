function lshade(fun,foutput,nloop,nsol,nvar,nbit,a,b)
%% ==========================================================================
% LSHADE-Neurodynamic
% K. Sallam, R. Sarker, D. Essam S. Elsayed. Neurodynamic Differential Evolution Algorithm and Solving CEC2015 Competition Problems.
% IEEE Congress on Evolutionary Computation, Sendai, Japan, 2015, in press
% parts of this code come from LSHADE source code
% Should you have any queries, please contact
% Mr. Karam Sallam.
% University of New South Wales at Canberra
% karam.sallam@student.adfa.edu.au
% ==========================================================================
format long;
format compact;
global xmin xmax n I_fno nfes result iter PND pop pop_size max_nfes all_res
n = nvar; %% problem dimension change to 30 50 and 100
% num_prbs=15; %% number of problems
I_fno=0;     %% problem number
max_nfes = nloop*nsol; %% max fitness evaluations
ND_FES=(1/2)*max_nfes; %% limit to apply ND
% totalTime = 51*num_prbs; % The total number of runs of all problems
% list_data = [0.0001 .001 0.01 .02 0.03 0.04 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]; %% to rercord the fitness values at these levels
% all_results=zeros(num_prbs*51,size(list_data,2)); %% an array to save the fitness values at each time step, i.e. 0.0001 MaxFES, 0.0001 * MaxFES, ...
% optimality=100:100:1500; %% optimal solutions known
xmax = a'; xmin = b'; %% variables limits

nfes=0; %% current fitness evaluations
m2=100; %% for resizing the memory size
%% main runs

    f_optimal=0;
    %% Record the best results
    iter=0; %% current generation
   
    %% parameter settings for L-SHADE
    p_best_rate = 0.11;
    arc_rate =1.4;
    memory_size = m2;
    PND=1; %% initial prob of ND
    min_memory_size=5; %% memory size
    pop_size=nsol; %% pop size
    max_pop_size = pop_size; %% max pop size
    min_pop_size =4; %% min pop size
    
    
    %% Initialize the main population
    popold = repmat(xmin.*ones(1,n), pop_size, 1) + rand(pop_size,n) .* (repmat(xmax.*ones(1,n) - xmin.*ones(1,n), pop_size, 1));
% %     pop = popold; % the old population becomes the current population
    
%     fitness = cec15_func(pop',I_fno);  %% fitness value
%     load(fundat)
%     popold=x0';
%     fitness=fp0;
%     fitness = fitness';

pop = popold; % the old population becomes the current population

for i=1:nsol
    [fitness(i,1),f0(i),g0(:,i)]= feval(fun,pop(i,:)');
end



    %% to save every fitness value to print all results as required--uncomment it if it is  required
     
%     all_res = fitness;
%     all_res (1:pop_size) = min(fitness);
    %% sort
    [valBest, ~] = sort(fitness, 'ascend'); %% sort fitness values
    [~,loc]=min(fitness);
    bestold=abs(f_optimal-valBest(1)); %% calculate the best value
    %% store the location of the best individual
    nfes = 0;
    for i = 1 : pop_size
        
        if nfes > max_nfes;
            break;
        end
    end
    memory_sf = 0.5 .* ones(memory_size, 1); 
    memory_cr = 0.5 .* ones(memory_size, 1);
    memory_pos = 1;
    archive.NP = arc_rate * pop_size; % the maximum size of the archive
    archive.pop = zeros(0, n); % the solutions stored in te archive
    archive.funvalues = zeros(0, 1); % the function value of the archived solutions
    %% a stoping criterion
    if nfes<max_nfes
        stop_con=0;
    end
    %% main loop
    while ~stop_con
        iter=iter+1;
        %% apply ND or LSHADE based on adaptive mechanism
        if rand<PND && nfes<ND_FES 
            [popold,fitness,z,fod]=ND(loc,popold,fitness,fun); %% calling ND
            pop = popold; % the old population becomes the current population
            [temp_fit, sorted_index] = sort(fitness, 'ascend');
        else
         
            pop = popold; % the old population becomes the current population
            [temp_fit, sorted_index] = sort(fitness, 'ascend');
            mem_rand_index = ceil(memory_size * rand(pop_size, 1));
            mu_sf = memory_sf(mem_rand_index);
            mu_cr = memory_cr(mem_rand_index);
            %% for generating crossover rate
            if rand <0.8
                % generating crossover rate using normal distribution
                cr = normrnd(mu_cr, 0.1);
                term_pos = find(mu_cr == -1);
                cr(term_pos) = 0;
                cr = min(cr, 1);
                cr = max(cr, 0);
            else
                %  generating crossover rate using cauchy distribution
                cr = mu_cr + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
                pos = find(cr <= 0);
                while ~ isempty(pos)
                    cr(pos) = mu_cr(pos) + 0.1* tan(pi * (rand(length(pos), 1) - 0.5));
                    pos = find(cr <= 0);
                end
                cr = min(cr, 1);
                cr = max(cr, 0);
            end
            %% for generating scaling factor
            sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
            pos = find(sf <= 0);
            while ~ isempty(pos)
                sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
                pos = find(sf <= 0);
            end
            sf = min(sf, 1);
            r0 = 1 : pop_size;
            popAll = [pop; archive.pop];
            [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);
            pNP = max(round(p_best_rate * pop_size), 2); %% choose at least two best solutions
            randindex = ceil(rand(1, pop_size) .* pNP); %% select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions
            vi = pop + sf(:, ones(1, n)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
            vi = boundConstraint(vi, pop);
            mask = rand(pop_size, n) > cr(:, ones(1, n)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * n)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([pop_size n], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop(mask);
            children_fitness=[];
            for i=1:size(ui,1)
                children_fitness(i,1) = feval(fun, ui(i,:)');
            end
%             children_fitness = children_fitness';
            for i = 1 : pop_size
                nfes = nfes + 1;
                if nfes > max_nfes
                    stop_con=1;
                    break;
                end
            end
            %% to save every fitness value to print all results as required-- uncomment it if it is  required
%             for i=1:size(children_fitness,1)
%                 all_res(size(all_res,1)+1)= min(min(all_res) ,children_fitness(i));
%             end
           
            dif = abs(fitness - children_fitness);
            %% I == 1: the parent is better; I == 2: the offspring is better
            I = (fitness > children_fitness);
            goodCR = cr(I == 1);
            goodF = sf(I == 1);
            dif_val = dif(I == 1);
            % isempty(popold(I == 1, :))
            archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));
            [fitness, I] = min([fitness, children_fitness], [], 2);
            popold = pop;
            popold(I == 2, :) = ui(I == 2, :);
            num_success_params = numel(goodCR);
            if num_success_params > 0
                sum_dif = sum(dif_val);
                dif_val = dif_val / sum_dif;
                %% for updating the memory of scaling factor
                memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);
                %% for updating the memory of crossover rate
                if max(goodCR) == 0 || memory_cr(memory_pos) == -1
                    memory_cr(memory_pos) = -1;
                else
                    memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
                end
                memory_pos = memory_pos + 1;
                if memory_pos > memory_size;
                    memory_pos = 1;
                end
            end
            %% for resizing the memory size from 100 to 5
            plan_pop_size2 = round((((min_memory_size - 100) / max_nfes) * nfes) + 100);
            if m2 > plan_pop_size2
                reduction_ind_num2 = m2 - plan_pop_size2;
                if m2 - reduction_ind_num2 < min_memory_size;
                    reduction_ind_num2 = m2 - min_memory_size;
                end
                m2=m2-reduction_ind_num2;
                memory_size=m2;
            end
            %% for resizing the population size
            plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * nfes) + max_pop_size);
            if pop_size > plan_pop_size
                reduction_ind_num = pop_size - plan_pop_size;
                if pop_size - reduction_ind_num < min_pop_size;
                    reduction_ind_num = pop_size - min_pop_size;
                end
                pop_size = pop_size - reduction_ind_num;
                for r = 1 : reduction_ind_num
                    [~, indBest] = sort(fitness, 'ascend');
                    worst_ind = indBest(end);
                    popold(worst_ind,:) = [];
                    pop(worst_ind,:) = [];
                    fitness(worst_ind,:) = [];
                end
                archive.NP = round(arc_rate * pop_size);
                if size(archive.pop, 1) > archive.NP
                    rndpos = randperm(size(archive.pop, 1));
                    rndpos = rndpos(1 : archive.NP);
                    archive.pop = archive.pop(rndpos, :);
                end
            end
        end
        [valBest, indBest] = sort(fitness, 'ascend');
        loc=indBest(1); %% location of best individual
        best_sol=popold(loc,:);
       result(1,iter)=min(fitness);  
        bestold= abs(f_optimal-valBest(1));
        %% a stopping criterion is met --> stop and go to the next run
        if abs(bestold<=1e-8)
            bestold=0;
            stop_con=1;
%             best_obt =f_optimal;%% just for printing -- uncomment if it is  required
%             all_res (size(all_res,1):max_nfes)= best_obt; %% just for printing -- uncomment if it is  required
        end
        if nfes > max_nfes
            stop_con=1;
        end
    end
    
 xmin=best_sol';
[fpmin,fmin,gmin]=feval(fun,xmin);
maxeval=nloop*nsol;
save(foutput,'xmin','fpmin','fmin','gmin','maxeval')  
    %% just for printing -- uncomment if it is  required
%     v_limits=max_nfes*list_data;
%     all_results(time,:) = abs(f_optimal-all_res(v_limits))';
%     all_res=[];

function vi = boundConstraint (vi, pop)
global xmin xmax 
% if the boundary constraint is violated, set the value to be the middle
% of the previous value and the bound
%
% Version: 1.1   Date: 11/20/2007
% Written by Jingqiao Zhang, jingqiao@gmail.com

[NP, D] = size(pop);  % the population size and the problem's dimension

%% check the lower bound
xl = repmat(xmin, NP, 1);
pos = vi < xl;
vi(pos) = (pop(pos) + xl(pos)) / 2;

%% check the upper bound
xu = repmat(xmax, NP, 1);
pos = vi > xu;
vi(pos) = (pop(pos) + xu(pos)) / 2;

function archive = updateArchive(archive, pop, funvalue)
% Update the archive with input solutions
%   Step 1: Add new solution to the archive
%   Step 2: Remove duplicate elements
%   Step 3: If necessary, randomly remove some solutions to maintain the archive size
%
% Version: 1.1   Date: 2008/04/02
% Written by Jingqiao Zhang (jingqiao@gmail.com)

if archive.NP == 0, return; end

if size(pop, 1) ~= size(funvalue,1), error('check it'); end

% Method 2: Remove duplicate elements
popAll = [archive.pop; pop ];
funvalues = [archive.funvalues; funvalue ];
[dummy IX]= unique(popAll, 'rows');
if length(IX) < size(popAll, 1) % There exist some duplicate solutions
  popAll = popAll(IX, :);
  funvalues = funvalues(IX, :);
end

if size(popAll, 1) <= archive.NP   % add all new individuals
  archive.pop = popAll;
  archive.funvalues = funvalues;
else                % randomly remove some solutions
  rndpos = randperm(size(popAll, 1)); % equivelent to "randperm";
  rndpos = rndpos(1 : archive.NP);
  
  archive.pop = popAll  (rndpos, :);
  archive.funvalues = funvalues(rndpos, :);

end

function [r1, r2] = gnR1R2(NP1, NP2, r0)

% gnA1A2 generate two column vectors r1 and r2 of size NP1 & NP2, respectively
%    r1's elements are choosen from {1, 2, ..., NP1} & r1(i) ~= r0(i)
%    r2's elements are choosen from {1, 2, ..., NP2} & r2(i) ~= r1(i) & r2(i) ~= r0(i)
%
% Call:
%    [r1 r2 ...] = gnA1A2(NP1)   % r0 is set to be (1:NP1)'
%    [r1 r2 ...] = gnA1A2(NP1, r0) % r0 should be of length NP1
%
% Version: 2.1  Date: 2008/07/01
% Written by Jingqiao Zhang (jingqiao@gmail.com)

NP0 = length(r0);

r1 = floor(rand(1, NP0) * NP1) + 1;
%for i = 1 : inf
for i = 1 : 99999999
    pos = (r1 == r0);
    if sum(pos) == 0
        break;
    else % regenerate r1 if it is equal to r0
        r1(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
    end
    if i > 1000, % this has never happened so far
        error('Can not genrate r1 in 1000 iterations');
    end
end

r2 = floor(rand(1, NP0) * NP2) + 1;
%for i = 1 : inf
for i = 1 : 99999999
    pos = ((r2 == r1) | (r2 == r0));
    if sum(pos)==0
        break;
    else % regenerate r2 if it is equal to r0 or r1
        r2(pos) = floor(rand(1, sum(pos)) * NP2) + 1;
    end
    if i > 1000, % this has never happened so far
        error('Can not genrate r2 in 1000 iterations');
    end
end

%% ==========================================================================
% LSHADE-Neurodynamic
% K. Sallam, R. Sarker, D. Essam S. Elsayed. Neurodynamic Differential Evolution Algorithm and Solving CEC2015 Competition Problems.
% IEEE Congress on Evolutionary Computation, Sendai, Japan, 2015, in press
% parts of this code come from LSHADE source code
% Should you have any queries, please contact
% Mr. Karam Sallam.
% University of New South Wales at Canberra
% karam.sallam@student.adfa.edu.au
% ==========================================================================
function [popold,valParents,z,fod]=ND(loc1,popold,valParents,fun)
global nfes n I_fno PND iter result 
popold1=popold;
mode=popold(loc1,:);
mode_org=popold(loc1,:);
z=zeros(20,n);
f_loc1 = valParents(loc1);
mu=exp(-1);
for t=1:20
    z(t,:)=mu.*mode +(1-mu).*(griewankd2(mode,f_loc1,fun))';
    mode=z(t,:);
end
z(size(z,1)+1,:)=mode_org;
for i=1:size(z,1)
    fod(i)=feval(fun,z(i,:)');
end
nfes=nfes+size(fod,2);
fod=fod';
%% just for printing -- uncomment if it is  required
% for i=1:size(fod,1) 
%     all_res(size(all_res,1)+1)= min(min(all_res) ,fod(i));
% end
%% ---------------------------------------------------
[~,r]=min(fod);
if iter <=2        
    PND=1;    
else
    if r==21               
        PND=0.01;        
    else                     
        diff=result(iter-2)-result(iter-1);           
        PND= max(0.01,diff/(result(iter-2)));                
    end
end
popold(loc1,:) =z(r,:);
popold = boundConstraint(popold, popold1);
valParents(loc1) = feval(fun,popold(loc1,:)');
nfes=nfes+1;
%% just for printing -- uncomment if it is  required
% for i=1:1    
%     all_res(size(all_res,1)+1)= min(min(all_res) ,valParents(loc1));  
% end
%% ------------------------------------------------------
[valParents, sorted_index] = sort(valParents, 'ascend');

popold=popold(sorted_index,:);
[fod, sorted_index1] = sort(fod, 'ascend');

z=z(sorted_index1,:);
for i=1:5
    if valParents(i)>fod(i)             
        valParents(i)=fod(i);
        popold(i,:)=z(i,:);     
    end
end

%% ==========================================================================
% LSHADE-Neurodynamic
% K. Sallam, R. Sarker, D. Essam S. Elsayed. Neurodynamic Differential Evolution Algorithm and Solving CEC2015 Competition Problems.
% IEEE Congress on Evolutionary Computation, Sendai, Japan, 2015, in press
% parts of this code come from LSHADE source code
% Should you have any queries, please contact
% Mr. Karam Sallam.
% University of New South Wales at Canberra
% karam.sallam@student.adfa.edu.au
% ==========================================================================
function dydt=griewankd2(y,f,fun)

global n xmax xmin nfes
grd=zeros(1,n);
y=y';
 x=plf(y,xmin,xmax);
delta=1e-2;
for j=1:n
%     f_dum=0;
    x_dum= x';
    x_dum(1,j)= x_dum(j)+delta;
    f_dum = feval(fun,x_dum');
    nfes=nfes+1;
    grd (1,j)= (f_dum-f)/delta; %%for computing numerical gradient
end
%% just for printing -- uncomment if it is  required
% for i=1:n 
%     all_res(size(all_res,1)+1)= min(all_res);
% end
%% ------------------------------------------------

gy=grd';
vv1= gy>xmax';
vv2=gy<xmin';
for j=1:n
    if vv1(j)==1 || vv2(j)==1
        gy(j)=(gy(j)/sum(abs(gy)))*xmax(j);
    else
        gy(j)=gy(j);
    end
end
dydt=(plf(x-gy,xmin,xmax)); %% handling boundry


%% ==========================================================================
% LSHADE-Neurodynamic
% K. Sallam, R. Sarker, D. Essam S. Elsayed. Neurodynamic Differential Evolution Algorithm and Solving CEC2015 Competition Problems.
% IEEE Congress on Evolutionary Computation, Sendai, Japan, 2015, in press
% parts of this code come from LSHADE source code
% Should you have any queries, please contact
% Mr. Karam Sallam.
% University of New South Wales at Canberra
% karam.sallam@student.adfa.edu.au
% ==========================================================================

%% this function used to handle boundry
function lumta = plf (x,lmin, lmax)
        global n
        lumta=zeros(n,1);
        for i=1:n
            if x(i) < lmin(i)
                lumta(i) = lmin(i);
            elseif x(i) > lmax(i)
                lumta(i) = lmax(i);
            else
                lumta(i) = x(i);
            end
        end



