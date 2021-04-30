%  Sine Cosine Algorithm (SCA)  
%
%  Source codes demo version 1.0                                                                      
%                                                                                                     
%  Developed in MATLAB R2011b(7.13)                                                                   
%                                                                                                     
%  Author and programmer: Seyedali Mirjalili                                                          
%                                                                                                     
%         e-Mail: ali.mirjalili@gmail.com                                                             
%                 seyedali.mirjalili@griffithuni.edu.au                                               
%                                                                                                     
%       Homepage: http://www.alimirjalili.com                                                         
%                                                                                                     
%  Main paper:                                                                                        
%  S. Mirjalili, SCA: A Sine Cosine Algorithm for solving optimization problems
%  Knowledge-Based Systems, DOI: http://dx.doi.org/10.1016/j.knosys.2015.12.022
%_______________________________________________________________________________________________
% You can simply define your cost function in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of iterations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single numbers

% To run SCA: [Best_score,Best_pos,cg_curve]=SCA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%______________________________________________________________________________________________


function SOSCA(fun,foutput,nloop,nsol,nvar,nbit,a,b,gammaUBmh)
N=nsol;
Max_iteration=nloop;
lb=a';
ub=b';
dim=nvar;
fobj=fun;

%Initialize the set of random solutions
% X=initialization(N,dim,ub,lb);
XLHS= lhsamp(N*100, dim, lb', ub');
[L,C] = kmeans(XLHS,N);
X=C';

Destination_position=zeros(1,dim);
Destination_fitness=inf;

Convergence_curve=zeros(1,Max_iteration);
Objective_values = zeros(1,size(X,1));

% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    [Objective_values(1,i),f(i),g(:,i)]=feval(fobj,X(i,:)',gammaUBmh);
    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    end
    
    All_objective_values(1,i)=Objective_values(1,i);
end
fp=Objective_values;
[fpmin,npmin]=min(fp);
fpminhist=min(fp);fpavghist=mean(fp);fpmaxhist=max(fp);
fhist=f(npmin);ghist=g(:,npmin);



%Main loop
t=1; % start from the second iteration since the first iteration was dedicated to calculating the fitness
neval=0;
maxeval=Max_iteration*N;
iimprove=0;
pop_size=N;
min_pop_size=5;
max_pop_size=N;
while neval<Max_iteration*N
    
    % Eq. (3.4)
    a = 2;
    r1=a-neval*((a)/maxeval); % r1 decreases linearly from a to 0
    
    % Update the position of solutions with respect to destination
    for i=1:size(X,1) % in i-th solution
        for j=1:size(X,2) % in j-th dimension
            
            % Update r2, r3, and r4 for Eq. (3.3)
            r2=(2*pi)*rand();
            r3=2*rand;
            r4=rand();
            
            % Eq. (3.3)
            if r4<0.5
                % Eq. (3.1)
                X(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j)));
            else
                % Eq. (3.2)
                X(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j)));
            end
            
        end
    end
    iimprove=iimprove+1;
%     [rgramaU,cgrammaU]=find(Objective_values<gammaUBmh);
%     if isempty(cgrammaU)
%     else
%         gammaUBmh=max(Objective_values(1,cgrammaU));
%     end
    for i=1:size(X,1)
         
        % Check if solutions go outside the search spaceand bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        

        
        % Calculate the objective values
        [Objective_values(1,i),f(i),g(:,i)]=feval(fobj,X(i,:)',Objective_values(1,i));
        neval=neval+1;
        % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
            iimprove=0;
        end
        if Destination_fitness <1e-31; break; end
    end
    
    if iimprove>=10
        lbold=lb;
        ubold=ub;
        aU=ub<=1e8;
        if max(aU)
            ub(aU)=ub(aU)*10;
            pop_size02=10;
            XLHS=[];
            while size(XLHS,2)==0
                XLHS= lhsamp(pop_size02*10, dim, lb', ub');
                xoutL=XLHS>=lbold';
                xoutU=XLHS<=ubold';
                XLHS(:,all([xoutL;xoutU]))=[];
                
            end
            
            if size(XLHS,2)>pop_size02
                [L,C] = kmeans(XLHS,pop_size02);
                uiAdd = C';
            else
                uiAdd=XLHS';
            end
            
            for i=1:pop_size02
                [~,nmax]=max(Objective_values);
                [children_fitness02(i,1),f1(i),g1(:,i)] = feval(fobj, uiAdd(i,:)',gammaUBmh);
                neval=neval+1;
                if children_fitness02(i) < Objective_values(nmax)
                    Objective_values(nmax) = children_fitness02(i);
                    X(nmax,:) = uiAdd(i, :);
                    iimprove=0;
                end
                if children_fitness02(i) < Destination_fitness
                    Destination_fitness = children_fitness02(i);
                    Destination_position = uiAdd(i, :);
                end
            end
        end
    end
    
        %% for resizing the population size
    plan_pop_size = round((((min_pop_size - max_pop_size) / maxeval) * neval) + max_pop_size);
    
    if  pop_size > plan_pop_size
        reduction_ind_num = pop_size - plan_pop_size;
        if pop_size - reduction_ind_num <  min_pop_size; reduction_ind_num = pop_size - min_pop_size;end
        pop_size = pop_size - reduction_ind_num;
        for r = 1 : reduction_ind_num
            [valBest indBest] = sort(Objective_values, 'ascend');
            worst_ind = indBest(end);
            X(worst_ind,:) = [];
            Objective_values(worst_ind) = [];
        end
    end   
    
    
    

    Convergence_curve(t)=Destination_fitness;
    
    % Display the iteration and best optimum obtained so far
%     if mod(t,50)==0
%         display(['At iteration ', num2str(t), ' the optimum is ', num2str(Destination_fitness)]);
%     end
    fp=Objective_values;x=X;
    [fpmin,npmin]=min(fp);
    fpminhist=[fpminhist fpmin];fpavghist=[fpavghist mean(fp)];
    fpmaxhist=[fpmaxhist max(fp)];
    fhist=[fhist f(npmin)];ghist=[ghist g(:,npmin)];
    % Increase the iteration counter
    
    t=t+1;
    
end
xmin=Destination_position';
[fpmin,fmin,gmin]=feval(fun,xmin,gammaUBmh);
maxeval=neval;
save(foutput,'xmin','fpmin','fmin','gmin','maxeval',...
    'fpminhist','fpavghist','fpmaxhist','fhist','ghist')
% figure(1),clf,hold on
% plot((1:length(besthist))*nloop*nsol/length(besthist),besthist,'r')
% plot((1:length(besthist))*nloop*nsol/length(besthist),avghist,'b')
% plot((1:length(besthist))*nloop*nsol/length(besthist),worsthist,'g')
%  Sine Cosine Algorithm (SCA)  
%
%  Source codes demo version 1.0                                                                      
%                                                                                                     
%  Developed in MATLAB R2011b(7.13)                                                                   
%                                                                                                     
%  Author and programmer: Seyedali Mirjalili                                                          
%                                                                                                     
%         e-Mail: ali.mirjalili@gmail.com                                                             
%                 seyedali.mirjalili@griffithuni.edu.au                                               
%                                                                                                     
%       Homepage: http://www.alimirjalili.com                                                         
%                                                                                                     
%  Main paper:                                                                                        
%  S. Mirjalili, SCA: A Sine Cosine Algorithm for solving optimization problems
%  Knowledge-Based Systems, DOI: http://dx.doi.org/10.1016/j.knosys.2015.12.022

% This function creates the first random population of moths

function X=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    X=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end

