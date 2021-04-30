%%% ISCA
%Improved sine cosine algorithm with crossover scheme for global optimization
%ShubhamGupta?,KusumDeep 
%Knowledge-Based Systems 165(2019)374�406



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



function SOSCA(fun,foutput,nloop,nsol,nvar,nbit,a,b)
N=nsol;
Max_iteration=nloop;
lb=a';
ub=b';
dim=nvar;
fobj=fun;

%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);

Destination_position=zeros(1,dim);
Destination_fitness=inf;

Convergence_curve=zeros(1,Max_iteration);
Objective_values = zeros(1,size(X,1));

% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    [Objective_values(1,i),f(i),g(:,i)]=feval(fobj,X(i,:)');
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
Xbest=X(npmin,:);
%Main loop
t=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness
while t<=Max_iteration
    
    % Eq. (3.4)
    a = 2;
    Max_iteration = Max_iteration;
    r1=a-t*((a)/Max_iteration); % r1 decreases linearly from a to 0
    
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
                Xnew(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j)))+rand*(Xbest(j)-X(i,j));
            else
                % Eq. (3.2)
                Xnew(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j)))+rand*(Xbest(j)-X(i,j));
            end
            
            %%% crossover
            if rand<0.3  %% Cr=0.7
                Xnew(i,j)=X(i,j);
            end
            
        end
    end
    
    for i=1:size(Xnew,1)
         
        % Check if solutions go outside the search spaceand bring them back
        Flag4ub=Xnew(i,:)>ub;
        Flag4lb=Xnew(i,:)<lb;
        Xnew(i,:)=(Xnew(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate the objective values
        [fnew(1,i),f(i),g(:,i)]=feval(fobj,Xnew(i,:)');
        
        % Update the destination if there is a better solution
        if fnew(1,i)<Destination_fitness
            Destination_position=Xnew(i,:);
            Destination_fitness=fnew(1,i);
        end
        
        %%% greedy selection
        if fnew(1,i)<Objective_values(1,i)
            Objective_values(1,i)=fnew(1,i);
            X(i,:)=Xnew(i,:);
        end
       
    end
    
    [fbest,nbest]=min(fnew);
    Xbest=Xnew(nbest,:);
        
    
    Convergence_curve(t)=Destination_fitness;
    
%     % Display the iteration and best optimum obtained so far
%     if mod(t,50)==0
%         display(['At iteration ', num2str(t), ' the optimum is ', num2str(Destination_fitness)]);
%     end
    x=X;
    fp=Objective_values;
    [fpmin,npmin]=min(fp);
    fpminhist=[fpminhist fpmin];fpavghist=[fpavghist mean(fp)];
    fpmaxhist=[fpmaxhist max(fp)];
    fhist=[fhist f(npmin)];ghist=[ghist g(:,npmin)];
    % Increase the iteration counter
    t=t+1;
end
xmin=Destination_position';
[fpmin,fmin,gmin]=feval(fun,xmin);
maxeval=nloop*nsol;
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

