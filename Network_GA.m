clc;
clear vaiables;
close all;

global N A AA Edgedata OR CR RR SO SC SR
Montecarlo=6;
OwlsHeadServiceability=ones(Montecarlo,120);
ConeyIslandServiceability=ones(Montecarlo,120);
RockawayServiceability=ones(Montecarlo,120);

for Mc=1:Montecarlo
%% Defining the network
% Adjacency matrix and nodes' characteristics are defined
A=xlsread('Nodes.xlsx','Adjacency');
Nodedata=readtable('Nodes.xlsx','Sheet','Data');
Edgedata=readtable('Nodes.xlsx','Sheet','Edges');
NProb=table2array(Nodedata(:,2));
EProb=table2array(Edgedata(:,2));
N=size(A,1);
AA=zeros(N,N);

% Network is defined based on the adjacency matrix and it nodes' attributes
G=digraph(A,Nodedata);
G.Edges.Length=table2array(Edgedata(:,1));
G.Edges.Prob=table2array(Edgedata(:,2));
G.Edges.Recovery=table2array(Edgedata(:,3));

%% Initial Operability is assigned based on probability of damage

% For nodes:
for i=1:N 
    G.Nodes.Operability(i)=randsrc(1,1,[0,1;NProb(i),1-NProb(i)]);
end
G.Nodes.Serviceability=G.Nodes.Operability;
% For edges:
for i=1:size(EProb,1) 
    G.Edges.Operability(i)=randsrc(1,1,[0,1;EProb(i),1-EProb(i)]);
end

% Damaged nodes (and their links) and links operability is changed to zero
for n=1:N
    for m=1:N  
        if G.Nodes.Operability(n)==0
           G=rmedge(G,n,m);
        end
        idxOut=findedge(G,n,m);
        if idxOut~=0 && G.Edges.Operability(idxOut)==0
           G=rmedge(G,n,m);
        end
    end
end

%% Defining Recovery matrix for each time step

% Recovery matrix showing the Serviceability of nodes in each time step
R=ones(N,2);
% Before hazard (t=1)
R(:,1)=1;
% During Hazard (t=2)
R(:,2)=G.Nodes.Serviceability;
R(:,3)=G.Nodes.Serviceability;

%% Problem Definition

FitnessFunction=@(G,V,R) MaxR(G,V,R);     % Fitness Function

nVar=20;            % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

%% GA Parameters

MaxIt=10;      % Maximum Number of Iterations

nPop=6;        % Population Size

pc=0.8;                 % Crossover Percentage
nc=2*round(pc*nPop/2);  % Number of Offsprings (Parnets)

pm=0.3;                 % Mutation Percentage
nm=round(pm*nPop);      % Number of Mutants

%% Initialization

empty_individual.Sequence=[zeros(1,12);zeros(1,8),ones(1,4)*20];
empty_individual.Resilience=[];

pop=repmat(empty_individual,nPop,1);

% Generating first generation of the population
for i=1:nPop

    V=[randperm(12);randperm(8)+12,20,20,20,20];
    
  
  
    % Initialize Sequence
    pop(i).Sequence=V;
    
    % Evaluation
    pop(i).Resilience=FitnessFunction(G,pop(i).Sequence,R);
end

% Sort Population based on best resilience
Resiliences=[pop.Resilience];
[Resiliences, SortOrder]=sort(Resiliences,'descend');
pop=pop(SortOrder);

% Array to Hold Best Sequence
BestSequence=zeros(size(V,1),size(V,2),MaxIt);
BestBestSequence=zeros(size(V,1),size(V,2),Montecarlo);
% Array to Hold Best Resilience
BestResilience=zeros(MaxIt,1);

%% Main Loop

for it=1:MaxIt
    
    % Crossover
    popc=repmat(empty_individual,nc/2,2);
    for k=1:nc/2
        
        % Select Parents Indices
        p=randsample(nPop,2);

        % Select Parents
        p1=pop(p(1));
        p2=pop(p(2));
        
        % Apply Crossover
        [popc(k,1).Sequence,popc(k,2).Sequence]=DoublePointCrossover(p1.Sequence,p2.Sequence);
        
        % Evaluate Offsprings
        popc(k,1).Resilience=FitnessFunction(G,popc(k,1).Sequence,R);
        popc(k,2).Resilience=FitnessFunction(G,popc(k,2).Sequence,R);
        
    end
    popc=reshape(popc,[],1);   % Put all the offsprings in 1 column
    
    % Mutation
    popm=repmat(empty_individual,nm,1);
    for k=1:nm
        
        % Select Parent
        i=randi([1 nPop]);
        p=pop(i);
        
        % Apply Mutation
        popm(k).Sequence=Mutation(p.Sequence);
        
        % Evaluate Mutant
        popm(k).Resilience=FitnessFunction(G,popm(k).Sequence,R);
        
    end
    
    % Create Merged Population
    popt=[pop;popc;popm];
     
    % Sort Population based on best resilience
    Resiliences=[popt.Resilience];
    [Resiliences, SortOrder]=sort(Resiliences,'descend');
    popt=popt(SortOrder);
    
    % Truncation
    pop=popt(1:nPop);
    Resiliences=Resiliences(1:nPop);
    
    % Store Best Sequence Ever Found
    BestSequence(:,:,it)=pop(1).Sequence;
    
    % Store Best Resilience Ever Found
    BestResilience(it)=pop(1).Resilience;
    
    % Show Iteration Information
    disp(['Iteration ',num2str(it),': Best Resilience = ',num2str(BestResilience(it))]);
    
end

BestResult=FitnessFunction(G,popm(1).Sequence,R);
%% Results

figure;
plot(BestResilience,'LineWidth',2);
xlabel('Iteration');
ylabel('Resilience');

disp(['Owls Head WWTP resiliency is ',num2str(OR*100),'%']);
OwlsHeadServiceability(Mc,1:size(SO,2))=SO;

disp(['Coney Island WWTP resiliency is ',num2str(CR*100),'%']);
ConeyIslandServiceability(Mc,1:size(SC,2))=SC;

disp(['Rockaway WWTP resiliency is ',num2str(RR*100),'%']);
RockawayServiceability(Mc,1:size(SR,2))=SR;

disp('Optimized sequence for best resilience is ');
disp(BestSequence(:,:,MaxIt));
BestBestSequence(:,:,Mc)=BestSequence(:,:,MaxIt);

end