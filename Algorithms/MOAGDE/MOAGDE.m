function MOAGDE(Global)
% <algorithm> <A-G>
% fname problem ad?
% D dimension, variable say齭?
% n_obj ama? fonksiyon say齭?
% L = lbArray
% H = ubArray
% *************************** %
% ** ALGORITHM扴 VARIABLES ** %
% *************************** %
% GEN = ceil(maxIteration/NP);
% L = lbArray;
% H = ubArray;
GEN = 50;
NP=Global.N;
n_obj=Global.M;
D=Global.D;
L=Global.lower;
H=Global.upper;
X = zeros(D,1); % trial vector
Pop = zeros(D,NP); % population
Fit = zeros(NP,n_obj); % fitness of the population
r = zeros(3,1); % randomly selected indices
Archive_member_no = 0;
bestFitness=inf*ones(1,n_obj);
bestSolution=inf*ones(1,D);
Archive_X=zeros(NP,D);
Archive_F=ones(NP,n_obj)*inf;
ArchiveMaxSize = NP;
% *********************** %
% ** CREATE POPULATION ** %
% *********************** %
Population = Global.Initialization();
Pop=Population.decs;
Pop=Pop';
%[~, iBest]=min(Fit);
% ****************** %
% ** OPTIMIZATION ** %
% ****************** %

Cr_All=zeros(1,2);
NW=zeros(1,2);
global igd;
    igd=0.1;
    metric1=[];
    Score1=[];
    eva=[];
    PopObj=[];
    PopPF=[];
    
    global hvd;
    hvd = 0.1;
    metric2 =[];
    T=[];
    t=0; 
while Global.NotTermination(Population) % for each generation
    Population=INDIVIDUAL(Pop');
    OBJ=Population.objs;
    Global.Current;
    if Change(Global)
        t=t+1
        if t>=2      %除去第一次环境变化
                 score1=IGD(Population.objs,Global.PF)     %计算IGD
                 metric1=[metric1,score1];
                 migd=sum(metric1)/length(metric1); %计算MIGD
                 igd=migd;                         %将MIGD通过全局变量igd传给指标函数MIGD()
                 T=[T,t];
                 Score1=[Score1,log10(score1)];
                 
                 score2=HV(Population.objs,Global.PF)     %计算HV
                 metric2=[metric2,score2];
                 mhvd=sum(metric2)/length(metric2); %计算MGD
                 hvd=mhvd;                         %将Mhvd通过全局变量hvd传给指标函数HV() 
        end
        Population=Global.Initialization();
    else
    for i=1:NP %Calculate all the objective values first
        Fit(i,:)=OBJ(i,:);
        if dominates(Fit(i,:),bestFitness)
            bestFitness=Fit(i,:);
            bestSolution=Pop(:,i)';
        end
    end
    
    [Archive_X, Archive_F, Archive_member_no]=UpdateArchive(Archive_X, Archive_F, Pop', Fit, Archive_member_no);
    
    if Archive_member_no>ArchiveMaxSize
        Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, n_obj);
        [Archive_X, Archive_F, ~, Archive_member_no]=HandleFullArchive(Archive_X, Archive_F, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize);
    else
        Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, n_obj);
    end
    Archive_mem_ranks=RankingProcess(Archive_F, Archive_member_no, n_obj);
    allRank=RankingProcess(Fit, Pop', n_obj);
    index=RouletteWheelSelection(1./Archive_mem_ranks);
    if index==-1
        index=1;
    end
    % Update the best Organism
    bestSolution=Archive_X(index,:);
    bestFitness=Archive_F(index,:);
    
    
    CrPriods_Index=zeros(1,NP);
    Sr=zeros(1,2);
    CrPriods_Count=zeros(1,2);
    for j = 1:NP % for each individual
         
         %%%%%%%%ADAPTIVE CR RULE  %%%%%%%%%%%%%%%%%%%%%%%%%
            Ali = rand;
            if(Global.gen<=1) % Do for the first Generation
                if (Ali<=1/2)
                    CR=0.05+0.1*rand(1,1);
                    CrPriods_Index(j)=1;
                
                else
                    CR=0.9+0.1*rand(1,1);
                    CrPriods_Index(j)=2;    
                end
                CrPriods_Count(CrPriods_Index(j))=CrPriods_Count(CrPriods_Index(j)) + 1;
            else
                 if (Ali<=NW(1))
                    CR=0.05+0.1*rand(1,1);
                    CrPriods_Index(j)=1;
                
                 else
                    CR=0.9+0.1*rand(1,1);
                    CrPriods_Index(j)=2;    
                end
                CrPriods_Count(CrPriods_Index(j))=CrPriods_Count(CrPriods_Index(j)) + 1;
            end

            %%%%%%%%%%%%%%%%%END OF CR RULE%%%%%%%%%%%%%%%%%%%%%%%%%%%
                f=allRank;
                [~, in]=sort(f,'ascend');
                AA=[in(1) in(2) in(3)  in(4) in(5)  ];
                BB=[in(46) in(47) in(48) in(49) in(50) ];
                CC=[in(6) in(7) in(8) in(9) in(10) in(11) in(12) in(13) in(14) in(15) in(16) in(17) in(18) in(19) in(20) in(21) in(22) in(23) in(24) in(25) in(26) in(27) in(28) in(29) in(30) in(31) in(32) in(33) in(34) in(35) in(36) in(37) in(38) in(39) in(40) in(41) in(42) in(43) in(44) in(45)];
               
            % choose three random individuals from population,
            % mutually different 

           paraIndex=floor (rand(1,1)*length(AA))+1;
                paraIndex1=floor (rand(1,1)*length(BB))+1;
                paraIndex2=floor (rand(1,1)*length(CC))+1;
                r(1) = AA(paraIndex);
                r(2) = BB(paraIndex1);
                r(3) = CC(paraIndex2);

                F=0.1+0.9*rand(1,1);
               
                Rnd = floor(rand()*D) + 1;
                for i = 1:D
                    if ( rand()<CR ) || ( Rnd==i )
                        X(i)=Pop(i,r(3))+F*(Pop(i,r(1))-(Pop(i,r(2))));
                    else
                        X(i) = Pop(i,j);
                    end
                end
        % end%end of All cases
        % verify boundary constraints
        % verify boundary constraints
        for i = 1:D
            if (X(i)<L(i))||(X(i)>H(i))
                X(i) = L(i) + (H(i)-L(i))*rand();
            end
        end
        % select the best individual
        % between trial and current ones
        % calculate fitness of trial individual 
        PPP=INDIVIDUAL(X');
        
        XPf=PPP.objs;
        % if trial is better or equal than current
        % CRRatio(find(A==CRs(j)))=CRRatio(find(A==CRs(j)))+1-(min(f,Fit(j))/max(f,Fit(j)));
        Sr (CrPriods_Index(j)) = Sr(CrPriods_Index(j)) +1;

        Pop(:,j) = X; % replace current by trial
        Fit(j,:) = XPf ;
        % if trial is better than the best
    end
    CrPriods_Count(CrPriods_Count==0)=0.0001;
    Sr=Sr./CrPriods_Count;
%%%%%%%%%%%%%%%%USING SR ONLY%%%%%%%%%%5    
    if(sum(Sr)==0)
        W=[1/2 1/2];
    else
        W=Sr/sum(Sr);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%5
    NW=(NW*(Global.gen-1)+W)/Global.gen;
    Cr_All=Cr_All+CrPriods_Count;
    
    end
end
end

function isture = Change(Global)
    isture= rem(Global.Current,Global.tao_n) == 0;
end