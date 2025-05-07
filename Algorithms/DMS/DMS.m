function DMS(Global)
% <algorithm> <A>
% Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    global igd;
    igd=0.1;
    metric=[];
    %% Generate random population
%     Global.NT=10;
    Population = Global.Initialization();
    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N);
    %% Optimization
    %% Set time step t
    t=0;             %时间序列
    
    m=10;            %每个划分 Gra个体数量
    OP=Global.N/m;           %划分
    
    r1=0.5;
    r2=0.5;
    CenterPoint=[];  % 历史中心点
    Zlow=[];         % 历史极小点
    Zhigh=[];        % 历史极大点 
    PF=Global.PF;
    while Global.NotTermination(Population)
%         Cc=INDIVIDUAL(sum(Population.decs,1)/Global.N);
%         DF(func2str(Global.problem),Global.NT);    %画出对应PF
%         Draw(Global.PF,'ok','Markeredgecolor',[.9 .0 .0],'Markerfacecolor',[.9 .0 .0],'MarkerSize',2);
%         Draw(Cc.objs,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',5);
%         Draw(Pt1.objs,'k+','Markeredgecolor',[.6 .0 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',3);
        
        if Change(Global)  % 检测环境变化
            t=t+1 
            Decs=Population.decs;
            Population=INDIVIDUAL(Decs);
            
            Pt2=Population;
            [FrontNo,MaxFNo] = NDSort(Pt2.objs,Global.N);
            First = find(FrontNo==1);   
            Ct=sum(Pt2(First).decs,1)/length(First);  % 计算非支配层中心点
            CenterPoint=[CenterPoint;Ct];       % 将Ct加进历史中心点
            low = min(Population.decs,[],1);        %极小点
            high  = max(Population.decs,[],1);   %极大点
            Zlow =[Zlow;low];         
            Zhigh =[Zhigh;high];
            
            P1=GraSearch(Population,t,Zlow,Zhigh,m,OP,Global,CenterPoint);                       %Gra Search
            P2=prediction(Population,t,CenterPoint,r1*Global.N,Global);               %预测
            P3=DM(Zlow,Zhigh,r2*Global.N,t,CenterPoint,Global);                        %  random diversity maintenance strategy
            Pachive=[P1,P2,P3];

            
            [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Pachive,Population],Global.N);
            PF=Global.PF;
            if t>=2      %除去第一次环境变化
                 score=IGD(Population.objs,Global.PF);      %计算IGD
                 metric=[metric,score];
                 migd=sum(metric)/length(metric); %计算MIGD
                 igd=migd;                         %将MIGD通过全局变量igd传给指标函数MIGD()
            end
        else
            MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
            Offspring  = Global.Variation(Population(MatingPool));
           [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
        end
    end
end
function  P1 = GraSearch(Pop,t,Zlow,Zhigh,m,OP,Global,CenterPoint)    
    lowD=zeros(1,Global.D);
    highD=zeros(1,Global.D);
    if t>=2
        lowD=Zlow(t,:)-Zlow(t-1,:);   %极小点移动向量
        highD=Zhigh(t,:)-Zhigh(t-1,:); %极大点移动向量
    end
    if  norm(lowD) >  norm(highD)
        D=lowD;
    else
        D=highD;
    end
    
    [CrowdDis,index] = sort(CrowdingDistance1(Pop.decs),'descend');   %决策空间拥挤度距离
    Ob=Pop(index(1:m));           %找出POP中拥挤度最大的m个个体作为种群Ob
    Popb=Ob;                      
    Ob=Ob.decs;
    for k = 1:OP-1 
       for i =1:m 
            temp=Ob(i,:)+(D/OP)*k;
            for j = 1:Global.D
                if  temp(j)<Global.lower(j)
                    temp(j)=0.5*(CenterPoint(t,j)+Global.lower(j));
                end
                if  temp(j)>Global.upper(j)
                    temp(j)=0.5*(CenterPoint(t,j)+Global.upper(j));
                end  
            end
            p=INDIVIDUAL(temp);
            Popb=[Popb,p];
       end
    end
    P1=Popb;
end
function P2 = prediction(Pop,t,CenterPoint,OP,Global)          %预测              
     [CrowdDis,index] = sort(CrowdingDistance1(Pop.decs),'descend');  %决策空间拥挤度距离
     Ob=Pop(index(1:OP));                                             %找出POP中拥挤度最大的OP个个体作为种群O
     if t>=2
        Decs=Ob.decs+(CenterPoint(t,:)-CenterPoint(t-1,:));               %预测
        for i = 1 : OP              %修正decs
            for j =1 : Global.D
                  if  Decs(i,j)<Global.lower(j)
                          Decs(i,j)=0.5*(CenterPoint(t,j)+Global.lower(j));
                  end
                  if   Decs(i,j)>Global.upper(j)
                          Decs(i,j)=0.5*(CenterPoint(t,j)+Global.upper(j));
                  end  
            end
        end
        P2=INDIVIDUAL(Decs);
     else
        P2=INDIVIDUAL(Ob.decs);
     end
end
function P3 = DM(Zlow,Zhigh,n,t,CenterPoint,Global)                       %random diversity maintenance strategy
    lowD=zeros(1,Global.D);
    highD=zeros(1,Global.D);
    if t>=2
        lowD=Zlow(t,:)-Zlow(t-1,:);   %极小点移动向量
        highD=Zhigh(t,:)-Zhigh(t-1,:); %极大点移动向量
    end
    nextlow=Zlow(t)+lowD;
    nexthigh=Zhigh(t)+highD;
    D=nexthigh-nextlow;
    P3=[];
    for i = 1:n
        Decs=[];
        for j = 1:Global.D
            temp= nextlow(j)+rand(1,1)*D(j);
            if  temp<Global.lower(j)
                temp=0.5*(CenterPoint(t,j)+Global.lower(j));
            end
            if  temp>Global.upper(j)
                temp=0.5*(CenterPoint(t,j)+Global.upper(j));
            end  
            Decs=[ Decs,temp];
        end
        P=INDIVIDUAL(Decs);
        P3=[P3,P];
    end
end

function CrowdDis = CrowdingDistance1(PopObj)
% Calculate the crowding distance of each solution in the same front

    [N,M]    = size(PopObj);
    
    CrowdDis = zeros(1,N);
    Fmax     = max(PopObj,[],1);
    Fmin     = min(PopObj,[],1);
    for i = 1 : M
        [~,rank] = sortrows(PopObj(:,i));
        CrowdDis(rank(1))   = inf;
        CrowdDis(rank(end)) = inf;
        for j = 2 : N-1
            CrowdDis(rank(j)) = CrowdDis(rank(j))+(PopObj(rank(j+1),i)-PopObj(rank(j-1),i))/(Fmax(i)-Fmin(i));
        end
    end
end

function isture = Change(Global)
    isture= rem(Global.Current,Global.tao_n) == 0;
end
