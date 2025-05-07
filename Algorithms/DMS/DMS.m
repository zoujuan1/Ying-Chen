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
    t=0;             %ʱ������
    
    m=10;            %ÿ������ Gra��������
    OP=Global.N/m;           %����
    
    r1=0.5;
    r2=0.5;
    CenterPoint=[];  % ��ʷ���ĵ�
    Zlow=[];         % ��ʷ��С��
    Zhigh=[];        % ��ʷ����� 
    PF=Global.PF;
    while Global.NotTermination(Population)
%         Cc=INDIVIDUAL(sum(Population.decs,1)/Global.N);
%         DF(func2str(Global.problem),Global.NT);    %������ӦPF
%         Draw(Global.PF,'ok','Markeredgecolor',[.9 .0 .0],'Markerfacecolor',[.9 .0 .0],'MarkerSize',2);
%         Draw(Cc.objs,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',5);
%         Draw(Pt1.objs,'k+','Markeredgecolor',[.6 .0 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',3);
        
        if Change(Global)  % ��⻷���仯
            t=t+1 
            Decs=Population.decs;
            Population=INDIVIDUAL(Decs);
            
            Pt2=Population;
            [FrontNo,MaxFNo] = NDSort(Pt2.objs,Global.N);
            First = find(FrontNo==1);   
            Ct=sum(Pt2(First).decs,1)/length(First);  % �����֧������ĵ�
            CenterPoint=[CenterPoint;Ct];       % ��Ct�ӽ���ʷ���ĵ�
            low = min(Population.decs,[],1);        %��С��
            high  = max(Population.decs,[],1);   %�����
            Zlow =[Zlow;low];         
            Zhigh =[Zhigh;high];
            
            P1=GraSearch(Population,t,Zlow,Zhigh,m,OP,Global,CenterPoint);                       %Gra Search
            P2=prediction(Population,t,CenterPoint,r1*Global.N,Global);               %Ԥ��
            P3=DM(Zlow,Zhigh,r2*Global.N,t,CenterPoint,Global);                        %  random diversity maintenance strategy
            Pachive=[P1,P2,P3];

            
            [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Pachive,Population],Global.N);
            PF=Global.PF;
            if t>=2      %��ȥ��һ�λ����仯
                 score=IGD(Population.objs,Global.PF);      %����IGD
                 metric=[metric,score];
                 migd=sum(metric)/length(metric); %����MIGD
                 igd=migd;                         %��MIGDͨ��ȫ�ֱ���igd����ָ�꺯��MIGD()
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
        lowD=Zlow(t,:)-Zlow(t-1,:);   %��С���ƶ�����
        highD=Zhigh(t,:)-Zhigh(t-1,:); %������ƶ�����
    end
    if  norm(lowD) >  norm(highD)
        D=lowD;
    else
        D=highD;
    end
    
    [CrowdDis,index] = sort(CrowdingDistance1(Pop.decs),'descend');   %���߿ռ�ӵ���Ⱦ���
    Ob=Pop(index(1:m));           %�ҳ�POP��ӵ��������m��������Ϊ��ȺOb
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
function P2 = prediction(Pop,t,CenterPoint,OP,Global)          %Ԥ��              
     [CrowdDis,index] = sort(CrowdingDistance1(Pop.decs),'descend');  %���߿ռ�ӵ���Ⱦ���
     Ob=Pop(index(1:OP));                                             %�ҳ�POP��ӵ��������OP��������Ϊ��ȺO
     if t>=2
        Decs=Ob.decs+(CenterPoint(t,:)-CenterPoint(t-1,:));               %Ԥ��
        for i = 1 : OP              %����decs
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
        lowD=Zlow(t,:)-Zlow(t-1,:);   %��С���ƶ�����
        highD=Zhigh(t,:)-Zhigh(t-1,:); %������ƶ�����
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
