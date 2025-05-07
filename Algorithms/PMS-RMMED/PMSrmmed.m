function PMSrmmed(Global)
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
    metric1=[];
    Score1=[];
    eva=[];
    PopObj=[];
    PopPF=[];
    
    global hvd;
    hvd = 0.1;
    metric2 =[];
    T=[];
    %% Generate random population
%     Global.NT=10;
    Population = Global.Initialization();
%     [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N);
    %% Optimization
    k=0;         %��������������
    delta=2;    %��������10��
    %% Set time step t
    t=0;             %ʱ������
    OP=10;           %������������
    OB=10;            %����
    CenterPoint=[];  % ��ʷ���ĵ�
    MemoryPool=[];       %�����
    msize=100;    
    PF=Global.PF;
    while Global.NotTermination(Population)
%         Cc=INDIVIDUAL(sum(Population.decs,1)/Global.N);
%         DF(func2str(Global.problem),Global.NT);    %������ӦPF
%         Draw(Global.PF,'ok','Markeredgecolor',[.9 .0 .0],'Markerfacecolor',[.9 .0 .0],'MarkerSize',2);
%         Draw(Cc.objs,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',5);
%         Draw(Pt1.objs,'k+','Markeredgecolor',[.6 .0 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',3);
        if Global.evaluated>0.9*Global.evaluation     
             y=0;
        end
        if Change(Global)  % ��⻷���仯
            k=0;
            t=t+1 

            Decs=Population.decs;
            Population=INDIVIDUAL(Decs);
            
            Pt2=Population;
            [FrontNo,MaxFNo] = NDSort(Pt2.objs,Global.N);
            First = find(FrontNo==1);   
            Ct=sum(Pt2(First).decs,1)/length(First);  % �����֧������ĵ�
            CenterPoint=[CenterPoint;Ct];       % ��Ct�ӽ���ʷ���ĵ�
            while k<delta                %��������delta��
                
                
                Offspring  = Global.Variation(Pt2,Global.N,@RMMEDA_operator);
                [Pt2,FrontNo,CrowdDis] = EnvironmentalSelection([Pt2,Offspring],Global.N);
                k=k+1;
%                 Global.gen=Global.gen+1;
            end

            [FrontNo,MaxFNo] = NDSort(Pt2.objs,Global.N);
            First = find(FrontNo==1);
            St=sum(Pt2(First).decs,1)/length(First)-Ct ;              %������������
            
            P1=exploration(Pt2,St,Global,OB,OP);                       %��̽
            P2=exploiation(Population,t,CenterPoint,OP,Global);               %����
            P3=MemoryOp(Population,MemoryPool,msize);                        %����
            Pachive=[P1,P2,P3];
            for i = 1 : length(Pachive)              %����decs
                for j =1 : Global.D
                      temp=Pachive(i).decs;
                      if  temp(j)<Global.lower(j)
                             Decs1=0.5*(CenterPoint(t-1,j)+Global.lower(j));
                             Pachive(i)=INDIVIDUAL(Decs1);
                      end
                      if  temp(j)>Global.upper(j)
                             Decs1=0.5*(CenterPoint(t-1,j)+Global.upper(j));
                             Pachive(i)=INDIVIDUAL(Decs1);
                      end  
                end
            end
            
            [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Pachive,Population],Global.N);
            PF=Global.PF;
            if t>=2      %��ȥ��һ�λ����仯
                 score1=IGD(Population.objs,Global.PF);      %����IGD
                 metric1=[metric1,score1];
                 migd=sum(metric1)/length(metric1); %����MIGD
                 igd=migd;                         %��MIGDͨ��ȫ�ֱ���igd����ָ�꺯��MIGD()
                 T=[T,t];
                 Score1=[Score1,log10(score1)];
                 score2=HV(Population.objs,Global.PF);      %����HV
                 metric2=[metric2,score2];
                 mhvd=sum(metric2)/length(metric2); %����MGD
                 hvd=mhvd;                         %��Mhvdͨ��ȫ�ֱ���hvd����ָ�꺯��HV() 
                 temp=Population.objs+0.2*(t-1);
                 PopObj=[PopObj;temp];
                 PopPF=[PopPF;Global.PF+0.2*(t-1)];
            end
        else
            Offspring  = Global.Variation(Population,Global.N,@RMMEDA_operator);
           [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
        end
    end
end
function  P1 = exploration(Pop,St,Global,OB,OP)       %��̽����
    [CrowdDis,index] = sort(CrowdingDistance1(Pop.decs),'descend');   %���߿ռ�ӵ���Ⱦ���
    Ob=Pop(index(1:OP));           %�ҳ�POP��ӵ��������OP��������Ϊ��ȺOb
    Popb=Ob;                      %����ΪOB*OP�� ����
    Ob=Ob.decs;
    for k = 1:OB-1                %����OP*(OB-1)�� ������
        for j = 1:OP  
            Decs=[];
            for i = 1:Global.D
                if St(i)>0
                     temp=Ob(j,i)+((Global.upper(i)-Ob(j,i))*k)/OB;
                else
                     temp=Ob(j,i)-((Ob(j,i)- Global.lower(i))*k)/OB;
                end
                if  temp<Global.lower(i)
                    temp=0.5*(CenterPoint(t,i)+Global.lower(i));
                end
                if  temp>Global.upper(i)
                    temp=0.5*(CenterPoint(t,i)+Global.upper(i));
                end  
                Decs=[Decs,temp];
            end
            p=INDIVIDUAL(Decs);
            Popb=[Popb,p];
        end
    end
    [P1,FrontNo,CrowdDis]=EnvironmentalSelection(Popb,OP);
end
function P2 = exploiation(Pop,t,CenterPoint,OP,Global)                       %��������
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
function P3 = MemoryOp(Pop,MemoryPool,msize)                        %��������
    

    P=EnvironmentalSelection(Pop,msize);
    MemoryPool=[MemoryPool,P]; 
    P3=EnvironmentalSelection(MemoryPool,msize);                                %ѡ���������ĸ����������
 
    if length(MemoryPool)>2*length(Pop)
        MemoryPool=MemoryPool(msize:end)
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
