function HDMS(Global)
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
    k=0;         %自主进化计数器
    delta=2;    %自主进化10代
    %%                                                                                                                                                                                                                                                                                                             Set time step t
    t=0;             %时间序列
    OP=10;           %引导个体数量
    OB=10;            %划分
    CenterPoint=[];  % 历史中心点
    PF=Global.PF;

    while Global.NotTermination(Population)
%         Cc=INDIVIDUAL(sum(Population.decs,1)/Global.N);
%         DF(func2str(Global.problem),Global.NT);    %画出对应PF
%         Draw(Global.PF,'ok','Markeredgecolor',[.9 .0 .0],'Markerfacecolor',[.9 .0 .0],'MarkerSize',2);
%         Draw(Cc.objs,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if mod(Global.gen,1)
             S = IGD(Population.objs,Global.PF);
             Score1=[Score1,S];
             eva=[eva,Global.evaluated];
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if t==12
%             Pop=Population.objs;
%             y=1;
%         end
        
        if Global.evaluated>0.9*Global.evaluation
             y=0;
        end
        if Change(Global)  % 检测环境变化
            k=0;
            t=t+1 
            Decs=Population.decs;
            Population=INDIVIDUAL(Decs);
            
            Pt2=Population;
            [FrontNo,MaxFNo] = NDSort(Pt2.objs,Global.N);
            First = find(FrontNo==1);   
            Ct=sum(Pt2(First).decs,1)/length(First);  % 计算非支配层中心点
            CenterPoint=[CenterPoint;Ct];       % 将Ct加进历史中心点
            

            
            while k<delta                %自主进化delta代
                
%                  MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
%                  Offspring  = Global.Variation(Pt2(MatingPool));
                  Offspring  = Global.Variation(Pt2,Global.N,@RMMEDA_operator);
                [Pt2,FrontNo,CrowdDis] = EnvironmentalSelection([Pt2,Offspring],Global.N);
                k=k+1;
%                 Global.gen=Global.gen+1;
            end

            [FrontNo,MaxFNo] = NDSort(Pt2.objs,Global.N);
            First = find(FrontNo==1);
            St=sum(Pt2(First).decs,1)/length(First)-Ct;               %自主进化方向
            
            P1=exploration(Pt2,St,Global,OB,OP,CenterPoint,t);                %自主进化
            P2=exploiation(Population,t,CenterPoint,Global);               %中心点预测
              
            P3=PCM(Population,100,Global);                                       %变异
            Pachive=[P1,P2,P3];
            
            [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Pachive,Population],Global.N);
            PF=Global.PF;
            if t>=2      %除去第一次环境变化
                 score1=IGD(P2.objs,Global.PF);      %计算IGD
                 metric1=[metric1,score1];
                 migd=sum(metric1)/length(metric1); %计算MIGD
                 igd=migd;                         %将MIGD通过全局变量igd传给指标函数MIGD()
                 T=[T,t];
                 Score1=[Score1,log10(score1)];
                 score2=HV(Population.objs,Global.PF);      %计算HV
                 metric2=[metric2,score2];
                 mhvd=sum(metric2)/length(metric2); %计算MGD
                 hvd=mhvd;                          %将Mhvd通过全局变量hvd传给指标函数HV() 
                 
%                  if t==5||t==10||t==20||t==25||t==30
%                     PopObj=[PopObj;Population.objs];
%                     PopPF=[PopPF;Global.PF];
%                  end
            end

        else
%             MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
%             Offspring  = Global.Variation(Population(MatingPool));
             Offspring  = Global.Variation(Population,Global.N,@RMMEDA_operator);
            [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
           
        end
    end
end
function  P1 = exploration(Pop,St,Global,OB,OP,CenterPoint,t)       %开探算子
    [CrowdDis,index] = sort(CrowdingDistance1(Pop.decs),'descend');   %决策空间拥挤度距离
    Ob=Pop(index(1:OP));           %找出POP中拥挤度最大的OP个个体作为种群Ob
    Popb=Ob;                      %最终为OB*OP个 个体
    Ob=Ob.decs;  
    low = min(Pop.decs,[],1);        %极小点
    for k = 1:OB-1                %产生OP*(OB-1)个 检测个体
        for j = 1:OP  
            Decs=[];
            for i = 1:Global.D
                if St(i)>0
                     temp=Ob(j,i)+((Global.upper(i)-low(i))*k)/OB;
                else
                     temp=Ob(j,i)-((low(i)- Global.lower(i))*k)/OB;
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
    P1=Popb;
end
function P2 = exploiation(Pop,t,CenterPoint,Global)                       %开采算子
     [CrowdDis,index] = sort(CrowdingDistance1(Pop.decs),'descend');  %决策空间拥挤度距离
     Ob=Pop(index(1:Global.N));                                             %找出POP中拥挤度最大的OP个个体作为初始种群
     flag=1;
     if t>=2
        if  flag ==1                                            %如果环境变化相似不相同
            Decs=[];
            for i = 1 : Global.N
                temp=Ob(i).decs+(0.5+rand(1,Global.D)).*(CenterPoint(t,:)-CenterPoint(t-1,:));          %学习率为0.5-1.5 随机预测
                Decs=[Decs;temp];
            end
%             Decs=Ob.decs+(0.5+rand(1,Global.D)).*(CenterPoint(t,:)-CenterPoint(t-1,:));              
            
        else
            Decs=Ob.decs+(CenterPoint(t,:)-CenterPoint(t-1,:));                                       %学习率为1 随机预测
        end
        
        for i = 1 : Global.N              %修正decs
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
function P3 = PCM(Pop,OP,Global)                        %变异算子
     [CrowdDis,index] = sort(CrowdingDistance1(Pop.decs),'descend');  %决策空间拥挤度距离
     Ob=Pop(index(1:OP));                                             %找出POP中拥挤度最大的OP个个体作为种群
     Ob=Ob.decs;
     Decs=[];
     for i = 1 : OP
        Temp=[];
        for j = 1 : Global.D
%           newPop(i,j) = Ob(i,j);
            r=unidrnd(6);
            delta_alpha = 1/(10.^unidrnd(2))*unidrnd(9);    %delta_alpha的大小
            delta_beta = Ob(i,j) * delta_alpha - Ob(i,j);    %delta_beta的大小
            delta_gama = Ob(i,j) / delta_alpha - Ob(i,j);    %delta_gama的大小
            if r==1
                xtemp = Ob(i,j) + delta_alpha;
            elseif r==2
                    xtemp = Ob(i,j) - delta_alpha;
            elseif r==3
                    xtemp = Ob(i,j) + delta_beta;
            elseif r==4
                    xtemp = Ob(i,j) - delta_beta;
            elseif r==5
                    xtemp = Ob(i,j) + delta_gama;
            else
                 xtemp = Ob(i,j) - delta_gama;
            end
            if  xtemp < Global.lower(j)
                    xtemp = 0.5*(xtemp + Global.lower(j));
            end
            if  xtemp > Global.upper(j)
                    xtemp = 0.5*(xtemp + Global.upper(j));
            end  
            Temp=[Temp,xtemp];
        end 
        Decs=[Decs;Temp];
     end
     
      P3=INDIVIDUAL(Decs);
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
