function DNSGA2A(Global)
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
    T=[];
    Score1=[];
    global hvd;
    hvd = 0.1;
    metric2 =[];
    %% Generate random population
%     Global.NT=10;
    Population = Global.Initialization();
    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N);
    %% Optimization
    
    %% Set time step t
    t=0;             %时间序列
    PF=Global.PF;
    while Global.NotTermination(Population)
        Cc=INDIVIDUAL(sum(Population.decs,1)/Global.N);
        DF(func2str(Global.problem),Global.NT);    %画出对应PF
        Draw(Global.PF,'ok','Markeredgecolor',[.9 .0 .0],'Markerfacecolor',[.9 .0 .0],'MarkerSize',2);
        Draw(Cc.objs,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',5);
% %         Draw(Pt1.objs,'k+','Markeredgecolor',[.6 .0 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',3);
        if Global.evaluated>0.9*Global.evaluation     
             y=0;
        end
        if Change(Global)   % 检测环境变化
            t=t+1   
            num=floor(Global.N*0.2);   %用百分之20的个体随机
            Pt2=Population;
           
            
            Pp = Global.Initialization();
            Pop=[Pp(1:num),Pt2(num:Global.N)];
            Decs=Pop.decs;
            Population=INDIVIDUAL(Decs);
            PF=Global.PF;
            if t>=2      %除去第一次环境变化
                 score1=IGD(Population.objs,Global.PF);      %计算IGD
                 metric1=[metric1,score1];
                 migd=sum(metric1)/length(metric1); %计算MIGD
                 igd=migd;                         %将MIGD通过全局变量igd传给指标函数MIGD()
                 T=[T,t];
                 Score1=[Score1,log10(score1)];
                 
                 score2=HV(Population.objs,Global.PF);      %计算HV
                 metric2=[metric2,score2];
                 mhvd=sum(metric2)/length(metric2); %计算MGD
                 hvd=mhvd;                         %将Mhvd通过全局变量hvd传给指标函数HV() 
            end
        else
            MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
            Offspring  = Global.Variation(Population(MatingPool));
           [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
        end
    end
end

function isture = Change(Global)
    isture= rem(Global.Current,Global.tao_n) == 0;
end
