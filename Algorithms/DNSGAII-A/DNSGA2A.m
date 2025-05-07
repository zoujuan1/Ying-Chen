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
    t=0;             %ʱ������
    PF=Global.PF;
    while Global.NotTermination(Population)
        Cc=INDIVIDUAL(sum(Population.decs,1)/Global.N);
        DF(func2str(Global.problem),Global.NT);    %������ӦPF
        Draw(Global.PF,'ok','Markeredgecolor',[.9 .0 .0],'Markerfacecolor',[.9 .0 .0],'MarkerSize',2);
        Draw(Cc.objs,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',5);
% %         Draw(Pt1.objs,'k+','Markeredgecolor',[.6 .0 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',3);
        if Global.evaluated>0.9*Global.evaluation     
             y=0;
        end
        if Change(Global)   % ��⻷���仯
            t=t+1   
            num=floor(Global.N*0.2);   %�ðٷ�֮20�ĸ������
            Pt2=Population;
           
            
            Pp = Global.Initialization();
            Pop=[Pp(1:num),Pt2(num:Global.N)];
            Decs=Pop.decs;
            Population=INDIVIDUAL(Decs);
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
