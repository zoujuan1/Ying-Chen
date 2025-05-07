function MO_Ring_PSO_SCD(Global)
% <algorithm> <EMMO>
    global igd;
    igd=0.1;
    metric1=[];
    T=[];
    Score1=[];
    global hvd;
    hvd = 0.1;
    metric2 =[];
    %% Initialize parameters
    n_PBA = 5;          % Maximum size of PBA
    n_NBA = 3*n_PBA;	% Maximum size of NBA
    %% Generate random population
    mv   = 0.5*(Global.upper-Global.lower);
    Vmin = -mv;
    Vmax = mv;
    ParticleDec = Global.lower+(Global.upper-Global.lower).*rand(Global.N,Global.D);
    ParticleVel = Vmin+2.*Vmax.*rand(Global.N,Global.D);
    Population  = INDIVIDUAL(ParticleDec,ParticleVel);

    %% Initialize personal best archive PBA and Neighborhood best archive NBA
    PBA = cell(1,Global.N);
    NBA = cell(1,Global.N);
    for i = 1:Global.N
        PBA{i} = Population(i);
        NBA{i} = Population(i);
    end
    t=0;
    %% Optimization
    while Global.NotTermination(Population)
        if Change(Global)
            t=t+1
            Population=Global.Initialization();
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
        NBA = UpdateNBA(NBA,n_NBA,PBA);
        Population = Operator(Population,PBA,NBA,Global);
        PBA = UpdatePBA(Population,PBA,n_PBA);
        if Global.gen >= Global.maxgen
            tempNBA = [];
            for i = 1:Global.N
                tempNBA = [tempNBA,NBA{i}];
            end
            [tempNBA,FrontNo,~] = non_domination_scd_sort(tempNBA,Global.N);
            Global.NotTermination(tempNBA(FrontNo==1));
        end
        end
    end
end