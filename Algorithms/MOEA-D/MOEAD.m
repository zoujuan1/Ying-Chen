function MOEAD(Global)
% <algorithm> <H-N>
% MOEA/D: A Multiobjective Evolutionary Algorithm Based on Decomposition
% kind --- 1 --- The type of aggregation function

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
    
    global hvd;
    hvd = 0.1;
    metric2 =[];
    %% Parameter setting
    kind = Global.ParameterSet(4);

    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    T = ceil(Global.N/10);
    t=0;
    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);
    
    %% Generate random population
    Population = Global.Initialization();
    Z = min(Population.objs,[],1);

    %% Optimization
    while Global.NotTermination(Population)
        
     
        if Change(Global)  % 检测环境变化
            t=t+1 
            Decs=Population.decs;
            Population=INDIVIDUAL(Decs);
        end
        % For each solution
        for i = 1 : Global.N      
            % Choose the parents
            P = B(i,randperm(size(B,2)));

            % Generate an offspring
            Offspring = Global.Variation(Population(P(1:2)),1);

            % Update the ideal point
            Z = min(Z,Offspring.obj);

            % Update the neighbours
            switch kind
                case 1
                    % PBI approach
                    normW   = sqrt(sum(W(P,:).^2,2));
                    normP   = sqrt(sum((Population(P).objs-repmat(Z,T,1)).^2,2));
                    normO   = sqrt(sum((Offspring.obj-Z).^2,2));
                    CosineP = sum((Population(P).objs-repmat(Z,T,1)).*W(P,:),2)./normW./normP;
                    CosineO = sum(repmat(Offspring.obj-Z,T,1).*W(P,:),2)./normW./normO;
                    g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                    g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
                case 2
                    % Tchebycheff approach
                    g_old = max(abs(Population(P).objs-repmat(Z,T,1)).*W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z),T,1).*W(P,:),[],2);
                case 3
                    % Tchebycheff approach with normalization
                    Zmax  = max(Population.objs,[],1);
                    g_old = max(abs(Population(P).objs-repmat(Z,T,1))./repmat(Zmax-Z,T,1).*W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),T,1).*W(P,:),[],2);
                case 4
                    % Modified Tchebycheff approach
                    g_old = max(abs(Population(P).objs-repmat(Z,T,1))./W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z),T,1)./W(P,:),[],2);
            end
            Population(P(g_old>=g_new)) = Offspring;
        end
        
        if t>=2      %除去第一次环境变化
             score1=IGD(Population.objs,Global.PF);      %计算IGD
             metric1=[metric1,score1];
             migd=sum(metric1)/length(metric1); %计算MIGD
             igd=migd;                         %将MIGD通过全局变量igd传给指标函数MIGD()

             score2=HV(Population.objs,Global.PF);      %计算HVD
             metric2=[metric2,score2];
             mhvd=sum(metric2)/length(metric2); %计算MHVD
             hvd=mhvd;                         %将MHVD通过全局变量hvd传给指标函数MHVD()
        end
    end
end

function isture = Change(Global)
    isture= rem(Global.Current,Global.tao_n) == 0;
end