function NUMSMOEAD(Global)
% <algorithm> <H-N>
% NUMS-MOEA/D: integration of preference in MOEA/D
% point --- --- referenece point
% propotion --- 0.2 --- the range of ROI
% kind --- 1 --- The type of aggregation function

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%% Parameter setting
[point,propotion,kind] = Global.ParameterSet(zeros(1, Global.M)+0.5,0.2,1);
t =0.1* floor(Global.gen/Global.tao_n);
point(1) = (sin(pi*t)).^2 + 0.1;
for i = 2 : Global.M
    point(i) = (cos(pi*t)).^2 + 0.1;
end

%% Generate the weight vectors
[W,Global.N] = NUMSPoint(Global.N,Global.M,point,propotion);
T = ceil(Global.N/10);

%% Detect the neighbours of each solution
B = pdist2(W,W);
[~,B] = sort(B,2);
B = B(:,1:T);

Population = Global.Initialization();
Z = min(Population.objs,[],1);

%% Optimization
while Global.NotTermination(Population)  
    % For each solution
    for i = 1 : Global.N
        % Choose the parents
        P = B(i,randperm(size(B,2)));
        
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
end
end