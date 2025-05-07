function NSGAII(Global)
% <algorithm> <H-N>
% Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population = Global.Initialization();
    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N);
    PP=cell(1,1000);
    ind=1;
    %% Optimization
    while Global.NotTermination(Population)
        P = PPF(10000,Global);
 
        if mod(Global.gen,20)==0
            ind=ind+1;
            PP{ind}=P;
        end
        for i = 1 : ind 
           Draw( PP{i},'k+','Markeredgecolor',[.4 .9 .9],'Markerfacecolor',[.9 .8 .9]);
        end
        Draw( PP{ind},'k+','Markeredgecolor',[.1 .4 .6],'Markerfacecolor',[.9 .8 .9]);
%         Draw(P,'k+','Markeredgecolor',[.1 .4 .6],'Markerfacecolor',[.9 .8 .9]);
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
    end
end

function P = PPF(input,Global)
            f(:,1)    = (0:1/(input-1):1)';
            f(:,2)    = 1-f(:,1).^(10*Global.gen/Global.maxgen) ;
            P = f;
end