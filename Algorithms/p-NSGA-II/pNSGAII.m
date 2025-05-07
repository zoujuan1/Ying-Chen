function  pNSGAII(Global)
% <algorithm> <P>
% Nondominated sorting genetic algorithm II
% Point --- --- Preferred point
% delta  --- 0.1 --- Non-r-dominance threshold


%------------------------------- Reference --------------------------------
% K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, A fast and elitist
% multiobjective genetic algorithm: NSGA-II, IEEE Transactions on
% Evolutionary Computation, 2002, 6(2): 182-197.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [Point,delta] = Global.ParameterSet(zeros(1,Global.M)+0.5, 0.1);
    dpj = Inf;
    dpj_old = Inf;
    %% Generate random population
    Population = Global.Initialization();
    [~] = EnvironmentalSelectionA(Population,Global.M,Global.N, Point, delta,Global.gen,dpj_old,dpj);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,sum(max(0,Population.cons),2));
        Offspring  = Global.Variation(Population(MatingPool));
        [Population] = EnvironmentalSelectionA([Population,Offspring],Global.M,Global.N, Point, delta,Global.gen,dpj_old,dpj);
    end
end