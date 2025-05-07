function varargout = F8(Operation,Global,input)
% <problem> <F>
% Comparison of Dynamic Multiobjective Evolutionary Algorithms: Empirical Results
% operator --- EAreal

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    switch Operation
        case 'init'
            Global.M        = 3;
            Global.D        = 10;
            
            Global.lower    = ones(1,Global.D).*(-1);         
            Global.upper    = ones(1,Global.D).*2; 
            Global.lower(1) = 0;
            Global.upper(1) = 1;
            Global.lower(2) = 0;
            Global.upper(2) = 1;
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            y = [];
            PopDec = input;
            t = (1/Global.NT)* floor(Global.Current()/Global.tao_n);
            G = sin(0.5*pi*t);
            H = 1.25 + 0.75*sin(pi*t);
            g = sum((PopDec(:,3:end) - ((repmat(PopDec(:,1),1,Global.D-2)+ repmat(PopDec(:,2),1,Global.D-2))./2).^H - G).^2,2);
            PopObj(:,1) = (1 + g).*cos(0.5*pi.*PopDec(:,2)).*cos(0.5*pi.*PopDec(:,1));
            PopObj(:,2) = (1 + g).*cos(0.5*pi.*PopDec(:,2)).*sin(0.5*pi.*PopDec(:,1));
            PopObj(:,3) = (1 + g).*sin(0.5*pi.*PopDec(:,2));
            
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            % Generate a set of reference points on the true Pareto front
            f = UniformPoint(input,Global.M);
            f = f./repmat(sqrt(sum(f.^2,2)),1,Global.M);
            % Return the reference points
            varargout = {f};
    end
end