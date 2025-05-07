function varargout = FDA5(Operation,Global,input)
% <problem> <FDA>
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
            Global.D        = 12;
            
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D); 
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            t = (1/Global.NT) * floor(Global.Current()/Global.tao_n);
            G = abs(sin(0.5*pi*t));
            F = 1 + 100*(sin(0.5*pi*t)).^4;
            g = G + sum((PopDec(:,3:end)-G).^2,2);
            Y = PopDec(:,1:2).^F;
            
            PopObj(:,1) = (1.+g).*cos(0.5*pi.*Y(:,2)).*cos(0.5*pi.*Y(:,1));
            PopObj(:,2) = (1.+g).*cos(0.5*pi.*Y(:,1)).*sin(0.5*pi.*PopDec(:,2));
            PopObj(:,3) = (1.+g).*sin(0.5*pi.*Y(:,1));
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f(:,1)    = (0:1/(input-1):1)';
            f(:,2)    = 1-f(:,1).^0.5;
            varargout = {f};
    end
end