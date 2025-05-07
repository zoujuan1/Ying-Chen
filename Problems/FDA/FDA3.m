function varargout = FDA3(Operation,Global,input)
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
            Global.M        = 2;
            Global.D        = 20;
            
            Global.lower    = ones(1,Global.D);
            Global.lower    = - Global.lower; %È«ÖÃ-1
            Global.lower(1) = 0;
            
            Global.upper    = ones(1,Global.D); 
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            
            A = 1;
            
            t = (1/Global.NT)* floor(Global.Current()/Global.tao_n);
            G = abs(sin(0.5*pi*t));
            F = 10.^(2.*sin(0.5*pi*t));
            PopObj(:,1) = sum(PopDec(:,1:A).^F,2); 
            g = 1 + G + sum((PopDec(:,A+1:end)-G).^2,2);
            h = 1-(PopObj(:,1)./g).^0.5;
            PopObj(:,2) = g.*h;
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            t = (1/Global.NT)* floor(Global.Current()/Global.tao_n);
            F = 10^(2*sin(0.5*pi*t));
            G = abs(sin(0.5*pi*t));
            f(:,1)    = (0:1/(input-1):1)'.^F;
            f(:,2)    = (1-sqrt(f(:,1)./(1+G)))*(1+G);
            
  
            varargout = {f};
%             varargout = {f1,f};
    end
end