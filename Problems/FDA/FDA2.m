function varargout = FDA2(Operation,Global,input)
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
            Global.D        = 21;
            
            Global.lower    = ones(1,Global.D);
            Global.lower    = - Global.lower; %È«ÖÃ-1
            Global.lower(1) = 0;
            
            Global.upper    = ones(1,Global.D); 
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            
            A = ceil(Global.D/2);
            
            PopObj(:,1) = PopDec(:,1);
            t = (1/Global.NT)*floor(Global.Current()/Global.tao_n);
            H = 0.75+0.7*sin(0.5*pi*t);
            
           
            g = 1+sum(PopDec(:,2:A).^2,2);
            T = H + sum(((PopDec(:,A+1:end)-H)).^2,2);
            h = 1-(PopObj(:,1)./g).^(1./T);
            
            PopObj(:,2) = g.*h;
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            t = (1/Global.NT)*floor(Global.Current()/Global.tao_n);
            H = 0.75+0.7*sin(0.5*pi*t);
            f(:,1)    = (0:1/(input-1):1)';
            f(:,2)    = 1-f(:,1).^((H+10*(1+H)*(1+H))^-1);
            varargout = {f};
    end
end