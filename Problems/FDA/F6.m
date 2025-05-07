function varargout = F6(Operation,Global,input)
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
            Global.M        = 2;
            Global.D        = 10;
            
            Global.lower    = zeros(1,Global.D);         
            Global.upper    = ones(1,Global.D).*5; 
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            y = [];
            PopDec = input;
            t = (1/Global.NT)* floor(Global.Current()/Global.tao_n);
            a = 2*cos(1.5*pi*t)*sin(0.5*pi*t)+2;
            b = 2*cos(1.5*pi*t)*cos(0.5*pi*t)+2;
            H = 1.25 + 0.75*sin(pi*t);
            for i = Global.D -2 : Global.D
                y(:,i) = PopDec(:,i) - b - 1 + (abs(PopDec(:,1)-a)).^(H+i/Global.D);
            end

            PopObj(:,1) = (abs(PopDec(:,1) - a)).^H + sum(y(:,1:2:end).^2,2);
            PopObj(:,2) = (abs(PopDec(:,1) - a - 1)).^H + sum(y(:,2:2:end).^2,2);
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            t =(1/Global.NT) * floor(Global.Current()/Global.tao_n);
            H = 1.25 + 0.75*sin(pi*t);
            s = 0:1/(input-1):1;
            f(:,1)    = s'.^H;
            f(:,2)    = (1-s).^H;
            varargout = {f};
    end
end