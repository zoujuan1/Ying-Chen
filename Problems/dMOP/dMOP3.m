function varargout = dMOP3(Operation,Global,input)
% <problem> <AdMOP>
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

    A = unidrnd(Global.D);
    switch Operation
        case 'init'
            Global.M        = 2;
            Global.D        = 10;
            
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D); 
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            
            
            PopObj(:,1) = PopDec(:,A);
            
            t = (1/Global.NT)* floor(Global.Current()/Global.tao_n);
%             H = 0.75*sin(0.5*pi*t)+1.25;
            G = abs(sin(0.5*pi*t));
            if A ~= 1 && A ~= Global.D
                g = 1 + sum((PopDec(:,1:A-1)-G).^2,2)+sum((PopDec(:,A+1:end)-G).^2,2);
            elseif A ==1
                g = 1 + sum((PopDec(:,A+1:end)-G).^2,2);
            elseif A == Global.D
                g = 1+ sum((PopDec(:,1:A-1)-G).^2,2);
            end
            h = 1 - (PopObj(:,1)./g).^0.5;
            PopObj(:,2) = g.*h;
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f(:,1)    = (0:1/(input-1):1)';
            f(:,2)    = 1-f(:,1).^0.5;
            varargout = {f};
    end
end