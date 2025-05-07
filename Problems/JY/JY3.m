function varargout = JY3(Operation,Global,input)
% <problem> <A>
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
            y=PopDec;
            t =(1/Global.NT)* floor(Global.Current()/Global.tao_n);
            A=0.05;
            W=floor(6*sin(0.5*pi*(t-1)));
            a=floor(100*sin(0.5*pi*t)^2);
            y(:,1)=abs(PopDec(:,1).*sin((2*a+0.5)*pi.*PopDec(:,1)));
            
            g = 1+sum((y(:,2:end).^2-y(:,1:end-1)).^2,2);
            PopObj(:,1) =g.*( y(:,1) + A*sin(W*pi.*y(:,1)));
            PopObj(:,2) = g.*(1-y(:,1)+A*sin(W*pi.*y(:,1)));
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            t =(1/Global.NT)* floor(Global.Current()/Global.tao_n);
            W=floor(6*sin(0.5*pi*(t-1)));
            x=(0:1/(input-1):1)';
            f(:,1)=x+0.05*sin(W*pi*x);
            f(:,2)=1-x+0.05*sin(W*pi*x);
            varargout = {f};
    end
end