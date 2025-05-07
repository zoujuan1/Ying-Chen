function [W,N] = NUMSPoint(N,M,point,propotion)
%UniformPoint - Generate a set of uniformly distributed points on the unit
%hyperplane
%
%   [W,N] = UniformPoint(N,M) returns approximate N uniformly distributed
%   points with M objectives.
%
%   Example:
%       [W,N] = UniformPoint(275,10)

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

m_type = 1;    % 1: linear model; 2: exponential model
b_flag = 2;    %1：保留边界点； 2：去除边界点
ref_point = point;

vote_w = ref_point ./ sum(ref_point);

[W, N] = UniformPoint(N, M);

N = size(W, 1);
transformed_w = zeros(N, M);

if m_type == 1
    for i = 1 : N
        if sum(abs(vote_w - W(i, :))) == 0
            transformed_w(i, :) = W(i, :);
        else
            t_i = vote_w * norm(vote_w - W(i, :)) ./ (vote_w - W(i, :));
            idx = find(t_i > 0);
            [min_value, min_idx] = min(t_i(idx));
            temp = min_value - norm(vote_w - W(i, :));
            if(temp < 1.0e-04) && (b_flag ==1)
                t_d = norm(vote_w - W(i, :));
            elseif (temp < 1.0e-04) && (b_flag == 2)
                t_d = norm(vote_w - W(i, :)) * propotion;
            else
                t_d = norm(vote_w - W(i, :)) * propotion;
            end
            transformed_w(i, :) = vote_w + t_d * (W(i, :) - vote_w) / norm(vote_w - W(i, :));
        end
    end
else
    %% calculate the 'eta' parameter
    eta = eta_calculattion(M, H, propotion, b_flag);
    for i = 1: N
        if sum(abs(vote_w - W(i, :))) < 1.0e-04
            transformed_w(i, :) = W(i, :);
        else
            t_i = vote_w *norm(vote_w - W(i, :)) ./ (vote_w - W(i, :));
            idx = find(t_i > 0);
            [min_value, min_idx] = min(t_i(idx));
            temp = min_value - norm(vote_w - W(i, :));
            if (temp < 1.0e-04) && (b_flag == 1)
                temp = 0;
                t_d = min_value - min_value * (temp / min_value)^(1 / (eta + 1));
                bound(i, :) = W(i, :);
            elseif (temp < 1.0e-04) && (b_flag == 2)
                t_d = norm(vote_w - W(i, :)) * propotion;
            else
                t_d = min_value - min_value * (temp / min_value)^(1 / (eta + 1));
            end
            transformed_w(i, :) = vote_w + t_d * (W(i, :) - vote_w) / norm(vote_w - W(i, :));
        end
    end
end

W = [transformed_w; vote_w];
N = size(W, 1);
end