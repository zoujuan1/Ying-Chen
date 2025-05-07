function [fname,n_var, n_obj, xl, xu] = terminate_problem(index)
    switch index
        case 1
            fname='MMF1';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[1 -1];     % the low bounds of the decision variables
            xu=[3 1];      % the up bounds of the decision variables
    end
end
