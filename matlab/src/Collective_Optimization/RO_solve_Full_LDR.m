function[sol] = RO_solve_Full_LDR(si,DRO_param,jcc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Energy and Reserve Dispatch with Distributionally Robust Joint Chance Constraints
%   Christos Ordoudis, Viet Anh Nguyen, Daniel Kuhn, Pierre Pinson
%   
%   This function solves the full Linear Decision Rule problem where the
%   load shedding, wind spill and reserve increase are also modelled as LDR
%   
%   This is part of the collective optimization model.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    yalmip('clear')

    % Getting the number of thermals power plants, wind farms, scenarions,
    % transmission lines and nodes
    Nunits = length(si.Pmax);
    Nwind = length(si.Wmax);
    Nscen = size(si.Wscen,2);
    Nvertices_RO = size(si.Wscen_RO,2);
    Ndemand = length(si.D);
   

    % Definition of variables
    y0 = sdpvar(Nunits, 1); % Day-ahead power production from thermal power plants
    ru = sdpvar(Nunits, 1); % Day-ahead power production from thermal power plants
    rd = sdpvar(Nunits, 1); % Day-ahead power production from thermal power plants
    Y = sdpvar(Nunits, Nwind, 'full'); % Linear decision rule for real-time power production
    Y2 = sdpvar(Nunits, Nwind, 'full'); % Linear decision rule for reserve increase
    Y3 = sdpvar(Nwind , Nwind, 'full'); % Linear decision rule for wind spillage
    Y4 = sdpvar(Ndemand, Nwind, 'full'); % Linear decision rule for loadshedding
    
    % variables for the objective function
    z = sdpvar(Nwind, 1);
    s_obj = sdpvar(1, Nscen); % sigma variable for obj
    lambda_obj = sdpvar(1, 1); % lambda variable for obj
    
    % create x by stacking up y0, ru and rd
    x = [y0; ru; rd];
    
    % Constraints set
    CS = [];
    
    % Day-ahead constraints    
    CS = [CS, si.Pmin <= y0 - rd, y0 + ru <= si.Pmax, 0 <= ru <= si.ResCap, 0 <= rd <= si.ResCap];
    CS = [CS, sum(y0) + sum(si.Wmax.*si.mu) - sum(si.D) == 0];                              % day ahead
    CS = [CS, sum(Y, 1) + si.Wmax' - ones(Nwind,1)'*Y3 + ones(Ndemand,1)'*Y4 == 0];         % real time

    % Run a for-loop to add the constraints related to the robust
    % constraints
    % The set of code below is generic, it can be copied and paste for any
    % structure jcc of interest
    
    % find the number of Individual chance constraints we have
    nICC = 0;
    for j=1:size(jcc, 1)
        nICC = nICC + size(jcc{j, 1}, 1);
    end
    
    for j=1:size(jcc, 1)
        A_C{j,1} = jcc{j,1};
        B_C{j,1} = jcc{j,2};
        C_C{j,1} = jcc{j,3};
        b_C{j,1} = jcc{j,4};
        % the fifth element is for epsilon
        B_C2{j, 1} = jcc{j, 6};
        B_C3{j, 1} = jcc{j, 7};
        B_C4{j, 1} = jcc{j, 8};
    end
    A = cell2mat(A_C);
    B = cell2mat(B_C);
    C = cell2mat(C_C);
    b = cell2mat(b_C);
    
    B2 = cell2mat(B_C2);
    B3 = cell2mat(B_C3);
    B4 = cell2mat(B_C4);
    
    for j = 1:nICC
        CS = [CS, repmat(A(j,:)*x - b(j), 1, Nvertices_RO) + (B(j,:)*Y + B2(j,:)*Y2 + B3(j,:)*Y3 + B4(j,:)*Y4+C(j,:))*si.zi  <= 0 ];
    end
    
    % Add the last constraints
    CS = [CS, 0 <= Y2*si.zi, 0 <= Y3*si.zi <= si.DiagWmax*repmat(si.mu, 1, size(si.zi,2)), 0 <= Y4*si.zi <= repmat(si.D, 1, size(si.zi,2))];
    
    % Build the objective function 
    Obj = si.cru'*ru + si.crd'*rd + si.c'*y0 + DRO_param.rho*lambda_obj + 1/Nscen * sum(s_obj); 
    CS = [CS, z'*si.xi <= s_obj];
    CS = [CS, norm( z, DRO_param.dual_norm) <= lambda_obj];
    CS = [CS, z == Y'*si.c + si.cl*Y4'*ones(Ndemand, 1) + si.cw*Y3'*ones(Nwind, 1) + si.cr*Y2'*ones(Nunits, 1)];
    % Settings
    optim_options = sdpsettings('solver', 'gurobi','gurobi.TimeLimit',1000,'gurobi.NumericFocus',3,'verbose',0);

    % Solve
    sol = optimize(CS, Obj, optim_options);

    sol.y0 = value(y0);
    sol.Y = value(Y);
    sol.ru = value(ru);
    sol.rd = value(rd);
    sol.Y = value(Y);
    sol.fy = si.Qg*value(y0) + si.Qw*si.DiagWmax * si.mu - si.Qd*si.D;
    sol.fY = si.Qg*value(Y) + si.Qw*si.DiagWmax;
    sol.q = si.PG * value(y0);
    sol.qY = si.PG * value(Y);
    sol.Obj = value(Obj);
    sol.Flag = sol.problem;
    
end