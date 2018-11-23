function[sol] = SAA_StoProg_solve(si)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Energy and Reserve Dispatch with Distributionally Robust Joint Chance Constraints
%   Christos Ordoudis, Viet Anh Nguyen, Daniel Kuhn, Pierre Pinson
%   
%   This function the SAA problem where each training sample has its own
%   corresponding set of wind spill, load shedding and
%   reserve adjustments as variables.
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
    Nnodes = size(si.AG,1);
    Ndemand = length(si.D);
   

    % Definition of variables
    y1 = sdpvar(Nunits, 1); % Day-ahead power production from thermal power plants
    ru = sdpvar(Nunits, 1); % Day-ahead power production from thermal power plants
    rd = sdpvar(Nunits, 1); % Day-ahead power production from thermal power plants
    Y = sdpvar(Nunits, Nscen, 'full'); % real-time power production
    Y2 = sdpvar(Nunits, Nscen, 'full'); % reserve increase
    Y3 = sdpvar(Nwind , Nscen, 'full'); % wind spillage
    Y4 = sdpvar(Ndemand, Nscen, 'full'); % loadshedding
    fi_real = sdpvar(Nnodes, Nscen, 'full'); % for node injection
 
   
    % Constraints set
    CS = [];
    
    % Day-ahead constraints    
    CS = [CS, si.Pmin <= y1 - rd, y1 + ru <= si.Pmax, 0 <= ru <= si.ResCap, 0 <= rd <= si.ResCap];
    CS = [CS, sum(y1) + sum(si.Wmax.*si.mu) - sum(si.D) == 0];                              % day ahead
    CS = [CS, sum(Y, 1) + sum(si.DiagWmax*si.xi, 1) - sum(Y3, 1) + sum(Y4, 1) == 0];        % for each scenario

    
    % First constraint about generation
    CS = [CS, repmat(rd, 1, Nscen) <= Y <= repmat(ru, 1, Nscen)];
    % Second constraint about line flow
    CS = [CS, sum(fi_real, 1) == 0, repmat(-si.F, 1, Nscen) <= si.PTDF*fi_real <= repmat(si.F, 1, Nscen)];
    CS = [CS, si.AG*(repmat(y1, 1, Nscen) + Y) + si.AW * si.DiagWmax*(si.mu+si.xi) - si.AD*(repmat(si.D, 1, Nscen) - Y4) - si.AW * Y3 == fi_real];
    %CS = [CS, 0 <= si.AW * wsp <= si.AW * si.DiagWmax*si.Wreal];

    % Add the last constraints
    CS = [CS, 0 <= Y2 <= repmat(y1-rd, 1, Nscen) ,0 <= Y3 <= si.DiagWmax*(repmat(si.mu, 1, Nscen) + si.xi), 0 <= Y4 <= repmat(si.D, 1, Nscen)];
    
    Obj = si.cru'*ru + si.crd'*rd + si.c'*y1 + (sum(si.c'*Y) + si.cl*sum(Y4(:)) + si.cw*sum(Y3(:)) + si.cr*sum(Y2(:)))/Nscen;
    % Settings
    optim_options = sdpsettings('solver', 'gurobi','gurobi.TimeLimit',10000,'gurobi.NumericFocus',3,'verbose',0);

    % Solve
    sol = optimize(CS, Obj, optim_options);

    sol.y0 = value(y1);
    sol.Y = value(Y);
    sol.ru = value(ru);
    sol.rd = value(rd);
    sol.Y = value(Y);
    sol.Obj = value(Obj);
    sol.Flag = sol.problem;
    
end