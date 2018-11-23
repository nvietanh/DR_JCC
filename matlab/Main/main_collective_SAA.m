%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Energy and Reserve Dispatch with Distributionally Robust Joint Chance Constraints
%   Christos Ordoudis, Viet Anh Nguyen, Daniel Kuhn, Pierre Pinson
% 
%   This scripts generate the results for the collective optimization model
%   with rho = 0
%
%   In this case, we solve the SAA problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all; clc;

% call startup to add the necessary path
startup;

tic

%%

% Data input
RTS_Data2;
DRO_param.dual_norm = 'inf'; % dual norm
DRO_param.eps_joint_cvar = 0;

% Number of individual runs (number of coupled datasets in the numerical
% study)
IR_max = 100;
IR_sim = 100;

% Number of out of sample data for each individual run (N') for testing
% dataset
OOS_max = 1000;
OOS_sim = 1000;

% Number of maximum sample size
N_max = 200;

% Number of sample data in training dataset (N)
N_vect = [25, 50, 100, 200];

% Total number of data 
Nscen = IR_max * (N_max + OOS_max);

Load_All_Data;

% Initializing the matrices to gather final results
SAA_TC = NaN(length(N_vect), IR_sim);

for N_vect_idx = 1:length(N_vect)
    N = N_vect(N_vect_idx);
    for j = 1:IR_sim
        display(['out of sample iteration: ', num2str(j)]);
        

        % For each coupled dataset, we pick N and N' samples
        WPf_max = nWind(:,1:N_max,j)';
        WPr_max = nWind(:,N_max+1:N_max+OOS_max,j)';
        WPf = WPf_max(1:N,:);
        WPr = WPr_max(1:OOS_sim,:);

        % Build the corresponding data related to wind power production
        all = [1:N];
        system_info.Wscen = WPf(all,:)';
        system_info.mu = mean(system_info.Wscen,2); 
        system_info.xi = system_info.Wscen - repmat(system_info.mu, 1, size(system_info.Wscen,2));

        % Calculation of A,B,C,b matrices for joint chance constraints
        jcc = CC_matrices(system_info, DRO_param);     

        DRO_param.rho = 0;
        % Solve robust optimization model

        SAA_SP_sol = SAA_StoProg_solve(system_info);
        SAA_SP_y0_DA{j} = SAA_SP_sol.y0;
        SAA_SP_ru{j} = SAA_SP_sol.ru;
        SAA_SP_rd{j} = SAA_SP_sol.rd;
        TimeSAA(j) = toc;

        tic
        % Loop for each out-of-sample realization 
        for k = 1:OOS_sim
            system_info.Wreal = WPr(k,:)';
            system_info.DWreal = system_info.Wreal - system_info.mu;

            % Solve real-time optimal power flow
            RT_solution_SAA = RT_solve_R(system_info,SAA_SP_sol.y0,SAA_SP_sol.ru,SAA_SP_sol.rd);
            RT_Obj_IR(N_vect_idx,j,k) = RT_solution_SAA.Obj_RT;  
            RT_flag(j,k) = RT_solution_SAA.Flag;

        end

        TimeOOS(j) = toc;

        % Calculation of expected cost

        if SAA_SP_sol.Flag == 0
            SAA_TC(N_vect_idx, j) = system_info.cru'*SAA_SP_sol.ru + system_info.crd'*SAA_SP_sol.rd + mean(RT_Obj_IR(N_vect_idx,j,:));
        else
            SAA_TC(N_vect_idx, j) = NaN;
        end
        
    end % end of IR loop
end % end of N loop

toc
