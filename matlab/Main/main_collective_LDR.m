%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Energy and Reserve Dispatch with Distributionally Robust Joint Chance Constraints
%   Christos Ordoudis, Viet Anh Nguyen, Daniel Kuhn, Pierre Pinson
% 
%   This scripts generate the results for the collective optimization model
%   with rho > 0
%
%   In this case, we solve the problem with the robust constraints and
%   worst-case expectation objective function
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
rho_vectorC = [0 linspace(0.0001, 0.04, 24)];
% Initializing the matrices to gather final results
RO_TC = NaN(length(N_vect), length(rho_vectorC), IR_sim);

for N_vect_idx = 4:length(N_vect)
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

        % Build the corresponding data for RO
        system_info.Wscen_RO = de2bi(0:2^size(system_info.Wmax,1)-1)';
        system_info.Wexp_RO = mean(system_info.Wscen_RO,2); 
        system_info.zi = system_info.Wscen_RO - repmat(system_info.mu, 1, size(system_info.Wscen_RO,2));

        % Calculation of A,B,C,b matrices for joint chance constraints
        jcc = CC_matrices(system_info, DRO_param);
        

        tic
        % Loop for each value of \rho in P vector
        for i = 1:length(rho_vectorC) 
            DRO_param.rho = rho_vectorC(i);
            % Solve robust optimization model

            RO_sol = RO_solve_Full_LDR(system_info, DRO_param, jcc);
            RO_y0_DA{j,i} = RO_sol.y0;
            RO_ru{j,i} = RO_sol.ru;
            RO_rd{j,i} = RO_sol.rd;
            RO_obj{j,i} = RO_sol.Obj;
            RO_flag{j,i} = RO_sol.Flag;
            RO_Y{j,i} = RO_sol.Y * system_info.zi;
            RO_Qy{j,i} = RO_sol.q;
            RO_QY{j,i} = RO_sol.qY * system_info.zi;
            RO_Fy{j,i} = RO_sol.fy;
            RO_FY{j,i} = RO_sol.fY * system_info.zi;

            TimeRO(j,i) = toc;

            tic
            % Loop for each out-of-sample realization 
            for k = 1:OOS_sim
                system_info.Wreal = WPr(k,:)';
                system_info.DWreal = system_info.Wreal - system_info.mu;

                % Solve real-time optimal power flow
                RT_solution_CVaR = RT_solve_R(system_info,RO_sol.y0,RO_sol.ru,RO_sol.rd);
                RT_Obj_IR(N_vect_idx,j,k,i) = RT_solution_CVaR.Obj_RT;  
                RT_lshed{j,k,i} = RT_solution_CVaR.lshed_RT;
                RT_flow{j,k,i} = RO_sol.fy + RO_sol.fY * system_info.DWreal;
                RT_p{j,k,i} = RO_sol.Y * system_info.DWreal;
                RT_q{j,k,i} = RO_sol.q + RO_sol.qY * system_info.DWreal;
                RT_flag(j,k,i) = RT_solution_CVaR.Flag;

            end

            TimeOOS(j) = toc;

            % Calculation of expected cost

            if RO_sol.Flag == 0
                RO_TC(N_vect_idx, i, j) = system_info.cru'*RO_sol.ru + system_info.crd'*RO_sol.rd + mean(RT_Obj_IR(N_vect_idx,j,:, i));
            else
                RO_TC(N_vect_idx, i, j) = NaN;
            end
        end
    end % end of IR loop
end % end of N loop


Time = toc
