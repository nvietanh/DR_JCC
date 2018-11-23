%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Energy and Reserve Dispatch with Distributionally Robust Joint Chance Constraints
%   Christos Ordoudis, Viet Anh Nguyen, Daniel Kuhn, Pierre Pinson
% 
%   This scripts generate the results for the joint chance constraint
%   approach
%   It runs with both the Bonferroni and CVaR approximation and the
%   Optimized CVaR approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all; clc;

% call startup to add the necessary path
startup;

tic

%%

% Data input
RTS_Data2;

% DRO input

% Definition of dual norm used in the Wasserstein metric, value of \epsilon
% and the parameters for the sequential optimization algorithm

DRO_param.dual_norm = 'inf'; % dual norm
DRO_param.eps_joint_cvar = 0.05; % \epsilon
DRO_param.CVaR_max_iter = 40; % MaxIter
DRO_param.tolerance = 1e-1; % \eta
DRO_param.alpha_max = 1000; % \overline{\delta}
DRO_param.alpha_min = 1e-4; % \undeline{\delta}
DRO_param.penalty = 1e6; % BigM

% Number of individual runs (number of coupled datasets in the numerical
% study)

IR_max = 100;
IR_sim = 100;

% Number of out of sample data for each individual run (N') for testing
% dataset
OOS_max = 1000;
OOS_sim = 1000;

% Number of maximum sample size (N)
N_max = 1000;

% Number of sample data in training dataset (N)
N_vect = [25, 50, 100, 200];
% Vector for Bonferroni and CVaR approximation
rho_vectorC = [0 linspace(0.0001, 0.0025, 24)];

% Vector for optimized CVaR approximation
rho_vectorJC = [0 0.0001 linspace(0.001, 0.04, 23)]; 

% Total number of data 
Nscen = IR_max * (N_max + OOS_max);

% Run the script to generate data
Load_All_Data;

% Initializing the matrices to gather final results

Joint_CVaR_Obj_IR = zeros(length(N_vect),IR_sim, OOS_sim, length(rho_vectorC));
CVaR_Obj_IR = zeros(length(N_vect),IR_sim, OOS_sim, length(rho_vectorJC));
ICC_TC = NaN(length(N_vect),IR_sim,length(rho_vectorC));
JCC_TC = NaN(length(N_vect),IR_sim,length(rho_vectorJC));

% Loop for each individual run for 100 coupled datasets
for N_vect_idx = 3:3
    N = N_vect(N_vect_idx);
    for j = 1:IR_sim
        display('out of sample iteration:');
        j

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
        

        % Loop for each value of \rho in P vector
        for i = 1:length(rho_vectorC) 

            % optimize for each value of rho for Bonferroni and CVaR approximation

            DRO_param.rho = rho_vectorC(i);

            DRO_ICC_CVaR = DRO_CVaR_ICC(system_info, DRO_param, jcc);
            ICC_p_DA{j, i} = DRO_ICC_CVaR.p;
            ICC_ru{j, i} = DRO_ICC_CVaR.ru;
            ICC_rd{j, i} = DRO_ICC_CVaR.rd;
            ICC_obj{j, i} = DRO_ICC_CVaR.Obj;
            ICC_flag{j, i} = DRO_ICC_CVaR.Flag;
            CVaR_Y{j,i} = DRO_ICC_CVaR.Y * system_info.xi;
            CVaR_Qy{j,i} = DRO_ICC_CVaR.q;
            CVaR_QY{j,i} = DRO_ICC_CVaR.qY * system_info.xi;
            CVaR_Fy{j,i} = DRO_ICC_CVaR.fy;
            CVaR_FY{j,i} = DRO_ICC_CVaR.fY * system_info.xi;

            % optimize for each value of rho for Optimized CVaR approximation

            DRO_param.rho = rho_vectorJC(i);

            DRO_JCC_CVaR = DRO_JCVaR_All(system_info, DRO_param, jcc);
            JCC_p_DA{j, i} = DRO_JCC_CVaR.y1;
            JCC_ru{j, i} = DRO_JCC_CVaR.ru;
            JCC_rd{j, i} = DRO_JCC_CVaR.rd;
            JCC_obj{j, i} = DRO_JCC_CVaR.Obj;
            JCC_flag{j, i} = DRO_JCC_CVaR.Flag;            
            Joint_CVaR_Y{j,i} = DRO_JCC_CVaR.Y * system_info.xi;
            Joint_CVaR_Qy{j,i} = DRO_JCC_CVaR.q;
            Joint_CVaR_QY{j,i} = DRO_JCC_CVaR.qY * system_info.xi;
            Joint_CVaR_Fy{j,i} = DRO_JCC_CVaR.fy;
            Joint_CVaR_FY{j,i} = DRO_JCC_CVaR.fY * system_info.xi;

            % Loop for each out-of-sample realization 
            for k = 1:OOS_sim
                system_info.Wreal = WPr(k,:)';
                system_info.DWreal = system_info.Wreal - system_info.mu;

                % Solve real-time optimal power flow for the solution of Bonferroni
                % approximation
                RT_solution_CVaR = RT_solve_R(system_info,DRO_ICC_CVaR.p,DRO_ICC_CVaR.ru,DRO_ICC_CVaR.rd);
                CVaR_Obj_IR(N_vect_idx,j,k,i) = RT_solution_CVaR.Obj_RT;  
                CVaR_lshed{j,k,i} = RT_solution_CVaR.lshed_RT;
                CVaR_flow{j,k,i} = DRO_ICC_CVaR.fy + DRO_ICC_CVaR.fY * system_info.DWreal;
                CVaR_p{j,k,i} = DRO_ICC_CVaR.Y * system_info.DWreal;
                CVaR_q{j,k,i} = DRO_ICC_CVaR.q + DRO_ICC_CVaR.qY * system_info.DWreal;
                CVaR_flag(j,k,i) = RT_solution_CVaR.Flag;


                % Solve real-time optimal power flow for the solution of Zymler
                % approximation
                RT_solution_Joint_CVaR = RT_solve_R(system_info,DRO_JCC_CVaR.y1,DRO_JCC_CVaR.ru,DRO_JCC_CVaR.rd);
                Joint_CVaR_Obj_IR(N_vect_idx,j,k,i) = RT_solution_Joint_CVaR.Obj_RT;  
                Joint_CVaR_lshed{j,k,i} = RT_solution_Joint_CVaR.lshed_RT;
                Joint_CVaR_flow{j,k,i} = DRO_JCC_CVaR.fy + DRO_JCC_CVaR.fY * system_info.DWreal;
                Joint_CVaR_p{j,k,i} = DRO_JCC_CVaR.Y * system_info.DWreal;     
                Joint_CVaR_q{j,k,i} = DRO_JCC_CVaR.q + DRO_JCC_CVaR.qY * system_info.DWreal;
                Joint_CVaR_flag(j,k,i) = RT_solution_Joint_CVaR.Flag;


            end

            % Calculation of expected cost

            if DRO_ICC_CVaR.Flag == 0
                ICC_TC(N_vect_idx,j,i) = system_info.cru'*DRO_ICC_CVaR.ru + system_info.crd'*DRO_ICC_CVaR.rd + mean(CVaR_Obj_IR(N_vect_idx,j,:,i));
            else
                ICC_TC(N_vect_idx,j,i) = NaN;
            end

            if DRO_JCC_CVaR.Flag == 0
                JCC_TC(N_vect_idx,j,i) = system_info.cru'*DRO_JCC_CVaR.ru + system_info.crd'*DRO_JCC_CVaR.rd + mean(Joint_CVaR_Obj_IR(N_vect_idx,j,:,i));
            else
                JCC_TC(N_vect_idx,j,i) = NaN;
            end

        end
    end
    save('JCC_resultN100.mat');
end

Time = toc
