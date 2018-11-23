function [ jcc ] = CC_matrices(si, DRO_param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Energy and Reserve Dispatch with Distributionally Robust Joint Chance Constraints
%   Christos Ordoudis, Viet Anh Nguyen, Daniel Kuhn, Pierre Pinson
%
%   This function calculate the matrices of the joint chance constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Matrices for generation 

A_g = [zeros(length(si.Pmax)),-eye(length(si.Pmax)),zeros(length(si.Pmax));
       zeros(length(si.Pmax)),zeros(length(si.Pmax)),-eye(length(si.Pmax))];
   
B_g = [eye(length(si.Pmax));
       -eye(length(si.Pmax))];  
   
C_g = [zeros(length(si.Pmax),length(si.Wmax));
       zeros(length(si.Pmax),length(si.Wmax))]; 
   
d_g = [zeros(length(si.Pmax),1);
       zeros(length(si.Pmax),1)];   

% Matrices for transmission lines 

A_l = [si.Qg,zeros(length(si.F),length(si.Pmax)),zeros(length(si.F),length(si.Pmax));
       -si.Qg,zeros(length(si.F),length(si.Pmax)),zeros(length(si.F),length(si.Pmax))];
   
B_l = [si.Qg;
       -si.Qg];  
   
C_l = [si.Qw*si.DiagWmax;
       -si.Qw*si.DiagWmax]; 
   
d_l = [si.F - si.Qw*si.DiagWmax*si.mu + si.Qd*si.D;
       si.F + si.Qw*si.DiagWmax*si.mu - si.Qd*si.D];  
   
% Matrices for pipelines 

A_m = [si.PG,zeros(length(si.FP),length(si.Pmax)),zeros(length(si.FP),length(si.Pmax));
       -si.PG,zeros(length(si.FP),length(si.Pmax)),zeros(length(si.FP),length(si.Pmax))];
   
B_m = [si.PG;
       -si.PG];  
   
C_m = [zeros(length(si.FP),length(si.Wmax));
       zeros(length(si.FP),length(si.Wmax))]; 
   
d_m = [si.FP;
       zeros(length(si.FP),1)];  
   
   
% In the hybrid formulation found in the e-companion, we need more matrices
% for reserve increase, windspill and loadshedding

% for the generator constraints
B_g2 = zeros(2*length(si.Pmax), length(si.Pmax));
B_g3 = zeros(2*length(si.Pmax), length(si.Wmax));
B_g4 = zeros(2*length(si.Pmax), length(si.D));

% for the grid constraints
B_l2 = zeros(2*length(si.F), length(si.Pmax));
B_l3 = [-si.Qw; si.Qw]; 
B_l4 = [si.Qd; -si.Qd];

% for the gas constraints
B_m2 = zeros(2*length(si.PG), length(si.Pmax));
B_m3 = zeros(2*length(si.PG), length(si.Wmax));
B_m4 = zeros(2*length(si.PG), length(si.D));
   
jcc{1,1} = A_g;
jcc{2,1} = A_l;
jcc{3,1} = A_m;

jcc{1,2} = B_g;
jcc{2,2} = B_l;
jcc{3,2} = B_m;

jcc{1,3} = C_g;
jcc{2,3} = C_l;
jcc{3,3} = C_m;

jcc{1,4} = d_g;
jcc{2,4} = d_l;
jcc{3,4} = d_m;

jcc{1,5} = DRO_param.eps_joint_cvar;
jcc{2,5} = DRO_param.eps_joint_cvar;
jcc{3,5} = DRO_param.eps_joint_cvar;

jcc{1,6} = B_g2;
jcc{2,6} = B_l2;
jcc{3,6} = B_m2;

jcc{1,7} = B_g3;
jcc{2,7} = B_l3;
jcc{3,7} = B_m3;

jcc{1,8} = B_g4;
jcc{2,8} = B_l4;
jcc{3,8} = B_m4;

end

