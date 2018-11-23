%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Energy and Reserve Dispatch with Distributionally Robust Joint Chance Constraints
%   Christos Ordoudis, Viet Anh Nguyen, Daniel Kuhn, Pierre Pinson
%
%   This script is part of the data generation process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script generates the dataset of the Wind power production


% Generation of data 
rng(123456,'twister')

% Number of wind farms
wf=[1:6];

% Loading the historical data for wind farms
wff=AV_AEMO2(:,wf);

% Cutting off very extreme values
cut_off_eps = 1e-2;
wff(wff<cut_off_eps) = cut_off_eps;
wff(wff>(1-cut_off_eps)) = 1 - cut_off_eps;

% Logit-normal transformation (Eq. (1) in ref. [31])
yy=log(wff./(1-wff));

% Calculation of mean and variance, note that we increase the mean to have
% higher wind penetration in our test-case
mu = mean(yy)+1.5;
sigma_m=cov(yy);
sigma_m=sigma_m./(std(yy)'*std(yy));

% Inverse of logit-normal transformation (Eq. (2) in ref. [31])
R = chol(sigma_m);
y = repmat(mu,Nscen,1) + randn(Nscen,size(WindDATA,1))*R;
Wind = (1+exp(-y)).^(-1);

% Checking correlation, mean and true mean of data
corrcoef(Wind);
mean(Wind);
true_mean_Wind = (1+exp(-mu)).^(-1);

% Reshaping the data structure
nWind = Wind';
nWind = reshape(nWind,size(WindDATA,1), N_max+OOS_max, IR_max);