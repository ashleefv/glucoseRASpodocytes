%% PARAMETER INITIALIZATION
% set up max and mix matrices

pmin=[1e-2, % VmaxoverKm 
1e-4, % k_cat_Renin
1e-3, % k_feedback
1e-7, % feedback_capacity
1e-5, % k_cons_AngII
1]; % dummy

pmax=[1e-2, % VmaxoverKm 
1e-4, % k_cat_Renin
1e-3, % k_feedback
1e-7, % feedback_capacity
1e-5, % k_cons_AngII
1]; % dummy

% Parameter Labels 
efast_var={'V_{max}/K_M','k_R','k_f','f','k_{AII}','dummy'};%,

% PARAMETER BASELINE VALUES
VmaxoverKm = 1.45e-02;
k_cat_Renin = 6.56e+04;
k_feedback = 6.00e-02;
feedback_capacity = 1.00e+02;
k_cons_AngII = 6.12e-01;
dummy=1;

%% TIME SPAN OF THE SIMULATION
t_end=4000; % length of the simulations
tspan=(0:1:t_end);   % time points where the output is calculated
time_points=[2000 4000]; % time points of interest for the US analysis

% INITIAL CONDITION FOR THE ODE MODEL
T0=1e3;
T1=0;
T2=0;
V=1e-3;

y0=[T0,T1,T2,V];

% Variables Labels
y_var_label={'T','T*','T**','V'};
