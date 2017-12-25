% First order and total effect indices for a given
% model computed with Extended Fourier Amplitude
% Sensitivity Test (EFAST).
% Andrea Saltelli, Stefano Tarantola and Karen Chan.
% 1999. % "A quantitative model-independent method for global
% sensitivity analysis of model output". % Technometrics 41:39-56
close all;
%% INPUT
NR = 5; %: no. of search curves - RESAMPLING
%k = 5 + 1; % # of input factors (parameters varied) + dummy parameter
k = 8 + 1;
NS = 65; % # of samples per search curve
wantedN=NS*k*NR; % wanted no. of sample points

% OUTPUT
% SI[] : first order sensitivity indices
% STI[] : total effect sensitivity indices
% Other used variables/constants:
% OM[] : vector of k frequencies
% OMi : frequency for the group of interest
% OMCI[] : set of freq. used for the compl. group
% X[] : parameter combination rank matrix
% AC[],BC[]: fourier coefficients
% FI[] : random phase shift
% V : total output variance (for each curve)
% VI : partial var. of par. i (for each curve)
% VCI : part. var. of the compl. set of par...
% AV : total variance in the time domain
% AVI : partial variance of par. i
% AVCI : part. var. of the compl. set of par.
% Y[] : model output

MI = 4; %: maximum number of fourier coefficients
% that may be retained in calculating the partial
% variances without interferences between the
% assigned frequencies

%% PARAMETERS AND ODE SETTINGS (they are included in the following file)
% Parameter_settings_EFAST;

time_points=1; %save only 24 hours
% Computation of the frequency for the group
% of interest OMi and the # of sample points NS (here N=NS)
OMi = floor(((wantedN/NR)-1)/(2*MI)/k);
NS = 2*MI*OMi+1
if(NS*NR < 65)
    fprintf(['Error: sample size must be >= ' ...
    '65 per factor.\n']);
    return;
end


%% Pre-allocation of the output matrix Y
%% Y will save only the points of interest specified in
%% the vector time_points
number_output_vars = 1;
Y(NS,length(time_points),number_output_vars,length(pmin),NR)=0;
% Y(NS,length(time_points),length(y0),length(pmin),NR)=0;  % pre-allocation
% Loop over k parameters (input factors)
for i=1:k % i=# of replications (or blocks)
    % Algorithm for selecting the set of frequencies.
    % OMci(i), i=1:k-1, contains the set of frequencies
    % to be used by the complementary group.
    OMci = SETFREQ(k,OMi/2/MI,i);   
    % Loop over the NR search curves.
    for L=1:NR
        % Setting the vector of frequencies OM
        % for the k parameters
        cj = 1;
        for j=1:k
            if(j==i)
                % For the parameter (factor) of interest
                OM(i) = OMi;
            else
                % For the complementary group.
                OM(j) = OMci(cj);
                cj = cj+1;
            end
        end
        % Setting the relation between the scalar
        % variable S and the coordinates
        % {X(1),X(2),...X(k)} of each sample point.
        FI = rand(1,k)*2*pi; % random phase shift
        S_VEC = pi*(2*(1:NS)-NS-1)/NS;
        OM_VEC = OM(1:k);
        FI_MAT = FI(ones(NS,1),1:k)';
        ANGLE = OM_VEC'*S_VEC+FI_MAT;
        
        X(:,:,i,L) = 0.5+asin(sin(ANGLE'))/pi; % between 0 and 1
        
        % Transform distributions from standard
        % uniform to general.
        X(:,:,i,L) = parameterdist(X(:,:,i,L),pmax,pmin,0,1,NS,'unif'); %%this is what assigns 'our' values rather than 0:1 dist
        % Do the NS model evaluations.
        for run_num=1:NS
            [i run_num L] % keeps track of [parameter run NR]
%             % ODE system file
%             f=@ODE_efast;
%             % ODE solver call
%             [t,y]=ode15s(@(t,y)f(t,y,X(:,:,i,L),run_num),tspan,y0,[]); 
        %% transform X into coefficients
        
            coefficients(1) = X(run_num,1,i,L);
            coefficients(2) = X(run_num,2,i,L);
            coefficients(3) = X(run_num,3,i,L);
            coefficients(4) = X(run_num,4,i,L);
            coefficients(5) = X(run_num,5,i,L);
            coefficients(6) = X(run_num,6,i,L);
            coefficients(7) = X(run_num,7,i,L);
            coefficients(8) = X(run_num,8,i,L);
           % coefficients(5) = X(run_num,5,i,L);
           GLU = 5;
           scenario = 0;
           printoutput = 0;
            Ycalc = glucoseRASss(coefficients,GLU,baseline,scenario,printoutput);%run for 2 hours
            %Ycalc = Ycalc_weighted.*sigma;
            % It saves only the output at the time points of interest
            Y(run_num,:,:,i,L)=Ycalc;%(:,time_points+1)';

        end %run_num=1:NS
    end % L=1:NR
end % i=1:k
save Model_efast.mat;
% CALCULATE Si AND STi for each resample (1,2,...,NR) [ranges]
[Si,Sti,rangeSi,rangeSti] = efast_sd(Y,OMi,MI,time_points,1:number_output_vars)
% Si and Sti are 3D matrices: (input variables, time points, output
% variables).
% Calculate Coeff. of Var. for Si and STi for output variables. See
% online Supplement A.5 for details.
[CVsiAngII CVstiAngII]=CVmethod(Si, rangeSi,Sti,rangeSti,1)
% [CVsiAngI CVstiAngI]=CVmethod(Si, rangeSi,Sti,rangeSti,2)
% [CVsiPRA CVstiPRA]=CVmethod(Si, rangeSi,Sti,rangeSti,3)
% T-test on Si and STi for output variables
y_var_label={'Ang II'}%,'Ang I','PRA'};
%p<0.01
s_AngII = efast_ttest(Si,rangeSi,Sti,rangeSti,1:length(time_points),efast_var,1,y_var_label,0.01)
% s_AngI = efast_ttest(Si,rangeSi,Sti,rangeSti,1:length(time_points),efast_var,2,y_var_label,0.01)
% s_PRA = efast_ttest(Si,rangeSi,Sti,rangeSti,1:length(time_points),efast_var,3,y_var_label,0.01)

save Model_efast.mat;

