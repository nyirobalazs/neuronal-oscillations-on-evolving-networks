clear
rng('default')

%% Read csv and build matrices

filename = 'budapest_connectome_3.0_209_0_median_coun_and_length.csv'; % connectivity database
budapest_dataset = readtable(filename);
filename2 = 'budapest_connectome_3.0_209_0_median.csv'; % distance database
budapest_dist_dataset = readtable(filename2);

N = 83;

%%% Assemble the weight matrix by taking the mean weight between two parent
%%% nodes

edgeWeight = table2array(grpstats(budapest_dataset, ["parentIdNode1", "parentIdNode2"], "mean", "DataVars", "edgeWeight_fibreCount_lengthMedian_")); % group by from, to and take mean
WEIGHT_MATRIX = sparse(edgeWeight(:,1), edgeWeight(:,2), edgeWeight(:,4)); % only need the from, to and weight
WEIGHT_MATRIX = WEIGHT_MATRIX + WEIGHT_MATRIX' - diag(diag(WEIGHT_MATRIX)); % ensure symmetry
%WEIGHT_MATRIX = WEIGHT_MATRIX - diag(diag(WEIGHT_MATRIX)); % remove self-weights

%%% Now get the distance matrix between nodes by taking the mean fibre
%%% length between them

edgeDelays = grpstats(budapest_dist_dataset, ["parentIdNode1", "parentIdNode2"], "mean", "DataVars", "edgeWeight_medFlm_"); % group by from, to and take mean
[tau_levels, tau_bins] = discretize(edgeDelays.mean_edgeWeight_medFlm_/150, 40); % convert distance in cm to time delay in s (using transmission speed of 150 cm/s), then discretize accoridng to instructions in Appendix
edgeDelays.tau_levels = tau_levels; % use the discretised buckets for the delays (tau)
edgeDelaysMat = table2array(edgeDelays(:, [1,2,5])); % only need the from, to, and bucket number
edgeDelaysRevMat = table2array(edgeDelays(:, [2,1,5])); % get the reverse edges (symmetric)
edgeDelaysRevMat = edgeDelaysRevMat(edgeDelaysRevMat(:,1) ~= edgeDelaysRevMat(:,2),:); % remove the self-edges frpom the reversa to avoid double-counting
edgeDelays = [edgeDelaysMat; edgeDelaysRevMat]; % all the edges between pairs

%% Calculate reatve change [Original code to compute network decay over time, with only beta = 0.25]

%define variables
c_array         = zeros([1,N]);
c_array(26)     = 0.025;
c_array(68)     = 0.025;
q_array         = zeros([1,N]);
tspan           = 1:30;
beta            = [0.25]; % only kept the value of beta used for the other images

% ODE function
input           = [c_array; q_array;WEIGHT_MATRIX];
RESULTS         = {};

for BETA = beta
    [t, out]        = ode45(@(t,C0) yearly_chang_ode(t,C0,BETA),tspan,input');
    RESULTS{end+1}  = out;
end

save(sprintf('RESULTS-Beta=%s-FixedWithSelfWeights.mat', num2str(0.25)), 'RESULTS', '-v7.3')

%% Create plots

% FIGURE 2
ind = 1;
linestyles = {'-'; ':'; '--'};

for OUT = RESULTS
    %Separate the outputs
    C_OUT                = OUT{1,1}(:,1:N);
    Q_OUT                = OUT{1,1}(:,N+1:2*N);
    W_OUT                = OUT{1,1}(:,(2*N)+1:2*N+N^2);
    
    %Calculate and plot average concentration
    C_T                  = (1/N)*sum(normalize(C_OUT,'range'),2);
    fig                  = plot(t,C_T);
    fig.LineStyle        = linestyles{ind};
    fig.Color            = 'b';
    fig.LineWidth        = 1;
    hold on
    
    %Calculate and plot average damage
    Q_T                  = (1/N)*sum(normalize(Q_OUT,'range'),2);
    fig                  = plot(t,Q_T);
    fig.LineStyle        = linestyles{ind};
    fig.Color            = 'r';
    fig.LineWidth        = 1;
    hold on
    
    %Calculate and plot scaled average connection weight
    w_t                  = (1/(N^2)).*sum(normalize(W_OUT,'range'),2);
    w_0                  = (1/(N^2)).*sum(sum(normalize(WEIGHT_MATRIX,'range'),2),1);
    W_T                  = normalize(w_t/w_0,'range');
    fig                  = plot(t,W_T);
    fig.LineStyle        = linestyles{ind};
    fig.Color            = 'g';
    fig.LineWidth        = 1;
    hold on

    ind = ind + 1;
end

saveas(gcf,'Fig1-FixedWithSelfWeights.png')

% few extra setting for fig. 2
title('Evolution of averaged toxic concentration and damage (\\beta = 0.25)')
xlabel('Time T (yr)')
xlim([1 30])
ylim([0.005 1.05])


% FIGURE 3

%Defining brain areas
frontal         = [1, 4:11, 42, 45:54];
parietal        = [16:20, 22, 55:61];
occipital       = [23, 63:65];
temporal        = [25:27, 29:32, 34:35,40, 66:76, 81];
limbic          = [12:15, 82];
basa_ganglia    = [36:39, 77:80];
brain_stem      = 83;

% Calculate Q in every brain regions
Q_FRONT         = normalize(sum(RESULTS{1,1}(:,N + frontal),2),'range');
Q_PARI          = normalize(sum(RESULTS{1,1}(:,N + parietal),2),'range');
Q_OCC           = normalize(sum(RESULTS{1,1}(:,N + occipital),2),'range');
Q_TEMP          = normalize(sum(RESULTS{1,1}(:,N + temporal),2),'range');
Q_LIMB          = normalize(sum(RESULTS{1,1}(:,N + limbic),2),'range');
Q_BASG          = normalize(sum(RESULTS{1,1}(:,N + basa_ganglia),2),'range');
Q_BS            = normalize(sum(RESULTS{1,1}(:,N + brain_stem),2),'range');

% Make the plot
figure(2)
plot(tspan,Q_FRONT,tspan,Q_PARI,tspan,Q_OCC,tspan,Q_TEMP,tspan,Q_LIMB,tspan,Q_BASG,tspan,Q_BS)
title('Evolution of damage in different brain regions')
xlabel('Time T (yr)')
ylabel('Q_{s}')
xlim([1 30])
ylim([0.005 1.05])
legend('Frontal','Parietal','Occipital','Temporal','Limbic','Basal ganglia','Brain stem')

saveas(gcf,'Fig2-FixedWithSelfWeights.png')

%% Calculate biomarkers [NEW CODE STARTS HERE]

%% Parameters for Independent Dynamics

params.nmodules = N;
params.lambda = -0.01;

params.omega = 40;

%% Parameters for Connection Dynamics

params.kappa = 10;
params.tau = mean(tau_bins(edgeDelays(:,3))); % get the average of all delays (not used)
params.edgeDelays = edgeDelays; % store the delays

t_sim = 10; % time to run simulation for
t_warmup = 3; % time to warm up the simulation for
niter = 12; % number of iterations

betalevel = beta; % dummy

W_OUT = RESULTS{1,1}(:,(2*N)+1:2*N+N^2); % get all the network data across all 30 years (1 row for each year)

A_VECTORS = zeros(niter, tspan(end), N); % initialise matrix to store all A values (1 column for each year and 1 row for each iteration, third dimension is for each node)
sols = cell(niter, tspan(end)); % store the ODE solutions here (1 column for each year and 1 row for each iteration)

A_ALLs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)
A_FRONTSs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)
A_PARIs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)
A_OCCs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)
A_TEMPs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)
A_LIMBs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)
A_BASGs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)
A_BSs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)

Bs = zeros(niter, tspan(end)); % initialise matrix to store B values (1 column for each year and 1 row for each iteration)

B_ALLs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)
B_FRONTSs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)
B_PARIs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)
B_OCCs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)
B_TEMPs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)
B_LIMBs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)
B_BASGs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)
B_BSs = zeros(niter, tspan(end)); % initialise matrix to store A values (1 column for each year and 1 row for each iteration)

tic % time the code
for iter = 1:niter
    parfor T = tspan
        fprintf('Iter = %s, Year = %s\n', num2str(iter), num2str(T));

        W_matrix = sparse(reshape(W_OUT(T,:), N, N)); % turn network into matrix
        delta_k = normrnd(0, sqrt(0.1), 1, N); % get random values for "intrinsic frequencies"
        omega_k = params.omega + delta_k; % add to the baseline frequency

        [A, B, sol] = getBiomarkers(W_matrix, params, omega_k, tau_bins, t_warmup, t_sim); % get biomarkers by running dde23
        
        % calculate means and save all the values

        sols{iter, T} = sol;

        A_ALL = mean(A);
        A_FRONT = mean(A(frontal));
        A_PARI = mean(A(parietal));
        A_OCC = mean(A(occipital));
        A_TEMP = mean(A(temporal));
        A_LIMB = mean(A(limbic));
        A_BASG = mean(A(basa_ganglia));
        A_BS = mean(A(brain_stem));

        A_VECTORS(iter, T, :) = A;
        A_ALLs(iter, T) = A_ALL;
        A_FRONTs(iter, T) = A_FRONT;
        A_PARIs(iter, T) = A_PARI;
        A_OCCs(iter, T) = A_OCC;
        A_TEMPs(iter, T) = A_TEMP;
        A_LIMBs(iter, T) = A_LIMB;
        A_BASGs(iter, T) = A_BASG;
        A_BSs(iter, T) = A_BS;

        B_ALL = mean(B);
        B_FRONT = mean(B(frontal));
        B_PARI = mean(B(parietal));
        B_OCC = mean(B(occipital));
        B_TEMP = mean(B(temporal));
        B_LIMB = mean(B(limbic));
        B_BASG = mean(B(basa_ganglia));
        B_BS = mean(B(brain_stem));

        B_VECTORS(iter, T, :) = B;
        B_ALLs(iter, T) = B_ALL;
        B_FRONTs(iter, T) = B_FRONT;
        B_PARIs(iter, T) = B_PARI;
        B_OCCs(iter, T) = B_OCC;
        B_TEMPs(iter, T) = B_TEMP;
        B_LIMBs(iter, T) = B_LIMB;
        B_BASGs(iter, T) = B_BASG;
        B_BSs(iter, T) = B_BS;

    end

    % Save available results after each iterations in case it crashes

    save(sprintf('A-VECTOR-Beta=%s-FixedWithSelfWeights.mat', num2str(betalevel)), 'A_VECTORS')
    save(sprintf('B-VECTOR-Beta=%s-FixedWithSelfWeights.mat', num2str(betalevel)), 'B_VECTORS')

    writematrix(A_ALLs,sprintf('As-Beta=%s-FixedWithSelfWeights.csv', num2str(betalevel)))
    writematrix(B_ALLs,sprintf('Bs-Beta=%s-FixedWithSelfWeights.csv', num2str(betalevel)))
    
    if (iter > 1)

        % Plot and save available results after each iterations in case it crashes
        
        % A(T), unnormalized

        figure
        hold on
        
        thick = true;
        
        for A_MATRIX = {A_ALLs, A_FRONTs, A_PARIs, A_OCCs, A_TEMPs, A_LIMBs, A_BASGs}
            
            A_MATRIX_NONZERO = A_MATRIX{1};
            A_MATRIX_NONZERO = A_MATRIX_NONZERO(1:iter,:);
            A_mean = mean(A_MATRIX_NONZERO, 1);
            A_sd = std(A_MATRIX_NONZERO, 1);
            if thick
                errorbar(tspan, A_mean, A_sd, 'LineWidth', 3, 'Color', 'k')
                thick = false;
            else
                errorbar(tspan, A_mean, A_sd)
            end
            hold on
        end
        
        legend('Overall','Frontal','Parietal','Occipital','Temporal','Limbic','Basal ganglia')
        title(sprintf('Evolution of average oscillatory activity over time, \\beta=%s',num2str(betalevel)))
        xlabel('Time T (yr)')
        ylabel('A(T)')
        
        saveas(gcf,sprintf('AT-Beta=%s-FixedWithSelfWeights.png',num2str(betalevel)))

        % A(T), normalized

        figure
        hold on
        
        thick = true;

        for A_MATRIX = {A_ALLs, A_FRONTs, A_PARIs, A_OCCs, A_TEMPs, A_LIMBs, A_BASGs}
            
            A_MATRIX_NONZERO = A_MATRIX{1};
            A_MATRIX_NONZERO = A_MATRIX_NONZERO(1:iter,:);
            initMean = mean(A_MATRIX_NONZERO(:,1));
            A_MATRIX_NONZERO = A_MATRIX_NONZERO/initMean;

            A_mean = mean(A_MATRIX_NONZERO, 1);
            A_sd = std(A_MATRIX_NONZERO, 1);
            if thick
                errorbar(tspan, A_mean, A_sd, 'LineWidth', 3, 'Color', 'k')
                thick = false;
            else
                errorbar(tspan, A_mean, A_sd)
            end
            hold on
        end
        
        ylim([0 1.5])
        legend('Overall','Frontal','Parietal','Occipital','Temporal','Limbic','Basal ganglia')
        title(sprintf('Evolution of average oscillatory activity over time, \\beta=%s',num2str(betalevel)))
        xlabel('Time T (yr)')
        ylabel('A(T)')
        
        saveas(gcf,sprintf('AT-Normed-Beta=%s-FixedWithSelfWeights.png',num2str(betalevel)))

        % B(T), unnormalized
        
        figure
        hold on
        
        thick = true;
        
        for B_MATRIX = {B_ALLs, B_FRONTs, B_PARIs, B_OCCs, B_TEMPs, B_LIMBs, B_BASGs}
            
            B_MATRIX_NONZERO = B_MATRIX{1};
            B_mean = mean(B_MATRIX_NONZERO(1:iter,:), 1);
            B_sd = std(B_MATRIX_NONZERO(1:iter,:), 1);
            if thick
                errorbar(tspan, B_mean, B_sd, 'LineWidth', 3, 'Color', 'k')
                thick = false;
            else
                errorbar(tspan, B_mean, B_sd)
            end
            hold on
        end
        
        legend('Overall','Frontal','Parietal','Occipital','Temporal','Limbic','Basal ganglia')
        title(sprintf('Evolution of metastability index over time, \\beta=%s',num2str(betalevel)))
        xlabel('Time T (yr)')
        ylabel('B(T)')
        
        saveas(gcf,sprintf('BT-Beta=%s-FixedWithSelfWeights.png',num2str(betalevel)))

        % B(T), normalized

        figure
        hold on
        
        thick = true;

        for B_MATRIX = {B_ALLs, B_FRONTs, B_PARIs, B_OCCs, B_TEMPs, B_LIMBs, B_BASGs}
            
            B_MATRIX_NONZERO = B_MATRIX{1};
            B_MATRIX_NONZERO = B_MATRIX_NONZERO(1:iter,:);
            initMean = mean(B_MATRIX_NONZERO(:,1));
            B_MATRIX_NONZERO = B_MATRIX_NONZERO/initMean;
            
            B_mean = mean(B_MATRIX_NONZERO, 1);
            B_sd = std(B_MATRIX_NONZERO, 1);
            if thick
                errorbar(tspan, B_mean, B_sd, 'LineWidth', 3, 'Color', 'k')
                thick = false;
            else
                errorbar(tspan, B_mean, B_sd)
            end
            hold on
        end
        
        ylim([0 1.5])
        legend('Overall','Frontal','Parietal','Occipital','Temporal','Limbic','Basal ganglia')
        title(sprintf('Evolution of metastability index over time, \\beta=%s',num2str(betalevel)))
        xlabel('Time T (yr)')
        ylabel('B(T)')
        
        saveas(gcf,sprintf('BT-Normed-Beta=%s-FixedWithSelfWeights.png',num2str(betalevel)))

    end
end
toc

save(sprintf('sols-Beta=%s-FixedWithSelfWeights.mat', num2str(betalevel)), 'sols', '-v7.3')

%% Helper functions

function output = yearly_chang_ode(t,input,beta)

    % Constants
    alpha = 0.75;
    gamma = beta/2;
    p = 0.01;
    
    N = 83;
    C = input(1:N);
    Q = input(N+1:2*N);
    W = input(2*N+1:2*N+N^2);
    Q_PLUS = Q + Q';
    W=reshape(W,N,N);
    
    %Calculate Laplacian matrix
    diag_matrix = diag(sum(W,2));
    L = p*(diag_matrix-W);
    
    %Calculations
    
    CDOT     = -L*C+alpha*C.*(1-C);
    QDOT     = beta*C.*(1-Q);
    WDOT     = -gamma*W.*Q_PLUS;
    
    output     = [CDOT;QDOT;reshape(WDOT,N^2,1)];

end

function [A,B,sol] = getBiomarkers(W_matrix, params, omega_k, tau_bins, t_warmup, t_sim)

    params.W_matrix = W_matrix;
    params.omega_k = omega_k;
    
    X0 = 2*rand(1, params.nmodules) - 1; % random initial values of excitatory population (array of [nmodules])
    Y0 = 2*rand(1, params.nmodules) - 1; % random intiial values of inhibitory population (array of [nmodules])
    Z0 = X0' + Y0'*1i; % put them together into a complex number (array of [nmodules] complex)

    sol = dde23(@(t, Z, Xhist) dXdt(t, Z, Xhist, params), tau_bins, @(t) history(t, Z0), [0,t_sim+t_warmup]); % run DDE (required delays comes from the tau_bins vector)

    t = sol.x;
    X = sol.y;
    Z = X(1:params.nmodules,:); % get the complex values
    G = X(params.nmodules+1:end,:); % get the integrals of the absolute values of the complex values

    % Calculate A(T)

    G_final_value = G(:,end); % get the final integral value
    G_warmup = G(:,t < t_warmup); 
    G_warmup_final_value = G_warmup(:,end);
    G_steadystate = G_final_value - G_warmup_final_value; % the integral over the simulation should be between the end and from after the warmup
    A = G_steadystate/t_sim; % get average across time

    % Calculate B(T)

    aZ = abs(Z(:,t > t_warmup)); % discard the warmup activity
    B = var(aZ, 0, 2); % get variation across time
    
end


function dXdt = dXdt(t, X, Xhist, params)
   
    Z = X(1:params.nmodules); % complex numbers
    G = real(X(params.nmodules+1:end)); % absolute values of complex numbers

    Zhist = Xhist(1:params.nmodules,:); % get the history of the complex numbers
    
    % independent dynamics
    
    lambda = params.lambda;
    omega_k = params.omega_k;
    
    F = Z.*(lambda + omega_k'*1i - abs(Z).^2);

    % network dynamics

    kappa = params.kappa;
    W_matrix = params.W_matrix;
    edgeDelays = params.edgeDelays;

    z_kt = Zhist(sub2ind(size(Zhist), edgeDelays(:,2), edgeDelays(:,3))); % this is a shortcut to calculate the complex value of a neighbour k at some specified time t reflecting a delay
    S = sparse(edgeDelays(:,1), edgeDelays(:,2), z_kt); % turns the above into an N x N matrix (where i, j is the value of Z for neighbour j when it arrives to i after the delay)
    
    D1 = sum(W_matrix.*S, 2); % row sums of the element-wise product
    D2 = real(D1); % get reals
    D3 = sigmoid(D2); % pass through sigmoid

    dZdt = F + kappa.*D3; % combine the dynamics
    dGdt = abs(Z); % record the absolute value for the integral trick

    dXdt = [dZdt; dGdt];

end

function Xhist = history(t, Z0)
    Xhist = [Z0; zeros(size(Z0))]; % for X < 0, use the initial values for the complex variables and zeros for the absolute values
end

function S = sigmoid(x)
    
    S = 1./(1+exp(-x));

end