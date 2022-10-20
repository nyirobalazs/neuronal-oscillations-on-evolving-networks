RESULTS = load('RESULTS-Networks-Sum.mat');

N = 83;
tspan = 1:30;
left_hemisphere = 1:36;
right_hemisphere = 37:75;
brain_stem = 76;

alpha = 0.95;
sigma_0 = 1;
sigma_f = 1;
sigma_ratio = sigma_0/sigma_f;

betasize = 3;
betas = [0 0.25 1];

PHIS = NaN(betasize, 30);
NET_STATS = NaN(3, betasize, 30);
DEG_T = NaN(N, betasize, 30);

W_OUT = RESULTS.RESULTS{1}(:,(2*N)+1:2*N+N^2);
W_matrix = sparse(reshape(W_OUT(1,:), N, N));
[i,j,v] = find(W_matrix);
EDGE_T = NaN(length(i),betasize,30);
%%
for betalevel = 1:betasize
    
    W_OUT = RESULTS.RESULTS{betalevel}(:,(2*N)+1:2*N+N^2);

    for T = tspan

        W_matrix = sparse(reshape(W_OUT(T,:), N, N)); % turn network into matrix
        G = graph(W_matrix,'upper');

        degrees = sum(W_matrix, 2);
        DEG_T(:,betalevel, T) = degrees;
        [i,j,v] = find(W_matrix);
        EDGE_T(:,betalevel,T) = v;
        
        Q = W_matrix(degrees>0, degrees>0); % connected matrix
        D = sum(Q, 2); % degrees of connected matrix

        DEG = diag(1./D);
        SQRTDEG = diag(1./sqrt(D));
    
        I = diag(ones(1,length(Q)));
        A = alpha * SQRTDEG * Q * SQRTDEG; % normalised adjacency matrix
        L = diag(D) - Q;
        NL = I - A; % normalised laplacian

        % Network statistics
        
%         if (ismember(T,[1,10,20,30]))
%             
%             G_Q = graph(Q,'omitselfloops');
%             LWidths = 0.1*(G_Q.Edges.Weight);%/max(G_Q.Edges.Weight);
%             
%             figure
%             %set(gca,'FontSize',28) % Creates an axes and sets its FontSize to 18
%             p = plot(G_Q,'LineWidth',LWidths);
%             theta = -45; % to rotate 90 counterclockwise
%             R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
%             RX = R*[p.XData;
%                     p.YData];
%             p.XData = RX(1,:);
%             p.YData = RX(2,:);
%             ylim([min(p.YData)*1.1, max(p.YData)*1.1])
%             xlim([min(p.XData)*1.1, max(p.XData)*1.1])
%             highlight(p, left_hemisphere,'NodeColor','r') % subgraph A
%             highlight(p, right_hemisphere,'NodeColor','k') % subgraph B
%             highlight(p, brain_stem,'NodeColor','green') % subgraph B
%             %title(sprintf('$\beta$ = %s, T = %s', num2str(betas(betalevel)), num2str(T)), 'Interpreter','latex')
%             t = title(sprintf('$\\beta = %s, T = %s$', num2str(betas(betalevel)), num2str(T)), 'Interpreter','latex');
%             t.FontSize = 18;
%             %axis equal
% 
%             axis off
%             ax = gca;
%             outerpos = ax.OuterPosition;
%             ti = ax.TightInset; 
%             left = outerpos(1) + ti(1);
%             bottom = outerpos(2) + ti(2);
%             ax_width = outerpos(3) - ti(1) - ti(3);
%             ax_height = outerpos(4) - ti(2) - ti(4);
%             ax.Position = [left bottom ax_width ax_height];
% 
%             saveas(gcf,sprintf('BrainGraph-Beta=%s-T=%s.png', num2str(betas(betalevel)), num2str(T)))
% 
%         end

        k_mean = mean(degrees);
        nlambdas = eig(NL);
        lambdas = eig(L);

        NET_STATS(1, betalevel, T) = k_mean;
        NET_STATS(2, betalevel, T) = nlambdas(2);
        NET_STATS(3, betalevel, T) = lambdas(2);

        % Calculate Phi
        
        ATA = A*A';
        
        if (max(abs(eig(ATA)))>1)
            disp('Numerical error')
            break
        end

        IATA = I - ATA;
        SIATA = inv(IATA);
        
        V_X0XT = (I - A'*(IATA)*A*sigma_ratio);
        DET_V_X0XT = det(V_X0XT);
    
        M0 = diag(SIATA);
        M1 = diag(A);
        M2 = M1.*M1./M0;
        M3 = 1 - sigma_ratio*(M2);
        
        numerator = sum(log(abs(M3)));
        denom = log(DET_V_X0XT);
    
        PHI = (numerator - denom)/2;
    
        PHIS(betalevel, T) = PHI;
    
    end
end
%%
DEG1 = squeeze(DEG_T(:,1,:));

figure
plot(tspan, DEG1)
xlabel('Time, T (years)')
ylabel('Node connectivity ($d_i(T)$)','Interpreter','latex')
t=title('Node connectivity over time, $\beta = 0$', ...
    'Interpreter','latex');
t.FontSize = 18;
meanDEG1 = mean(DEG1);
hold on
plot(tspan, meanDEG1, 'LineWidth',3,'Color','r','LineStyle','--')

saveas(gcf,sprintf('Connectivity-Beta=0.png'))

DEG1norm = DEG1./meanDEG1;

figure
plot(tspan, DEG1norm)
xlabel('Time, T (years)')
ylabel('Relative node connectivity ($d_i(T)/\bar{d}(T)$)','Interpreter','latex')
t=title('Relative node connectivity over time, $\beta = 0$', ...
    'Interpreter','latex');
t.FontSize = 18;
hold on
plot(tspan,ones(size(tspan)), 'LineWidth',3, 'Color','r', 'LineStyle','--')

saveas(gcf,sprintf('NormedConnectivity-Beta=0.png'))
%%
DEG2 = squeeze(DEG_T(:,2,:));

figure
plot(tspan, DEG2)
xlabel('Time, T (years)')
ylabel('Node connectivity ($d_i(T)$)','Interpreter','latex')
t=title('Node connectivity over time, $\beta = 0.25$', ...
    'Interpreter','latex');
t.FontSize = 18;
meanDEG2 = mean(DEG2);
hold on
plot(tspan, meanDEG2, 'LineWidth',3,'Color','r','LineStyle','--')

saveas(gcf,sprintf('Connectivity-Beta=0.25.png'))

DEG2norm = DEG2./meanDEG2;

figure
plot(tspan, DEG2norm)
xlabel('Time, T (years)')
ylabel('Relative node connectivity ($d_i(T)/\bar{d}(T)$)','Interpreter','latex')
t=title('Relative node connectivity over time, $\beta = 0.25$', ...
    'Interpreter','latex');
t.FontSize = 18;
hold on
plot(tspan,ones(size(tspan)), 'LineWidth',3, 'Color','r', 'LineStyle','--')

saveas(gcf,sprintf('NormedConnectivity-Beta=0.25.png'))
%%
DEG3 = squeeze(DEG_T(:,3,:));

figure
plot(tspan, DEG3)
xlabel('Time, T (years)')
ylabel('Node connectivity ($d_i(T)$)','Interpreter','latex')
t=title('Node connectivity over time, $\beta = 1$', ...
    'Interpreter','latex');
t.FontSize = 18;
meanDEG3 = mean(DEG3);
hold on
plot(tspan, meanDEG3, 'LineWidth',3,'Color','r','LineStyle','--')

saveas(gcf,sprintf('Connectivity-Beta=1.png'))

DEG3norm = DEG3./meanDEG3;

figure
plot(tspan, DEG3norm)
xlabel('Time, T (years)')
ylabel('Relative node connectivity ($d_i(T)/\bar{d}(T)$)','Interpreter','latex')
t=title('Relative node connectivity over time, $\beta = 1$', ...
    'Interpreter','latex');
t.FontSize = 18;
hold on
plot(tspan,ones(size(tspan)), 'LineWidth',3, 'Color','r', 'LineStyle','--')

saveas(gcf,sprintf('NormedConnectivity-Beta=1.png'))
%%
% For some reason it freezes here, just run each section manually

EDGE1 = squeeze(EDGE_T(:,1,:));

figure
plot(tspan, EDGE1)
xlabel('Time, T (years)')
ylabel('Edge weight ($w_{ij}(T)$)','Interpreter','latex')
t=title('Edge weight over time, $\beta = 0$', ...
    'Interpreter','latex');
t.FontSize = 18;
meanEDGE1 = mean(EDGE1);
hold on
plot(tspan, meanEDGE1, 'LineWidth',3,'Color','r','LineStyle','--')

saveas(gcf,'Edges-Beta=0.png')

EDGE1norm = EDGE1./meanEDGE1;

figure
plot(tspan, EDGE1norm)
xlabel('Time, T (years)')
ylabel('Relative edge weight ($w_{ij}(T)/\bar{w}(T)$)','Interpreter','latex')
t=title('Relative edge weight over time, $\beta = 0$', ...
    'Interpreter','latex');
t.FontSize = 18;
hold on
plot(tspan,ones(size(tspan)), 'LineWidth', 3, 'Color','r', 'LineStyle','--')

saveas(gcf,'NormedEdges-Beta=0.png')
%%
% For some reason it freezes here, just run each section manually

EDGE2 = squeeze(EDGE_T(:,2,:));

figure
plot(tspan, EDGE2)
xlabel('Time, T (years)')
ylabel('Edge weight ($w_{ij}(T)$)','Interpreter','latex')
t=title('Edge weight over time, $\beta = 0.25$', ...
    'Interpreter','latex');
t.FontSize = 18;
meanEDGE2 = mean(EDGE2);
hold on
plot(tspan, meanEDGE2, 'LineWidth',3,'Color','r','LineStyle','--')

saveas(gcf,'Edges-Beta=0.25.png')

EDGE2norm = EDGE2./meanEDGE2;

figure
plot(tspan, EDGE2norm)
xlabel('Time, T (years)')
ylabel('Relative edge weight ($w_{ij}(T)/\bar{w}(T)$)','Interpreter','latex')
t=title('Relative edge weight over time, $\beta = 0.25$', ...
    'Interpreter','latex');
t.FontSize = 18;
hold on
plot(tspan,ones(size(tspan)), 'LineWidth', 3, 'Color','r', 'LineStyle','--')

saveas(gcf,'NormedEdges-Beta=0.25.png')
%%
EDGE3 = squeeze(EDGE_T(:,3,:));

figure
plot(tspan, EDGE3)
xlabel('Time, T (years)')
ylabel('Edge weight ($w_{ij}(T)$)','Interpreter','latex')
t=title('Edge weight over time, $\beta = 1$', ...
    'Interpreter','latex');
t.FontSize = 18;
meanEDGE3 = mean(EDGE3);
hold on
plot(tspan, meanEDGE3, 'LineWidth',3,'Color','r','LineStyle','--')

saveas(gcf,'Edges-Beta=1.png')

EDGE3norm = EDGE3./meanEDGE3;

figure
plot(tspan, EDGE3norm)
xlabel('Time, T (years)')
ylabel('Relative edge weight ($w_{ij}(T)/\bar{w}(T)$)','Interpreter','latex')
t=title('Relative edge weight over time, $\beta = 1$', ...
    'Interpreter','latex');
t.FontSize = 18;
hold on
plot(tspan,ones(size(tspan)), 'LineWidth', 3, 'Color','r', 'LineStyle','--')

saveas(gcf,'NormedEdges-Beta=1.png')
%%
figure
LAMBDA2 = squeeze(NET_STATS(3,:,:));
plot(tspan, LAMBDA2, 'LineWidth',3)
xlabel('Time, T (years)')
ylabel('\lambda_2(T)')
t=title('Algebraic connectivity over time', ...
    'Interpreter','latex');
t.FontSize = 18;
legend('\beta = 0', '\beta = 0.25', '\beta = 1')

saveas(gcf,'BaseSpectralGap.png')
%%
figure
NLAMBDA2 = squeeze(NET_STATS(2,:,:));
plot(tspan, NLAMBDA2, 'LineWidth',3)
xlabel('Time, T (years)')
ylabel('\lambda_2(T)')
ax = gca;
ax.FontSize = 15;  % Consistent font size for axis title, ticks, etc.
t=title('Algebraic connectivity over time', ...
    'Interpreter','latex');
t.FontSize = 18;
legend('\beta = 0', '\beta = 0.25', '\beta = 1')

saveas(gcf,'SpectralGap.png')
%%
PHIS(betalevel, T) = PHI;

figure
plot(tspan, PHIS, 'LineWidth',3)
xlabel('Time, T (years)')
ylabel('\Phi(T)')
ax = gca;
ax.FontSize = 15;  % Consistent font size for axis title, ticks, etc.
t=title('Phi decays over time as damage accumulates', ...
    'Interpreter','latex');
t.FontSize = 20;
legend('\beta = 0', '\beta = 0.25', '\beta = 1')


saveas(gcf,'Phi.png')
