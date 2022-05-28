%% NON LINEAR BREAKAGE CODE (2-6-21) (ASHOK DAS) - 2D
% clear all
close all
format longE

example = 3;
%% Inputs
x1_min = 0; x1_max = 1; I1 = 20; % 1st property limits
x2_min = 0; x2_max = 1; I2 = 20; % 1st property limits

T = 0.75;  % [sec] Process time
len_T = 11;
time  = linspace(0,T,len_T); % Time discretization

K_index = 1; % 1-> K=1; 2-> K=x1*x2*y1*y2

%% Initial PSD
N_ini    = zeros(I1*I2,1); % Initialization
N_ini(I1*I2) = 1;

tic
%% Discretizations and other functions
[x1,R1,del_x1] = Grids2(x1_min, x1_max, I1);
[x2,R2,del_x2] = Grids2(x2_min, x2_max, I2); % x-> pivot pts; R-> boundary pts; del_x-> grid length
% [x1,R1,del_x1] = Lin_Grids(x1_min, x1_max, I1);
% [x2,R2,del_x2] = Lin_Grids(x2_min, x2_max, I2); % x-> pivot pts; R-> boundary pts; del_x-> grid length

K = K_Fun(K_index,x1,x2,I1,I2); % K-function (collision freq function)

p1 = p_Fun_mat(x1,R1,I1); % p(i,m) matrix form (1st variable)
p2 = p_Fun_mat(x2,R2,I2); % p(i,m) matrix form (2nd variable)

%% Breakage dist fun related funcitons
B = B_Fun(p1,p2,x1,x2,R1,R2); % matrix version of \int_{R1(i1)}^{p1(i1,m1)} \int_{R2(i2)}^{p2(i2,m2)} b(x1,x2;x1_m1,x2_m2) dx1 dx2; where b=2/x1(m1)*x2(m2)

% breakage distribution function  (FOR CONSERVATIVE APPROACH-1: MID-PT. RULE)
beta = zeros(I1,I2,I1,I2); % Initialization
for i=1:I1
    for j=1:I2
        beta(:,:,i,j) = 2/(x1(i)*x2(j));
    end
end

% Int. of b(x,y;z) = 2/y1*y2   (FOR CONSERVATIVE APPROACH-2: ANALYTIC INT.)
beta_cons = zeros(I1,I2,I1,I2); % Initialization
for i=1:I1
    for j=1:I2
        beta_cons(:,:,i,j) = 2*log(R1(i+1)/R1(i))*log(R2(j+1)/R2(j)) /(del_x1(i)*del_x2(j));
    end
end

%% Functions related to schemes
[w1,w2_b,w2_d] = weights(x1,x2,B); % Weight functions for MC and NPMC

t_sim_0= toc;
%% Solution part
options = odeset('RelTol',1e-6, 'AbsTol',1e-6);
%options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'MaxSTep',1000,'Stat','on');
tic
[T1,N1] = ode45(@discrete_MC, time, N_ini, options, K,B,w1,x1,x2); % Mass conserving technique
t_sim_1 = toc;

tic
[T2,N2] = ode45(@discrete_NPMC, time, N_ini, options, K,B,w2_b,w2_d,x1,x2); % Number conserve + Mass conserving technique
t_sim_2 = toc;

tic
 [T3,N3] = ode45(@discrete_conserve, time, N_ini, options, x1,x2,del_x1,del_x2,K,beta); % Conservative approach-1
t_sim_3 = toc;

tic
[T4,N4] = ode45(@discrete_conserve, time, N_ini, options, x1,x2,del_x1,del_x2,K,beta_cons); % Conservative approach-2
t_sim_4 = toc;


%% Final PSD
% PSD_MC = vec2mat(N1,I1,I2);  PSD_NPMC = vec2mat(N2,I1,I2); PSD_cons = vec2mat(N4,I1,I2);
% [x,y]=meshgrid(x1,x2);
%
% figure
% surf(log(x),log(y),PSD_MC)
% figure
% surf(x,y,PSD_NPMC)
% figure
% surf(x,y,PSD_cons)

%% Figure plot
%% Total particle - M_00
N_tot1 = sum(N1,2); N_tot2 = sum(N2,2); N_tot3 = sum(N3,2); N_tot4 = sum(N4,2);
figure
plot(T1,N_tot1,'bo--','linewidth',1.5,'markersize',11)
hold on
plot(T2,N_tot2,'rs--','linewidth',1.5,'markersize',11)
plot(T3,N_tot3,'m^--','linewidth',1.5,'markersize',11)
% plot(T4,N_tot4,'m^--','linewidth',1.5,'markersize',11)
% Analytical result
plot(time,1./(1-time),'k-','markersize',16,'linewidth',2.5)

xlim([0 T])
legend({'WMC','WMNP','CF','Exact'},'fontsize',18,'Location','northwest')
xlabel('Time (\tau)','fontsize',25);
ylabel('M_{0,0}(\tau)','fontsize',25);
 %savePDF(['Ex_',num2str(example),'_M00'])

%% Total mass - M_1,1
area_mat = x1'*x2; area_mat_vec = mat2vec(area_mat);
M_tot_1 = N1*area_mat_vec; M_tot_2 = N2*area_mat_vec;
M_tot_3 = N3*area_mat_vec; M_tot_4 = N4*area_mat_vec;

figure
plot(T1,M_tot_1,'bo--','linewidth',1.5,'markersize',11)
hold on
plot(T2,M_tot_2,'rs--','linewidth',1.5,'markersize',11)
 plot(T3,M_tot_3,'m^--','linewidth',1.5,'markersize',11)
% plot(T4,M_tot_4,'m^','linewidth',1.5,'markersize',23)

plot(time,sqrt(1-time),'k-','markersize',16,'linewidth',2.5)
% title('Total mass')
xlim([0 T])
legend({'WMC','WMNP','CF','Exact'},'fontsize',18,'Location','best')
%ylim([0.998 1.003])
xlabel('Time (\tau)','fontsize',25);
ylabel('M_{1,1}(\tau)','fontsize',25);
% savePDF(['Ex_',num2str(example),'_M11'])

%% M_{1,0} + M_{0,1}
for i=1:len_T
    N_MC = vec2mat(N1(i,:),I1,I2); N_NPMC = vec2mat(N2(i,:),I1,I2);
    N_cons1 = vec2mat(N3(i,:),I1,I2); N_cons2 = vec2mat(N4(i,:),I1,I2); 
    
    Mix_1(i) = x1*sum(N_MC,2) + sum(N_MC)*x2'; Mix_2(i) = x1*sum(N_NPMC,2) + sum(N_NPMC)*x2';
    Mix_3(i) = x1*sum(N_cons1,2) + sum(N_cons1)*x2'; Mix_4(i) = x1*sum(N_cons2,2) + sum(N_cons2)*x2';
end
figure
plot(T1,Mix_1,'bo','linewidth',1.5,'markersize',18)
hold on
plot(T2,Mix_2,'rs','linewidth',1.5,'markersize',24)
  plot(T3,Mix_3,'m^','linewidth',1.5,'markersize',23)
% plot(T4,Mix_4,'m^--','linewidth',1.5,'markersize',11)

plot(time, 2*ones(1,length(time)),'k-','markersize',16,'linewidth',2.5)
xlim([-0.01 T+0.01])
legend({'WMC','WMNP','CF','Exact'},'fontsize',18,'Location','best')
ylim([1.99 2.01])
xlabel('Time (\tau)','fontsize',25);
ylabel('M_{1,0}(\tau) + M_{0,1}(\tau)','fontsize',25);
% savePDF(['Ex_',num2str(example),'_M10_M01'])

%% Average hypervolume
figure
plot(T1,M_tot_1./N_tot1,'bo--','linewidth',1.5,'markersize',11)
hold on
plot(T2,M_tot_2./N_tot2,'rs--','linewidth',1.5,'markersize',11)
  plot(T3,M_tot_3./N_tot3,'m^--','linewidth',1.5,'markersize',11)
% plot(T4,M_tot_4./N_tot4,'m^--','linewidth',1.5,'markersize',11)
xlim([0 T])
plot(time,(1-time).^(1.5),'k-','markersize',16,'linewidth',2.5)
% title('Total mass')
legend({'WMC','WMNP','CF','Exact'},'fontsize',18,'Location','best')
xlabel('Time (\tau)','fontsize',25);
ylabel('Average hypervolume','fontsize',25);
% savePDF(['Ex_',num2str(example),'_avg_hypervol'])


%% Error table

M00_err_1 = (1./(1-time) - N_tot1')./ (1./(1-time)); M00_err_2 = (1./(1-time) - N_tot2')./ (1./(1-time));
M00_err_3 = (1./(1-time) - N_tot3')./ (1./(1-time)); M00_err_4 = (1./(1-time) - N_tot4')./ (1./(1-time));

M11_err_1 = (sqrt(1-time) - M_tot_1')./ (sqrt(1-time)); M11_err_2 = (sqrt(1-time) - M_tot_2')./ (sqrt(1-time));
M11_err_3 = (sqrt(1-time) - M_tot_3')./ (sqrt(1-time)); M11_err_4 = (sqrt(1-time) - M_tot_4')./ (sqrt(1-time));

Mix_err_1 = (2 - Mix_1)/2; Mix_err_2 = (2 - Mix_2)/2;
Mix_err_3 = (2 - Mix_3)/2; Mix_err_4 = (2 - Mix_4)/2;