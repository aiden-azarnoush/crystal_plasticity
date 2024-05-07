%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CP_Project.m
% Authors:
%         Aiden Azarnoush 
%         Michael  Garcia
%         Aisha    Gautam
% E-mail: 
%         mazarnou@purdue.edu
%         michael.e.garcia@asu.edu
%         aisha.s.gautam@gmail.edu
% Description:
% This project involves the development of a point integrator capable of 
% evaluating evolution of flow stress in copper single crystals under 
% uniaxial strain conditions for different crystal orientations: <100>,
% <110>, <310>, and <580> and large strain kinematics. This needs to be
% done for different hardening formulations, including the ones used by 
% Asaro, Luscher et al, and Bassani and Wu.
%
% -------------------------    Functions    -------------------------------
%
% FlowRule    : Compute gamma dot based on Luscher et al. definition
% Hab_asaro   : Compute Asaro and Needleman hardening moduli
% Hab_bassani : Compute Bassani and wu hardening moduli
% Hab_luscher : Compute Luscher et al.  hardening moduli
% macaulay    : Calculate Maculauy bracket i.e. <x> = 1/2(x+|x|)
% rotate_matrix   : Create a orthonrmal basis along specified direction
% rotate_c        : Rotate Tangent Moduli C
% to_voigt_stress : Change input stress tensor to voigt vector
% to_voigt_strain : Change input strain tensor to voigt vector
% to_tensor       : Change input stress vector in voigt to tensor
% plot_tau_gamma     : plot tau versus gamma for 12 each slip system
% plot_stress_strain : plot normal stresses vs logarithmic strain
% plot_Sp            : plot Sp vs time for 12 each slip system
% plot_gamma_dot     : plot gamma_dot vs time for each 12 slip system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CP_Project

close all;
clc;
fprintf('Crystal Orientation: \n ')
fprintf('(1) <100> \n (2) <110>\n (3) <310>\n (4) <580>\n (5) other \n')
x = input('> ');
if x == 1
    direction = [1,0,0];
elseif x == 2
    direction = [1,1,0];
elseif x == 3
    direction = [3,1,0];
elseif x == 4
    direction = [5,8,0];
elseif x == 5
    fprintf('Enter Crystal Orientation: (e.g. [1,1,1]) \n ')
    direction = input('>');
end
R = rotate_matrix(direction);
% Slip Planes
n0 = [  1   1  -1;
        1   1  -1;
        1   1  -1;
        1  -1   1;
        1  -1   1; 
        1  -1   1; 
        1  -1  -1;
        1  -1  -1;
        1  -1  -1;
        1   1   1;
        1   1   1;
        1   1   1];
n0 = n0 ./ norm(n0(1,:));

% Slip direction
m0 = [ 1   0   1;
       0   1   1;
       1  -1   0;
       1   1   0;
       0   1   1;
       1   0  -1;
       1   0   1;
       1   1   0;
       0   1  -1;
       1   0  -1;
       0  -1   1;
       1  -1   0];
 m0 = m0 ./norm(m0(1,:));

 n0 = n0 * R';
 m0 = m0 * R';

%------------------------ Materials parameters----------------------------%
const.gamma_dot0 = 1d7;    % 1/sec
const.E0 = 1d-18;          % J 
const.Sl = 20;             % MPa
const.m1 = 0.33;           % Dimensionless 
const.m2 = 1.66;           % Dimensionless   
%------------------------ Stifness parameters-----------------------------%
const.C11 = 169d3;         % MPa
const.C12 = 122d3;         % MPa 
const.C44 = 75.3d3;        % MPa

 % Rotate Tangent Stiffness moduli
C = rotate_c(R,const);

fprintf('Hardening model: \n ')
fprintf(' 1-Asaro and Needleman Hardening \n')
fprintf('  2-Bassani and Wu Hardening \n')
fprintf('  3-Luscher et al. Hardening \n')
load_case = input('> ');

% Discertize time domain
dt = 1d-12;
t  = 10e-9;
step = round(t/dt);

% Define total deformation gradient
max = 1.1;
lambda = linspace(1, max, step+1);

% Intial conditions
tau    = zeros(12,1);
Sp     = zeros(12,1);
Fp     = eye(3,3);
gamma  = zeros(12,1);
Sigma  = zeros(6,step); 
E      = zeros(6,step);
store.Sp        = zeros(12,step);
store.tau       = zeros(12,step);
store.gamma_dot = zeros(12,step);
store.gamma     = zeros(12,step);

for ts = 1 : step
    
    gamma_dot = FlowRule(tau,Sp,const);
    
    if load_case ==1
      Hab = Hab_asaro(gamma);
    elseif load_case ==2
      Hab = Hab_bassani(gamma);
    else
      Hab = Hab_luscher(gamma_dot,Sp);
    end
    Sp_dot = Hab * abs(gamma_dot);
    
    % Update gamma_alpha and Sp_alpha
    Sp    = Sp + Sp_dot * dt;
    gamma = gamma + gamma_dot * dt;
    
    Lp = zeros(3,3);
    for j=1:12
        Lp = gamma_dot(j) * m0(j,:)' * n0(j,:) + Lp; 
    end
    
    Fp   = expm(Lp * dt) * Fp;
    F = eye(3,3); 
    F(1,1) = lambda(ts+1);

    % Compute Fe
    Fe   = F / Fp;
    Ee   = 0.5 * (Fe'* Fe - eye(3,3));
    
    %PK2 @ t+dt in voigt form
    S_voigt = C * to_voigt_strain(Ee);
    
    % Convert to tensor form
    S = to_tensor(S_voigt);
    
    % Push forward PK2 to Cauchy Stress
    Sigma_3x3 = Fe * S * Fe' /det(Fe);
    
    % Now we can finally compute tau_alpha t+dt
    n = n0 / Fe;
    m = m0 * Fe';
    for j=1:12   
        tau(j) = sum(sum(Sigma_3x3.*(m(j,:)'*n(j,:))));
    end

    
    % Store current log strain @ t+dt
    B = F * F';
    e_log = - 0.5 * logm(inv(B));
    E(:,ts+1) = to_voigt_strain(e_log);
    
    % Store other parameters @ t+dt
    Sigma(:,ts+1) = to_voigt_stress(Sigma_3x3);  
    store.Sp(:,ts+1)  = Sp;
    store.tau(:,ts+1) = tau;
    store.gamma_dot(:,ts+1) = gamma_dot;
    store.gamma(:,ts+1) = gamma;
    
end
   
% Plot Stress vs strain
plot_stress_strain(E(1,:),Sigma)

% Plot shear Stress vs shear strain
plot_tau_gamma(abs(store.gamma),abs(store.tau))

% Plot shear strain rate
plot_gamma_dot(0:dt:t,store.gamma_dot)

% Flow shear stress
plot_Sp(0:dt:t,store.Sp)

end

function y = FlowRule(tau,Sp,const)

% Include Flow Rule Parameters
gamma_dot0 = const.gamma_dot0; E0 = const.E0; Sl = const.Sl;
m1 = const.m1 ; m2 = const.m2; kT = 4.11e-21;
%-------------------------------------------------------------------------%

y = gamma_dot0 * sign(tau) .* exp(-E0/(kT) * ...
  ((macaulay(1-(macaulay((abs(tau)-Sp)/Sl)).^m1)).^m2));

end

function hab = Hab_asaro(gamma)

%------------------------ Hardening parameters----------------------------%
q = 1.4;
tau0  = 16;
tau_s = 4.4 * tau0;    % Mpa
h0 = 8.25 * tau0;      % Mpa
hs = 0.48 * tau0;      % Mpa
%-------------------------------------------------------------------------%

A   = ones(3,3);
Qab = [ A     q.*A   q.*A   q.*A;
       q.*A     A    q.*A   q.*A;
       q.*A   q.*A     A    q.*A;
       q.*A   q.*A   q.*A     A];

gamma_bar = sum(abs(gamma));

h = hs + (h0 - hs) * (sech (gamma_bar*(h0 - hs)/(tau_s - tau0)))^2;
hab = h * Qab;

end

function hab = Hab_bassani(gamma)

%------------------------ Hardening parameters----------------------------%

tau0 = 1;
tau_I = 1.3 * tau0;
h0 = 90 * tau0;
hs_I = 4.45 * tau0;
hs_III = 0.15 * tau0;
gamma0 = 1d-3;
gamma0_III = 1.75 * tau0;
q = 0;
%-------------------------------------------------------------------------%
n = 8; h = 8; c = 8; g = 15 ; s =20;
V = [0,c,c,s,g,h,n,g,g,h,s,g,0,c,g,n,g,g,s,h,s,h,g,0,h,g,s,g,h,s,g,g,n,...
     0,c,c,g,n,g,g,s,h,0,c,s,g,h,g,h,s,0,h,g,s,n,g,g,0,c,c,h,g,s,0,...
     c,s,g,h,0,g,n,g,0,c,c,0,c,0];
 
fab = tril(ones(12)); 
fab(fab==1) = V;
fab = fab + triu(fab.',1);
%-------------------------------------------------------------------------%
hs = hs_I + (hs_III - hs_I) * tanh(sum(abs(gamma))./gamma0_III);
LHS  = hs + (h0-hs) * sech((h0-hs) / (tau_I-tau0) * abs(gamma)).^2;


fab_new = fab - diag(diag(fab));
RHS = fab_new * tanh(abs(gamma)./gamma0) + 1 ;
hab  = diag(LHS .* RHS);


for i =1:11
    hab(i+1:end,i) = q*hab(i,i);
end

for i = 12:2
    hab(1:i-1,i)   = q* hab(i,i) ;
end


end

function hab = Hab_luscher(gamma_dot,Sp)

%%------------------------ Hardening parameters---------------------------%
h0 = 200; %132        % MPa
r = 1.4;              % Unitless
A = 1.5e-19;          % Joules 
KT = 4.11e-21;        % Joules
Ss_0 = 205;           % MPa
gamma_dot0 = 1e7;     % 1/s
S0 = 1;               % MPa
%-------------------------------------------------------------------------%

Ss = Ss_0*(abs(gamma_dot ./ gamma_dot0)).^(KT/A);
II = h0*r*ones(12,12) + h0*(1-r)*eye(12,12);

hab = II * diag((Ss - Sp)./(Ss-S0));




end

function f = macaulay(x)

f = 0.5 * (x + abs(x));

end  

function xyz = rotate_matrix(m)

% Create a orthonrmal basis
x  = m ./norm(m);
yz = null(x)';
xyz = [x;yz];  


end

function C_new = rotate_c(R,const)

c11 = const.C11; c44 = const.C44; c12 = const.C12;
C = zeros(3,3,3,3);

% Rewrite nonzero values;
C(1,1,1,1) = c11;
C(2,2,2,2) = c11;
C(3,3,3,3) = c11;

C(2,3,2,3) = c44;
C(2,3,3,2) = c44;
C(3,2,2,3) = c44;
C(3,2,3,2) = c44;

C(3,1,3,1) = c44;
C(3,1,1,3) = c44;
C(1,3,3,1) = c44;
C(1,3,1,3) = c44;

C(1,2,1,2) = c44;
C(1,2,2,1) = c44;
C(2,1,1,2) = c44;
C(2,1,2,1) = c44;

C(1,1,2,2) = c12;
C(2,2,1,1) = c12;
C(1,1,3,3) = c12;
C(3,3,1,1) = c12;
C(2,2,3,3) = c12;
C(3,3,2,2) = c12;

% Rotate C to C' using 4 rotations
a = zeros(3,3,3,3);
for p=1:3
    for q=1:3
         for r=1:3
             for s=1:3 
                for i=1:3
                     for j=1:3
                         for k=1:3
                             for l=1:3
                                 
            a(p,q,r,s) = R(p,i)*R(q,j)*R(r,k)*R(s,l)*C(i,j,k,l)+a(p,q,r,s);
                                 
                             end
                         end
                     end
                end                
             end
         end
    end
end


% Change to voigt notation
C_new =[a(1,1,1,1) a(1,1,2,2) a(1,1,3,3) a(1,1,2,3) a(1,1,3,1) a(1,1,1,2);
        a(2,2,1,1) a(2,2,2,2) a(2,2,3,3) a(2,2,2,3) a(2,2,3,1) a(2,2,1,2);
        a(3,3,1,1) a(3,3,2,2) a(3,3,3,3) a(3,3,2,3) a(3,3,3,1) a(3,3,1,2);
        a(2,3,1,1) a(2,3,2,2) a(2,3,3,3) a(2,3,2,3) a(2,3,3,1) a(2,3,1,2);
        a(3,1,1,1) a(3,1,2,2) a(3,1,3,3) a(3,1,2,3) a(3,1,3,1) a(3,1,1,2);
        a(1,2,1,1) a(1,2,2,2) a(1,2,3,3) a(1,2,2,3) a(1,2,3,1) a(1,2,1,2)];
    
end

function[]  = plot_tau_gamma(x,y)

figure(1)
set(0,'DefaultAxesFontSize',18)
set(0,'defaultaxeslinewidth',4)
plot(x',y','LineWidth',4)
xlabel('{\gamma}^{\alpha}','fontsize',24,'FontWeight','bold',...
       'FontName','Times New Roman')
ylabel('{\tau}^{\alpha} (MPa)','fontsize',24,'FontWeight','bold',...
        'FontName','Times New Roman')
legend('B4','B5','B2','C1','C5','C3','D4','D1','D6',...
         'A3','A6','A2','Location','BestOutside') 
grid on


end     

function[]  = plot_stress_strain(ep,x)

sigma0 = 117; %MPa
x = x./sigma0 ; 
figure(2)
 set(0,'DefaultAxesFontSize',18)
 set(0,'defaultaxeslinewidth',4)
plot(ep,x(1,:),ep,x(2,:),ep,x(3,:),ep,x(4,:),...
                                   ep,x(5,:),ep,x(6,:),'LineWidth',4)
xlabel('ln({\epsilon}_{xx})','fontsize',24,'FontWeight','bold',...
       'FontName','Times New Roman')
ylabel('{\sigma}/{\sigma_0}','fontsize',24,'FontWeight','bold',...
       'FontName','Times New Roman')
legend('{\sigma}_{xx}/{\sigma_0}','{\sigma}_{yy}/{\sigma_0}',...
       '{\sigma}_{zz}/{\sigma_0}','{\sigma}_{yz}/{\sigma_0}',...
       '{\sigma}_{xz}/{\sigma_0}','{\sigma}_{xy}/{\sigma_0}',...
       'Location','BestOutside')  

grid on  

end

function[]  = plot_Sp(t,x)

 figure(3)
 set(0,'DefaultAxesFontSize',18)
 set(0,'defaultaxeslinewidth',4)
 plot(t,x','linewidth',4)
 xlabel('Time(s)','fontsize',24,'FontWeight','bold',...
        'FontName','Times New Roman')
 ylabel('\boldmath${S_{\rho}}^{\alpha}$ (MPa)','interpreter','latex',...
        'fontsize',24,'FontWeight','bold','FontName','Times New Roman')
 legend('B4','B5','B2','C1','C5','C3','D4','D1','D6',...
         'A3','A6','A2','Location','BestOutside') 
 grid on
end

function[]  = plot_gamma_dot(t,x)

 figure(4)
 set(0,'DefaultAxesFontSize',18)
 set(0,'defaultaxeslinewidth',4)
 plot(t,x','linewidth',4)
 xlabel('Time(s)','fontsize',24,'FontWeight','bold',...
        'FontName','Times New Roman')
 ylabel('\boldmath$\dot{\gamma_{\alpha}}$','interpreter','latex',...
        'fontsize',24,'FontWeight','bold','FontName','Times New Roman')
 legend('B4','B5','B2','C1','C5','C3','D4','D1','D6',...
     'A3','A6','A2','Location','BestOutside') 

 grid on
end

function v = to_voigt_stress(V)
    v = [V(1,1); V(2,2); V(3,3); V(2,3);V(1,3);V(1,2)];
end

function v = to_voigt_strain(V)
    v = [V(1,1); V(2,2); V(3,3); 2*V(2,3);2*V(1,3);2*V(1,2)];
end

function v = to_tensor(V)
    
    v  = [V(1) , V(6) , V(5);
          V(6) , V(2) , V(4);
          V(5) , V(4) , V(3)];
end
    
