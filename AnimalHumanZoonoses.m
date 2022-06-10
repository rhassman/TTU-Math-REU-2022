%%%%%%%%%%%%%%%%%%%%%%%%%% 
% CTMC HUman-Animal Zoonoses Epidemic Model

clear all
close all
set(0,'DefaultAxesFontSize', 18)
set(gca,'fontsize',18);

%Parameters and Initial Condition
Lambda_a = 0.5; beta_aa = 0.3; beta_ah = 0.2; beta_hh = 0.2;
mu_a = 1/360; gamma_a = 0.1; alpha_a = 0; alpha_h = 0.05;
gamma_h = 0.1; % units in days
% Initial conditions
S_h0 = 500; I_h0 = 0; R_h0 = 0; 
N_h0 = 500; N_a0 = Lambda_a / mu_a;
S_a0 = N_a0 - 2; I_a0 = 2; R_a0 = 0;

% Number of sampe paths, time step, and end time
nsim=10^5; 
dt=0.01; time=50;  outbreak=50;

count_ext=0; % count when i==0;
% CTMC Model (Gillespie's algorithm)
for k=1:nsim
clear t I_h I_a S_h S_a R_h R_a
% Initialize ve tors
t(1)=0; 
I_h(1)= I_h0; S_h(1) = S_h0; R_h(1) = R_h0;
I_a(1) = I_a0; S_a(1) = S_a0; R_a(1) = R_a0;

% Reset starting index
j=1;

% While loop 
while    I_h(j) + I_a(j) > 0 && I_h(j) + I_a(j) < outbreak %calculate  Dis. Ext. sim=10^5
%for d = 1:10
  %disp(I_h)
  u1 = rand;  u2 = rand; % two uniform random numbers
  N_h = S_h(j) + I_h(j) + R_h(j);
  N_a = S_a(j) + I_a(j) + R_a(j);

  % Interevent Probabilities
  prob_a_birth = Lambda_a;
  prob_a_infect = beta_aa*S_a(j)*I_a(j)/N_a;
  prob_a_recover = gamma_a*I_a(j);
  prob_sus_a_death = mu_a*S_a(j);
  prob_inf_a_death = (mu_a + alpha_a)*I_a(j);
  prob_rec_a_death = mu_a*R_a(j);
  prob_human_infect = beta_ah*S_h(j)*I_a(j)/N_h + beta_hh*S_h(j)*I_h(j)/N_h;
  prob_human_recover = gamma_h*I_h(j);
  prob_human_death = alpha_h*I_h(j);

  sum = prob_a_birth + prob_a_infect + prob_a_recover + prob_sus_a_death + ...
      prob_inf_a_death + prob_rec_a_death + prob_human_infect + ...
      prob_human_recover + prob_human_death;
  
  a_birth_bound = prob_a_birth / sum;
  a_infect_bound = prob_a_infect / sum + a_birth_bound;
  a_recover_bound = prob_a_recover / sum + a_infect_bound;
  sus_a_death_bound = prob_sus_a_death / sum + a_recover_bound;
  inf_a_death_bound = prob_inf_a_death / sum + sus_a_death_bound;
  rec_a_death_bound = prob_rec_a_death / sum + inf_a_death_bound;
  human_infect_bound = prob_human_infect / sum + rec_a_death_bound;
  human_recover_bound = prob_human_recover / sum + human_infect_bound;
  human_death_bound = prob_human_death / sum + human_recover_bound; % check = 1


  % interevent time
  t(j+1) =  t(j) -log(u1)/sum;

  % which event
  if u2 <= a_birth_bound 
      S_a(j + 1) = S_a(j) + 1;

      I_a(j + 1) = I_a(j);
      R_a(j + 1) = R_a(j);
      S_h(j + 1) = S_h(j);
      I_h(j + 1) = I_h(j);
      R_h(j + 1) = R_h(j);
  elseif u2 <= a_infect_bound && u2 > a_birth_bound
      S_a(j + 1) = S_a(j) - 1;
      I_a(j + 1) = I_a(j) + 1;
      R_a(j + 1) = R_a(j);
      S_h(j + 1) = S_h(j);
      I_h(j + 1) = I_h(j);
      R_h(j + 1) = R_h(j);
  elseif u2 <= a_recover_bound && u2 > a_infect_bound
      S_a(j + 1) = S_a(j);

      I_a(j + 1) = I_a(j) - 1;
      R_a(j + 1) = R_a(j) + 1;

      S_h(j + 1) = S_h(j);
      I_h(j + 1) = I_h(j);
      R_h(j + 1) = R_h(j);
  elseif u2 <= sus_a_death_bound && u2 > a_recover_bound
      S_a(j + 1) = S_a(j) - 1;

      I_a(j + 1) = I_a(j);
      R_a(j + 1) = R_a(j);
      S_h(j + 1) = S_h(j);
      I_h(j + 1) = I_h(j);
      R_h(j + 1) = R_h(j);
  elseif u2 <= inf_a_death_bound && u2 > sus_a_death_bound
      S_a(j + 1) = S_a(j);

      I_a(j + 1) = I_a(j) - 1;

      R_a(j + 1) = R_a(j);
      S_h(j + 1) = S_h(j);
      I_h(j + 1) = I_h(j);
      R_h(j + 1) = R_h(j);
  elseif u2 <= rec_a_death_bound && u2 > inf_a_death_bound
      S_a(j + 1) = S_a(j);
      I_a(j + 1) = I_a(j);

      R_a(j + 1) = R_a(j) - 1;

      S_h(j + 1) = S_h(j);
      I_h(j + 1) = I_h(j);
      R_h(j + 1) = R_h(j);
  elseif u2 <= human_infect_bound && u2 > rec_a_death_bound
      S_a(j + 1) = S_a(j);
      I_a(j + 1) = I_a(j);
      R_a(j + 1) = R_a(j);

      S_h(j + 1) = S_h(j) - 1;
      I_h(j + 1) = I_h(j) + 1;

      R_h(j + 1) = R_h(j);
  elseif u2 <= human_recover_bound && u2 > human_infect_bound
      S_a(j + 1) = S_a(j);
      I_a(j + 1) = I_a(j);
      R_a(j + 1) = R_a(j);
      S_h(j + 1) = S_h(j);

      I_h(j + 1) = I_h(j) - 1;
      R_h(j + 1) = R_h(j) + 1;
  elseif u2 <= human_death_bound && u2 > human_recover_bound
      S_a(j + 1) = S_a(j);
      I_a(j + 1) = I_a(j);
      R_a(j + 1) = R_a(j);
      S_h(j + 1) = S_h(j);

      I_h(j + 1) = I_h(j) - 1;

      R_h(j + 1) = R_h(j);
  
  end
 
  j=j+1;
end  

if k<=5

stairs(t,I_h,'color',rand(1,3),'Linewidth',2);
hold on
end

if I_h(j) + I_a(j) ==0
    count_ext =count_ext + 1;
end
end

%ODE Model (Euler's Method)
%y(1)=i0; x(1)=N0-i0; z(1)=0; N=x(1)+y(1)+z(1);
%for k=1:time/dt
  %f1=-beta*x(k)*y(k)/N;
  %f2=beta*x(k)*y(k)/N0-(gam+alp)*y(k);
  %f3=gam*y(k);
  %x(k+1)=x(k)+dt*f1;
  %y(k+1)=y(k)+dt*f2;
  %z(k+1)=z(k)+dt*f3;
  %N=x(k+1)+y(k+1)+z(k+1);
%end

probDisExt=count_ext/nsim %estimate for no outbreak, disease extinction
%probDisExt_est=(1/R0)^i0

% plot([0:dt:time],y,'k--','LineWidth',2);
% axis([0,time,0,220]);
title('CTMC SIR Animal Human Zoonoses')
ylabel('I(t)');
xlabel('Time (days)')
hold off