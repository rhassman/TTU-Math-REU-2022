%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Royce CTMC Model
% %%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

% birth rates need to be multiplied by initial population

set(0,'DefaultAxesFontSize', 18)
set(gca,'fontsize',18);

% Initial Conditions
i0_w=1;i0_d=0; i0_h=0; t0_d=0; %initial infected
N0_w=1000;N0_d=1000; N0_h=1000; %Initial population
r0_w=0; r0_d=0; r0_h=0; %initial recovered
s0_w=N0_w-i0_w; s0_d=N0_d-i0_d-t0_d; s0_h=N0_h-i0_h; %initial susceptible

% Parameters
amp_bww=-0.9; beta_ww0 = 0.89;
 beta_wd=.51; beta_dd=.89; beta_td=.5; beta_th=.207; beta_hh=.078; %transmission rates
gamma_w=52/360; gamma_d=52/360; gamma_h=52/360; %recovery rates
mu_dt=.499; %mutation to human transmissable 
m_w=.0082; m_d=.0082; m_h=.009; %death rates
b_d=.0082*N0_d; b_h=N0_h*0.0118; %birth rates
dt=.00023844;


totext_sum=0; totext_fall=0; totext_winter=0; totext_spring=0;  % count when j==0;
% Number of sampe paths, time step, and end time
sim = 1000; outbreak=1000;
sea=4;

% CTMC Model (Gillespie's algorithm)
for g=1:sea 
TI=(365/4)*sea;

for k=1:sim
clear s_h s_w s_d i_w i_d i_h t_d r_h r_w r_d t N_h N_d N_w j
prealloc_size = 10;
s_h = zeros(prealloc_size, 1); s_w = zeros(prealloc_size,1); s_d = zeros(prealloc_size,1);
i_h = zeros(prealloc_size, 1); i_w = zeros(prealloc_size,1); i_d = zeros(prealloc_size,1); t_d = zeros(prealloc_size,1);
r_h = zeros(prealloc_size, 1); r_w = zeros(prealloc_size,1); r_d = zeros(prealloc_size,1);
t = zeros(prealloc_size, 1);

% Initialize vectors 
t(1)=0; i_w(1)=i0_w; i_d(1)=i0_d; i_h(1)=i0_h; t_d(1)=t0_d; 
s_w(1)=N0_w-i0_w; s_d(1)=N0_d; s_h(1)=N0_h;
r_w(1)=r0_w; r_d(1)=r0_d; r_h(1)=r0_h;
j=1;
N_w = 2; N_d = 2; N_h = 2; 

% While loop
while  i_h(j)+i_d(j)+i_w(j)+t_d(j) > 0 && i_h(j)+i_d(j)+i_w(j)+t_d(j) < outbreak && N_w > 1 && N_d > 1 && N_h > 1 
%while j < 5
  u2=rand; % Only need one random variable/ u2=rand; % two uniform random numbers
  
  beta_ww(j) = amp_bww*sin((t(j)+TI)*(3.14/180)) + beta_ww0;

  if sin((t(j)+TI)*(3.14/180))-.5>0
       b_w=.0025;
   elseif sin((t(j)+TI)*(3.14/180))-.5<=0
       b_w=0;
  end
count=0;
  N_w=s_w(j)+i_w(j)+r_w(j); % set population sizes
  N_d=s_d(j)+i_d(j)+r_d(j)+t_d(j);
  N_h=s_h(j)+i_h(j)+r_h(j);

  ev1=b_w*dt; % birth of wild bird
  ev2=m_w*s_w(j)*dt+ev1;% death of sw
  ev3=m_w*i_w(j)*dt+ev2;% death of Iw
  ev4=beta_ww(j)*s_w(j)*i_w(j)/(N_w)*dt+ev3; % infection of wild bird
  ev5=gamma_w*i_w(j)*dt+ev4;% recovery of wild bird
  ev6=m_w*r_w(j)*dt+ev5;% death of rw
  ev7=b_d*dt+ev6; %birth of chicken
  ev8=(beta_dd*s_d(j)*i_d(j)+beta_wd*s_d(j)*i_w(j))/(N_d)*dt+ev7;% infection of sd by id or iw
  ev9=beta_td*s_d(j)*t_d(j)/(N_d)*dt+ev8;% infection of sd by td
  ev10=gamma_d*i_d(j)*dt+ev9; % recovery of id
  ev11=gamma_d*t_d(j)*dt+ev10;% recovery of td
  ev12=m_d*s_d(j)*dt+ev11;% death of sd
  ev13=m_d*i_d(j)*dt+ev12; % death of id
  ev14=m_d*t_d(j)*dt+ev13;% death of td
  ev15=m_d*r_d(j)*dt+ev14;% death of rd
  ev16=mu_dt*i_d(j)*dt+ev15; % mutation of id
  ev17=b_h*dt+ev16;% birth of human
  ev18=(beta_th*s_h(j)*t_d(j)+beta_hh*s_h(j)*i_h(j))/(N_h)*dt+ev17;% infection of human
  ev19=m_h*s_h(j)*dt+ev18; %death of sh
  ev20=m_h*i_h(j)*dt+ev19;% death of ih
  ev21=m_h*r_h(j)*dt+ev20;% death of rh
  ev22=gamma_h*i_h(j)*dt+ev21;
  ev23=1;
   if ev22>1
       count=count+1;
   end

  t(j+1)=t(j)+dt; % interevent time

  if u2 <= ev1 % birth of wild bird
    %disp('birth of wild bird')
    s_w(j+1)=s_w(j)+1;
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev1 && u2<=ev2 %death of sw
    %disp('death of sw')
    s_w(j+1)=s_w(j)-1;
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev2 && u2<=ev3 %death of Iw
    %disp('death of Iw')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j)-1;
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev3 && u2<=ev4 %infection of wild bird
    %disp('infection of wild bird')
    s_w(j+1)=s_w(j)-1;
    i_w(j+1)=i_w(j)+1;
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev4 && u2<=ev5 %recovery of wild bird
    %disp('recovery of wild bird')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j)-1;
    r_w(j+1)=r_w(j)+1;
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev5 && u2<=ev6 %death of rw
    %disp('death of rw')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j)-1;
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev6 && u2<=ev7 %birth of chicken
    %disp('birth of chicken')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j)+1;
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev7 && u2<=ev8 %infection of sd by id or iw
    %disp('infection of sd by id or iw')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j)+1;
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev8 && u2<=ev9 %infection of sd by td
    %disp('infection of sd by td')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j)+1;
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev9 && u2<=ev10 %recovery of id
    %disp('recovery of id')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j)-1;
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j)+1;
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev10 && u2<=ev11 %recovery of td
    %disp('recovery of td')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j)-1;
    r_d(j+1)=r_d(j)+1;
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev11 && u2<=ev12 %death of sd
    %disp('death of sd')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j)-1;
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev12 && u2<=ev13 %death of id
    %disp('death of id')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j)-1;
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev13 && u2<=ev14 %death of td
    %disp('death of td')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j)-1;
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev14 && u2<=ev15 %death of rd
    %disp('death of rd')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j)-1;
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev15 && u2<=ev16 %mutation of id
    %disp('mutation of id')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j)-1;
    t_d(j+1)=t_d(j)+1;
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev16 && u2<=ev17 %birth of human
    %disp('birth of human')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j)+1;
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev17 && u2<=ev18 %infection of human
    %disp('infection of human')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j)+1;
    r_h(j+1)=r_h(j);
  elseif u2>ev18 && u2<=ev19 %death of sh
    %disp('death of sh')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j)-1;
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev19 && u2<=ev20 %death of ih
    %disp('death of ih')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j)-1;
    r_h(j+1)=r_h(j);
  elseif u2>ev20 && u2<=ev21 %death of rh
    %disp('death of rh')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j)-1;
  elseif u2>ev21 && u2<=ev22 %Recovery of Ih
    %disp('recovery of ih')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j)-1;
    r_h(j+1)=r_h(j)+1;
  elseif u2>ev22 && u2<=ev23 %No Event
    %disp('recovery of ih')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  end
  j=j+1;
  
end  
% resize vectors to remove zeros at end
s_h = s_h(1:j, 1); s_w = s_w(1:j, 1); s_d = s_d(1:j,1); 
i_h = i_h(1:j); i_w = i_w(1:j); i_d = i_d(1:j); t_d = t_d(1:j); 
r_h = r_h(1:j); r_w = r_w(1:j); r_d = r_d(1:j); 
t = t(1:j); 

figure(g)
if k==1
    disp(size(t))
    Y = [i_h i_w i_d];
    stairs(t, Y);
    hold on
end
if g==1
    
title('Stochastic Summer with Periodicity')
legend('ih', 'iw', 'id')
ylabel('Infections');
xlabel('Time (days)')
hold off

figure(g+4)
if k<=5
    disp(size(t))
    stairs(t,i_w,'color',rand(1,3),'LineWidth',2);
    hold on
end
title('Stochastic Summer with Periodicity')
ylabel('Infections of Wild Birds');
xlabel('Time (days)')
%hold off
elseif g==2
    title('Stochastic Fall with Periodicity')
legend('ih', 'iw', 'id')
ylabel('Infections');
xlabel('Time (days)')
hold off

figure(g+4)
if k<=5
    disp(size(t))
    stairs(t,i_w,'color',rand(1,3),'LineWidth',2);
    hold on
end
title('Stochastic Fall with Periodicity')
ylabel('Infections of Wild Birds');
xlabel('Time (days)')
elseif g==3
    title('Stochastic Winter with Periodicity')
legend('ih', 'iw', 'id')
ylabel('Infections');
xlabel('Time (days)')
hold off

figure(g+4)
if k<=5
    disp(size(t))
    stairs(t,i_w,'color',rand(1,3),'LineWidth',2);
    hold on
end
title('Stochastic Winter with Periodicity')
ylabel('Infections of Wild Birds');
xlabel('Time (days)')
elseif g==4
    title('Stochastic Spring with Periodicity')
legend('ih', 'iw', 'id')
ylabel('Infections');
xlabel('Time (days)')
hold off

figure(g+4)
if k<=5
    disp(size(t))
    stairs(t,i_w,'color',rand(1,3),'LineWidth',2);
    hold on
end
title('Stochastic Spring with Periodicity')
ylabel('Infections of Wild Birds');
xlabel('Time (days)')
end
if i_h(j)+i_w(j)+i_d(j)+t_d(j)==0 && g==1;
    totext_sum=totext_sum+1;
elseif i_h(j)+i_w(j)+i_d(j)+t_d(j)==0 && g==2;
    totext_fall=totext_fall+1;
elseif i_h(j)+i_w(j)+i_d(j)+t_d(j)==0 && g==3;
    totext_winter=totext_winter+1;
elseif i_h(j)+i_w(j)+i_d(j)+t_d(j)==0 && g==4;
    totext_spring=totext_spring+1;
end
end
end
Prob_ext_summer=totext_sum/sim
Prob_ext_fall=totext_fall/sim
Prob_ext_winter=totext_winter/sim
Prob_ext_spring=totext_spring/sim
