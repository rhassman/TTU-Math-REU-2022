%%%%%%%%%%%%%%%%%%%%%%%%%% 
% CTMC Model Environmental Transmission in Wild Birds Full Stochastic
% %%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
dbstop if error
% birth rates need to be multiplied by initial population

set(0,'DefaultAxesFontSize', 18);
set(gca,'fontsize',18);

% Initial Conditions
i0_w=1;i0_d=0; i0_h=0; t0_d=0; %initial infected
N0_w=1000;N0_d=1000; N0_h=1000; %Initial population
r0_w=0; r0_d=0; r0_h=0; %initial recovered
s0_w=N0_w-i0_w; s0_d=N0_d-i0_d-t0_d; s0_h=N0_h-i0_h; %initial susceptible


R0w=5; R0d=5; R0t=5; R0h=5;
gamma_w=1/11; gamma_d=1/11; gamma_h=1/11; %recovery rates
mu_dt=.499; %mutation to human transmissible 
m_w=1/(3*360); m_d=1/360; m_h=1/(70*360); %death rates
beta_ww0 = R0w*(m_w+gamma_w); amp_bww=0.25*beta_ww0; %seasonal transmission
beta_dd=R0d*(m_w+gamma_w+mu_dt); beta_td=R0t*(m_d+gamma_d); beta_hh=R0h*(m_h+gamma_h); %R0 transmission rates
beta_wd=0.25*beta_dd;  beta_th=0.5*beta_td; %contact based transmission
b_w=1/(3*360)*N0_w; b_d=m_d*N0_d; b_h=N0_h*m_h; %birth rates
dt=10e-4;



%%% Environmental transmissibility params pulled from Breban
% (eta just inspired by Breban)
eta = [3 5 4 2]; % Decay rate spring summer fall winter
alpha = 1; % 'Environmental infectiousness' % Maybe just here for units?
rho = 10e-4; % Exposure rate
rho_w = rho; rho_d = rho; rho_h = rho;
omega = 10e5; % virus shedding rate
omega_w = omega; omega_d = omega; omega_t = omega;
eta_w = eta; eta_d = eta; eta_t = eta;
dose = 10e5; % number of virions for infectious dose

  % count when j==0;
% Number of sampe paths, time step, and end time
sim = 1; outbreak=100;
sea=4;

% CTMC Model (Gillespie's algorithm)

for g=1:sea 
%g = 4;

TI=90*(g - 1) + 60;
totext(g)=0;
clear x

for k=1:sim
clear s_h s_w s_d i_w i_d i_h t_d r_h r_w r_d t N_h N_d N_w j x v_w 

prealloc_size = 52000;
s_h = zeros(4, prealloc_size); s_w = zeros(4,prealloc_size); s_d = zeros(4,prealloc_size);
i_h = zeros(4,prealloc_size); i_w = zeros(4,prealloc_size); i_d = zeros(4,prealloc_size); t_d = zeros(4,prealloc_size);
r_h = zeros(4,prealloc_size); r_w = zeros(4,prealloc_size); r_d = zeros(4,prealloc_size);
v_w = zeros(4,prealloc_size); 
beta_ww = zeros(1, prealloc_size); 
t = zeros(1,prealloc_size);

% Initialize vectors 
t(1)=0; i_w(g,1)=i0_w; i_d(g,1)=i0_d; i_h(g,1)=i0_h; t_d(g,1)=t0_d; 
s_w(g,1)=N0_w-i0_w; s_d(g,1)=N0_d; s_h(g,1)=N0_h;
r_w(g,1)=r0_w; r_d(g,1)=r0_d; r_h(g,1)=r0_h;
j=1;
N_w = 2; N_d = 2; N_h = 2; 

% While loop
while  i_h(g,j)+i_d(g,j)+i_w(g,j)+t_d(g,j) > 0 && i_h(g,j)+i_d(g,j)+i_w(g,j)+t_d(g,j) < outbreak && N_w > 1 && N_d > 1 && N_h > 1 
%while j < 100000
  u2=rand; % Only need one random variable/ u2=rand; % two uniform random numbers

  beta_ww(j) = amp_bww*sin((t(j) + TI)*3.14/180 - 3.14) + beta_ww0;
  b_w = (0.5*tanh(57*sin(((t(j) + TI)*pi/180)-106*2*pi/360)-40)+0.5)/270;

count=0;
  N_w=s_w(g,j)+i_w(g,j)+r_w(g,j); % set population sizes
  N_d=s_d(g,j)+i_d(g,j)+r_d(g,j)+t_d(g,j);
  N_h=s_h(g,j)+i_h(g,j)+r_h(g,j);

  ev1=b_w*dt; % birth of wild bird
  ev2=m_w*s_w(g,j)*dt+ev1;% death of sw
  ev3=m_w*i_w(g,j)*dt+ev2;% death of Iw
  ev4=dt*(beta_ww(j)*s_w(g,j)*i_w(g,j)/N_w + rho_w*s_w(g, j)*v_w(g, j)/dose)+ev3; % infection of wild bird
  ev5=gamma_w*i_w(g,j)*dt+ev4;% recovery of wild bird
  ev6=m_w*r_w(g,j)*dt+ev5;% death of rw
  ev7=b_d*dt+ev6; %birth of chicken
  ev8=((beta_dd*s_d(g,j)*i_d(g,j)+beta_wd*s_d(g,j)*i_w(g,j))/N_d)*dt+ev7;% infection of sd by id or iw
  ev9=(beta_td*s_d(g,j)*t_d(g,j)/N_d)*dt +ev8;% infection of sd by td
  ev10=gamma_d*i_d(g,j)*dt+ev9; % recovery of id
  ev11=gamma_d*t_d(g,j)*dt+ev10;% recovery of td
  ev12=m_d*s_d(g,j)*dt+ev11;% death of sd
  ev13=m_d*i_d(g,j)*dt+ev12; % death of id
  ev14=m_d*t_d(g,j)*dt+ev13;% death of td
  ev15=m_d*r_d(g,j)*dt+ev14;% death of rd
  ev16=mu_dt*i_d(g,j)*dt+ev15; % mutation of id
  ev17=b_h*dt+ev16;% birth of human
  ev18=((beta_th*s_h(g,j)*t_d(g,j)+beta_hh*s_h(g,j)*i_h(g,j))/N_h)*dt +ev17;% infection of human
  ev19=m_h*s_h(g,j)*dt+ev18; %death of sh
  ev20=m_h*i_h(g,j)*dt+ev19;% death of ih
  ev21=m_h*r_h(g,j)*dt+ev20;% death of rh
  ev22=gamma_h*i_h(g,j)*dt+ev21; % recovery of human
  ev23=omega_w*i_w(g, j)*dt/dose+ev22; % shed infectious dose v_w
  ev24=eta_w(g)*v_w(g, j)*dt/dose+ev23; % decay infectious dose v_w
  ev25=1; % nothing happens
   
  t(j+1)=t(j)+dt; % interevent time

  if u2 <= ev1 % birth of wild bird
    %disp('birth of wild bird')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j)+1;
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev2  %death of sw
    %disp('death of sw')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j)-1;
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev3 %death of Iw
    %disp('death of Iw')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j)-1;
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev4  %infection of wild bird
    %disp('infection of wild bird')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j)-1;
    i_w(g,j+1)=i_w(g,j)+1;
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev5  %recovery of wild bird
    %disp('recovery of wild bird')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j)-1;
    r_w(g,j+1)=r_w(g,j)+1;
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev6  %death of rw
    %disp('death of rw')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j)-1;
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev7  %birth of chicken
    %disp('birth of chicken')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j)+1;
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev8  %infection of sd by id or iw
    %disp('infection of sd by id or iw')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j)-1;
    i_d(g,j+1)=i_d(g,j)+1;
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev9  %infection of sd by td
    %disp('infection of sd by td')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j)-1;
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j)+1;
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev10  %recovery of id
    %disp('recovery of id')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j)-1;
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j)+1;
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev11  %recovery of td
    %disp('recovery of td')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j)-1;
    r_d(g,j+1)=r_d(g,j)+1;
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev12 %death of sd
    %disp('death of sd')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j)-1;
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev13  %death of id
    %disp('death of id')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j)-1;
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev14  %death of td
    %disp('death of td')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j)-1;
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev15  %death of rd
    %disp('death of rd')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j)-1;
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev16  %mutation of id
    %disp('mutation of id')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j)-1;
    t_d(g,j+1)=t_d(g,j)+1;
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev17  %birth of human
    %disp('birth of human')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j)+1;
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev18  %infection of human
    %disp('infection of human')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j)-1;
    i_h(g,j+1)=i_h(g,j)+1;
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev19   %death of sh
    %disp('death of sh')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j)-1;
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev20 %death of ih
    %disp('death of ih')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j)-1;
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev21  %death of rh
    %disp('death of rh')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j)-1;
  elseif u2<=ev22  %Recovery of Ih
    %disp('recovery of ih')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j)-1;
    r_h(g,j+1)=r_h(g,j)+1;
  elseif u2<=ev23 % shed infectious dose v_w
    v_w(g,j+1)=v_w(g,j)+1;
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  elseif u2<=ev24 % decay infectious dose v_w
    v_w(g,j+1)=v_w(g,j)-1;
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);

  else  %No Event
    %disp('no event')
    v_w(g,j+1)=v_w(g,j);
    s_w(g,j+1)=s_w(g,j);
    i_w(g,j+1)=i_w(g,j);
    r_w(g,j+1)=r_w(g,j);
    s_d(g,j+1)=s_d(g,j);
    i_d(g,j+1)=i_d(g,j);
    t_d(g,j+1)=t_d(g,j);
    r_d(g,j+1)=r_d(g,j);
    s_h(g,j+1)=s_h(g,j);
    i_h(g,j+1)=i_h(g,j);
    r_h(g,j+1)=r_h(g,j);
  end
  j=j+1;
  
end  

% resize vectors to remove zeros at end
s_h = s_h(g, 1:j); s_w = s_w(g, 1:j); s_d = s_d(g, 1:j); 
i_h = i_h(g, 1:j); i_w = i_w(g, 1:j); i_d = i_d(g, 1:j); t_d = t_d(g, 1:j); 
r_h = r_h(g, 1:j); r_w = r_w(g, 1:j); r_d = r_d(g, 1:j); 
v_w = v_w(g, 1:j); 
t = t(1, 1:j) + 90*(g - 1) + 60;


figure(g)
if k<=5
    %col_list = ["A2142F" "#EDB120" "#77AC30" "#0072BD"; ...
                %"#D95319" "yellow" "green" "#4DBEEE"];
    col_list = ["red" "green" "cyan" "blue" "magenta"];
    
    col1 = convertStringsToChars(col_list(k)); 
    plot(t, v_w, 'red', t, i_w, 'blue')
    legend('v_w(t)','i_w(t)')
    %plot(t, i_w, col1)
    hold on

end

if i_h(j)+i_w(j)+i_d(j)+t_d(j)==0
    totext(g)=totext(g)+1;
end
end
end
figure(1)
title('CTMC Spring with Periodicity')
ylabel('I_d(t) and v_d(t)');
xlabel('Time (days)')
hold off

figure(2)
title('CTMC Summer with Periodicity')
ylabel('I_d(t) and v_d(t)');
xlabel('Time (days)')
hold off

figure(3)
title('CTMC Fall with Periodicity')
ylabel('I_d(t) and v_d(t)');
xlabel('Time (days)')
hold off

figure(4)
title('CTMC Winter with Periodicity')
ylabel('I_d(t) and v_d(t)');
xlabel('Time (days)')
hold off

sim_prob_ext_summer=totext(1)/sim
sim_prob_ext_fall=totext(2)/sim
sim_prob_ext_winter=totext(3)/sim
sim_prob_ext_spring=totext(4)/sim

R0h=beta_hh/(gamma_h+m_h)
R0w=beta_ww0/(gamma_w+m_w)
R0t=beta_td/(gamma_d+m_d)
R0d=beta_dd/(gamma_d+m_d+mu_dt)
