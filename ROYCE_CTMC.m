%%%%%%%%%%%%%%%%%%%%%%%%%% 
% SIR Epidemic Model
% with disease-relatd deaths
% Sample Paths for CTMC  models
% ODE model:
% ds/dt=-beta is/N;
% di/dt= beta is/N-(gam +alp)*i;
% dr/dt= gam*i;
% %%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
set(0,'DefaultAxesFontSize', 18)
set(gca,'fontsize',18);
%Parameters and Initial Condition
beta_ww=.5; beta_wd=.5; beta_dd=.5; beta_td=.5; beta_th=.5; beta_hh=.5; %transmission rates
gamma_w=.2; gamma_d=.2; gamma_h=.2; %recovery rates
mu_dt=.2; %mutation to human transmissable 
m_w=.3; m_d=.2; m_h=.2; %death rates
b_w=.5; b_d=.5; b_h=.5; %birth rates
i0_w=2;i0_d=0; i0_h=0; t0_d=0; %initial infected
N0_w=1000;N0_d=1000; N0_h=1000; %Initial population
r0_w=0; r0_d=0; r0_h=0; %initial recovered
s0_w=N0_w; s0_d=N0_d; s0_h=N0_h; %initial susceptible

R0=beta_hh/(gamma_h); % basic reprod number
totext=0; % count when i==0;
%Number of sampe paths, time step, and end time
sim=10^5; 
dt=0.01; time=50;  outbreak=50;
%sim=10^5; %use for calc of dis extinction
% ODE Model (Euler's Method)
y(1)=i0_h; x(1)=s0_h; z(1)=r0_h; Nh=x(1)+y(1)+z(1);
w(1)=s0_d; u(1)=i0_d; n(1)=t0_d; v(1)=r0_d; Nd=w(1)+u(1)+n(1)+v(1);
q(1)=s0_w; s(1)=i0_w; a(1)=r0_w; Nw=q(1)+s(1)+a(1);
for k=1:time/dt
  w1=b_w-beta_ww*q(k)*s(k)/Nw-m_w*q(k);
  w2=beta_ww*q(k)*s(k)/Nw-(gamma_w+m_w)*s(k);
  w3=gamma_w*s(k)-m_w*a(k);
  d1=b_d-beta_dd*w(k)*u(k)/Nd-beta_wd*w(k)*s(k)/Nd-beta_td*w(k)*n(k)/Nd-m_d*w(k);
  d2=beta_dd*w(k)*u(k)/Nd+beta_wd*w(k)*s(k)/Nd-mu_dt*u(k)-(gamma_d+m_d)*u(k);
  d3=mu_dt*u(k)+beta_td*w(k)*n(k)/Nd-(gamma_d+m_d)*n(k);
  d4=gamma_d*(u(k)+n(k))-m_d*v(k);
  h1=b_h-beta_hh*y(k)*x(k)/Nh-beta_th*n(k)*x(k)/Nh-m_h*x(k);
  h2=beta_hh*y(k)*x(k)/Nh+beta_th*n(k)*x(k)/Nh-(gamma_h+m_h)*y(k);
  h3=gamma_h*y(k)-m_h*z(k);
  q(k+1)=q(k)+dt*w1;
  s(k+1)=s(k)+dt*w2;
  a(k+1)=a(k)+dt*w3;
  w(k+1)=w(k)+dt*d1;
  u(k+1)=u(k)+dt*d2;
  n(k+1)=n(k)+dt*d3;
  v(k+1)=v(k)+dt*d4;
  x(k+1)=x(k)+dt*h1;
  y(k+1)=y(k)+dt*h2;
  z(k+1)=z(k)+dt*h3;
  Nh=x(k+1)+y(k+1)+z(k+1);
  Nw=q(k+1)+s(k+1)+a(k+1);
  Nd=w(k+1)+u(k+1)+n(k+1)+v(k+1);
end

%CTMC Model (Gillespie's algorithm)
for k=1:sim
clear s_h s_w s_d i_w i_d i_h t_d r_h r_w r_d t N_h N_d N_w
t(1)=0; i_w(1)=i0_w; i_d(1)=i0_d; i_h(1)=i0_h; t_d(1)=t0_d; 
s_w(1)=N0_w-i0_w; s_d(1)=N0_d; s_h(1)=N0_h;
r_w(1)=r0_w; r_d(1)=r0_d; r_h(1)=r0_h;
j=1;
while  i_h(j)+i_d(j)+i_w(j)+t_d(j)>0 && i_h(j)<outbreak % graph sample paths sim=5
%while    i(j)>0 && i(j)<outbreak %calculate  Dis. Ext. sim=10^5
  u1=rand;  u2=rand; % two uniform random numbers
  N_w=s_w(j)+i_w(j)+r_w(j);
  N_d=s_d(j)+i_d(j)+r_d(j)+t_d(j);
  N_h=s_h(j)+i_h(j)+r_h(j);
  sum=b_w+m_w*s_w(j)+((beta_ww*s_w(j)*i_w(j))/N_w)+gamma_w*i_w(j)+m_w*i_w(j)+m_w*r_w(j)+b_d+((beta_dd*s_d(j)*i_d(j))/N_d)+((beta_td*s_d(j)*t_d(j))/N_d)+((beta_wd*s_d(j)+i_w(j))/N_d)+gamma_d*(i_d(j)+t_d(j))+m_d*(s_d(j)+i_d(j)+r_d(j))+mu_dt*i_d(j)+b_h+((beta_th*s_h(j)*t_d(j))/N_h)+((beta_hh*s_h(j)*i_h(j))/N_h)+m_h*(s_h(j)+i_h(j)+r_h(j))+gamma_h*i_h(j); %sum of all rates
  ev1=b_w/sum; % birth of wild bird
  ev2=m_w*s_w(j)/sum+ev1;% death of sw
  ev3=m_w*i_w(j)/sum+ev2;% death of Iw
  ev4=(beta_ww*s_w(j)*i_w(j))/(N_w*sum)+ev3; % infection of wild bird
  ev5=gamma_w*i_w(j)/sum+ev4;% recovery of wild bird
  ev6=m_w*r_w(j)/sum+ev5;% death of rw
  ev7=b_d/sum+ev6; %birth of chicken
  ev8=(beta_dd*s_d(j)*i_d(j)+beta_wd*s_d(j)+i_w(j))/(N_d*sum)+ev7;% infection of sd by id or iw
  ev9=((beta_td*s_d(j)*t_d(j))/N_d*sum)+ev8;% infection of sd by td
  ev10=(gamma_d*i_d(j))/sum+ev9; % recovery of id
  ev11=(gamma_d*t_d(j))/sum+ev10;% recovery of td
  ev12=m_d*s_d(j)/sum+ev11;% death of sd
  ev13=m_d*i_d(j)/sum+ev12; % death of id
  ev14=m_d*t_d(j)/sum+ev13;% death of td
  ev15=m_d*r_d(j)/sum+ev14;% death of rd
  ev16=mu_dt*i_d(j)/sum+ev15; % mutation of id
  ev17=b_h/sum+ev16;% birth of human
  ev18=(beta_th*s_h(j)*t_d(j)+beta_hh*s_h(j)+i_h(j))/(N_h*sum)+ev17;% infection of human
  ev19=m_h*s_h(j)/sum+ev18; %death of sh
  ev20=m_h*i_h(j)/sum+ev19;% death of ih
  ev21=m_h*r_h(j)/sum+ev20;% death of rh
  ev22=(gamma_h*i_h(j))/sum+ev21; % recovery of ih
  %check ev2=1;
  t(j+1)=t(j)-log(u1)/sum;
  if u2 <= ev1; % birth of wild bird
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
  elseif u2>ev1 && u2<=ev2; %death of sw
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
  elseif u2>ev2 && u2<=ev3; %death of Iw
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
  elseif u2>ev3 && u2<=ev4; %infection of wild bird
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
  elseif u2>ev4 && u2<=ev5; %recovery of wild bird
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
  elseif u2>ev5 && u2<=ev6; %death of rw
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
  elseif u2>ev6 && u2<=ev7; %birth of chicken
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
  elseif u2>ev7 && u2<=ev8; %infection of sd by id or iw
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
  elseif u2>ev8 && u2<=ev9; %infection of sd by td
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
    elseif u2>ev9 && u2<=ev10; %recovery of id
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
  elseif u2>ev10 && u2<=ev11; %recovery of td
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
  elseif u2>ev11 && u2<=ev12; %death of sd
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
    elseif u2>ev12 && u2<=ev13; %death of id
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
  elseif u2>ev13 && u2<=ev14; %death of td
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
  elseif u2>ev14 && u2<=ev15; %death of rd
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
     elseif u2>ev15 && u2<=ev16; %mutation of id
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
  elseif u2>ev16 && u2<=ev17; %birth of human
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
  elseif u2>ev17 && u2<=ev18; %infection of human
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
     elseif u2>ev18 && u2<=ev19; %death of sh
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
  elseif u2>ev19 && u2<=ev20; %death of ih
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
  elseif u2>ev20 && u2<=ev21; %death of rh
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
  elseif u2>ev21 && u2<=ev22; %Recovery of Ih
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
  end
  j=j+1;
end  
if k<=5
stairs(t,i_h,'color',rand(1,3),'Linewidth',2);
hold on
end
if i_h(j)+i_w(j)+i_d(j)==0
    totext=totext+1;
end
end
%probDisExt=totext/sim %estimate for no outbreak, disease extinction
probDisExt_est=(1/R0)^i0_w
plot([0:dt:time],y,'k--','LineWidth',2);
axis([0,time,0,220]);
title('CTMC and ODE SIR with disease-related deaths')
ylabel('I(t)');
xlabel('Time (days)')
hold off