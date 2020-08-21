clc
clear

pkg load control
pkg load signal

m1 = 10 %%kg
m2 = 30 %%kg
k1 = 2.5 %% KN/m
k2 = 1.5 %%KN/m
b1 = 30 %% Ns/m

A = [0 1 0 0; 
-k1/m1 -b1/m1 0 b1/m1;
 0 0 0 1; 
 0 b1/m2 -k2/m2 -b1/m2];
 
B = [0 0;1/m1 0;
0 0;0 1/m2];

C = [1 0 0 0;
0 0 1 0];

D = zeros(size(C,1),size(B,2));

sys = ss(A,B,C,D);

%%Matriz de transmision de la planta
Gs = tf(sys);

%%Polos y zeros de transmision
figure(1)
pzmap(sys);
title('Polos y Zeros de transmision')

%%Barreras de estabilidad y desempeño

w=logspace(-2,2,400);
a=80;b=120;c=200;
x1=25*ones(1,a); x2=60*zeros(1,b); x3=-25*ones(1,c);
xt=[x1 x2 x3];
x4=[5*ones(1,100) 0*zeros(1,300)];

semilogx(w, xt,'r')
ylim([-50,50])
grid;
hold on;
semilogx(w,x4,'b')
title('Barreras de Estabilidad')
xlabel('w(rad/s)')
ylabel('dB')

%%Hallamos los valores singular y las ploteamos

w=logspace(-2,2,100);
[sv,w]=sigma(Gs,w);
sv=20*log10(sv);
figure(3)
semilogx(w,sv)
xlabel('Frecuencia (rad/sec)')
ylabel('Amplitud (dB)')
title('Valores singulares de la planta aumentada con I/s')
grid

% Planta aumentada por integrador
[ns nc] = size(B); % ns = number of inputs; nc = number of controls;
Ai = [ A B
0*ones(nc,ns) 0*ones(nc,nc)]
Bi = [ 0*ones(ns,nc)
eye(nc) ]
Ci = [ C 0*ones(nc,nc)]
Di = 0*ones(nc,nc)

figure(4);
w = logspace(-2,2,100);
sv = sigma(ss(Ai, Bi, Ci, Di),w);
sv = 20*log10(sv);
semilogx(w, sv);

xlabel('Frecuencia (rad/sec)')
ylabel('Amplitud (dB)')
title('Valores singulares de la planta aumentada con I/s')
grid

%Controlador LQR
Q=Ci'*Ci;
R=0.1*eye(2);
[K, P, E] = lqr(Ai,Bi,Q,R)

%Filtro de Kalman
sys = ss(Ai,Bi,Ci,Di);
R=0.1*eye(2);
Q=1*eye(6);
[H,P,E] = lqe(sys,Q,R)

figure(5)
w = logspace(-2,2,100);
sv = sigma(ss(Ai,H,Ci,Di),w);
sv = 20*log10(sv);
semilogx(w,sv);
xlabel('Frecuencia (rad/sec)')
ylabel('Amplitud (dB)')
title('Malla del filtro (G_{KF} = C*inv(sI-A)*H)')
grid

% Estructura del compensador usando G=K y H

Aq = [Ai-Bi*K-H*Ci 0*ones(ns+nc,nc)
          K       0*ones(nc,nc)]
Bq = [H 
      0*ones(nc,nc)]
Cq = [0*ones(nc,ns+nc) eye(nc,nc)]
Dq = 0*ones(nc,nc)

[cpoles, czeros] = pzmap(ss(Ai, H, K, Di));
[polecheck, zerocheck] = pzmap(ss(Aq, Bq, Cq, Dq));

figure(6)
w = logspace(-2,2,100);
sv = sigma(ss(Aq, Bq, Cq, Dq),w);
sv = 20*log10(sv);
semilogx(w,sv)
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
title('Compensator Singular Values')
grid

% Analisis en lazo abierto

Al = [ A                     B*Cq
       0*ones(ns+nc+nc,ns)    Aq    ]

Bl = [ 0*ones(ns,nc)
       Bq ]
    
Cl = [ C  0*ones(nc,ns+nc+nc) ]

[olpoles, olzeros] = pzmap(ss(Al,Bl,Cl,0*ones(nc,nc)));

figure(7)
    
w = logspace(-2,2,100);
sv = sigma(ss(Al, Bl, Cl, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Open Loop Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%Sensibilidad (S)

figure(8)
sv = sigma(ss(Al-Bl*Cl, Bl, -Cl, eye(nc)),w);
sv1 = 20*log10(sv);
semilogx(w, sv1)
title('Sensitivity Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%Sensibilidad complementaria (T)

figure(9)
sv = sigma(ss(Al-Bl*Cl, Bl, Cl, 0*eye(nc)),w);
sv2 = 20*log10(sv);
semilogx(w, sv2)
title('Complementary Sensitivity Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%Graficas de S y T superpuestas

figure(10)
semilogx(w, sv1,w,sv2)
title('Complementary Sensitivity and Sensitivity Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

% Respuesta a una entrada escalon en lazo cerrado

figure(11)

subplot(221)
t=0:0.01:50;
[y,t] = step(ss(Al-Bl*Cl, Bl, Cl, 0*eye(nc)));
plot(t,y(:,1,1))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 for input 1')

subplot(222)
plot(t,y(:,1,2))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 for input 2')

subplot(223)
plot(t,y(:,2,1))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 2 for input 1')

subplot(224)
plot(t,y(:,2,2))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 2 for input 2')