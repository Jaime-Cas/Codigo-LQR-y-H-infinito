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

G = ss(A,B,C,D);

% Respuesta escalon del sistema

figure(1)

subplot(221)
t=0:0.01:50;
[y,t] = step(G);
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

% Funciones de ponderacion
Ms = 0.8; wb = 0.8; epsi = 0.01;
We = tf([1/Ms wb],[1 wb*epsi]);

Mu = 30; wbc = 200; epsi1 = 0.01;
Wu = tf ([1 wbc/Mu],[epsi1 wbc]);

We = [We 0;0 We];
Wu = [Wu 0;0 Wu];

figure(2)
w = logspace(-4,6,1000);
sv = sigma(We,w);
sv1 = 20*log10(sv);
semilogx(w, sv1)
title('We Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure(3)
w = logspace(-4,6,1000);
sv = sigma(Wu,w);
sv2 = 20*log10(sv);
semilogx(w, sv2)
title('Wu Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

% Planta extendida P

P=augw(G,We,Wu,[])

% Compensador robuzto K
[K,CL,gamma] = hinfsyn(P); 
gamma

% Funciones de sensibilidad

L1 = G*K;
I = eye(size(L1));
S1 = feedback(I,L1);
T1 = I-S1;

figure(4)
w = logspace(-4,6,1000);
sv = sigma(S1,w);
sv = 20*log10(sv);
semilogx(w, sv,'r--',w ,sv1,'g')
title('S y Ws Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
legend('Ws','S')

figure(5)
w = logspace(-4,6,1000);
sv = sigma(T1,w);
sv = 20*log10(sv);
semilogx(w, sv,'r--',w ,sv2,'g')
title('T y Wt Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
legend('Wt','T')

% Respuesta temporal en lazo cerrado

figure(6)

subplot(221)
t=0:0.01:50;
[y,t] = step(T1);

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

%Respuesta en frecuencia

figure(7)
sigma(I+L1,'--',T1,':',L1,'r--',We/gamma,'k--',gamma/Wu,'k-',{.1,100})
legend('1/\sigma(S) performance','\sigma(T) robustness','\sigma(L) loopshape',...
'\sigma(We) performance bound','\sigma(1/Wu) robutness bound')