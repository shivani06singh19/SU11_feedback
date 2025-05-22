clear all
ntime = 30;

r = 0.1;
%ph = pi/4;
th = 0;
time = zeros(ntime,1);
an = zeros(20,1);
In = zeros(ntime,20);
%time(1,1) = 0;
%In(1,1) = 0;

for k1 = 0:20
ph = (k1+10)*pi/40;
% ph = pi/4; 
sig = eye(4)/2;

U = [1 0 0 0;0 1 0 0;0 0 exp(i*ph) 0;0 0 0 exp(-i*ph)];
S = [cosh(r) 0 0 -exp(i*th)*sinh(r); 0 cosh(r) -exp(-i*th)*sinh(r) 0;0 -exp(-i*th)*sinh(r) cosh(r) 0; -exp(i*th)*sinh(r) 0 0 cosh(r)];
S1 = [cosh(r) 0 0 exp(i*th)*sinh(r); 0 cosh(r) exp(-i*th)*sinh(r) 0;0 exp(-i*th)*sinh(r) cosh(r) 0; exp(i*th)*sinh(r) 0 0 cosh(r)];


for k = 0:ntime  
    sig1 = S*sig*S';
    sig = (S*U*S*sig*S'*U'*S');
                  
    time(k+1,1) = k;
    In(k+1,k1+1) = sig1(1,1) + sig1(3,3) - 1  ;
    an(k1+1,1) = k1;
end
 end

 Y= time;
 X = an;
 Z = abs(In);

hold on;
%plot(time,In,'r')
% xlabel('No. of loops');
% ylabel('H(\phi)');
%ribbon(Z)
surf(X,Y,Z)
