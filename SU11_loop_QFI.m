clear all
ntime = 15;


r = 0.1;
ph = pi/4;
th = 0;

time = zeros(ntime,1);
Qf = zeros(ntime,1);
QH = zeros(ntime,1);
In = zeros(ntime,1);

h = zeros(ntime,1);
nr = zeros(ntime,1);

%for k = 1: ntime       
%r = k*0.1;
   
sig = eye(4)/2;
dsig = zeros(4);
 
U = [1 0 0 0;0 1 0 0;0 0 exp(i*ph) 0;0 0 0 exp(-i*ph)];
S = [cosh(r) 0 0 -exp(i*th)*sinh(r); 0 cosh(r) -exp(-i*th)*sinh(r) 0; 0 -exp(-i*th)*sinh(r) cosh(r) 0; -exp(i*th)*sinh(r) 0 0 cosh(r)];
Sq = [cosh(r) 0 0 exp(i*th)*sinh(r); 0 cosh(r) exp(-i*th)*sinh(r) 0; 0 exp(-i*th)*sinh(r) cosh(r) 0; exp(i*th)*sinh(r) 0 0 cosh(r)];
dU = [0 0 0 0;0 0 0 0;0 0 i*exp(i*ph) 0;0 0 0 -i*exp(-i*ph)];
 
% Nsig = (Sq*U*S*S'*U'*Sq')/2;
% dNsig = (Sq*dU*S*S'*U'*Sq' + Sq*U*S*S'*dU'*Sq')/2;
% hi = trace(dNsig*inv(Nsig)*dNsig*inv(Nsig))/4;
% 
%  %h(k+1,1) = trace(dNsig*inv(Nsig)*dNsig*inv(Nsig))/4;
%  h(k,1) = 1/sqrt(hi)
%  nr(k,1) = k;
%  end
 
for k = 1:ntime  
    dsig = (S*dU*S*sig*S'*U'*S' + S*U*S*dsig*S'*U'*S' + S*U*S*sig*S'*dU'*S');
    sig1 = S*sig*S';
    sig = (S*U*S*sig*S'*U'*S');

    % dsig = (dU*S*sig*S'*U' + U*S*dsig*S'*U' + U*S*sig*S'*dU');
    % sig = (U*S*sig*S'*U');          
    % X = S*sig*S';
    % dX = S*dsig*S';
    % H = trace(dX*inv(X)*dX*inv(X))/2;

    H = trace(dsig*inv(sig)*dsig*inv(sig))/4; 

    time(k,1) = k;
    In(k,1) = (sig(1,1) - 1/2)^2;
    Qf(k,1) = H/((k)^2); 
    %%QH(k+1,1) = 1/sqrt(h*(k+1));
    QH(k,1) = 1/sqrt(H);

end

hold on;
plot(time,Qf,'m')
%plot(nr,h,'r')
xlabel('N');
ylabel('H(\phi)/N^2');

