clear all;
%num = 10;
Nloop = 29;
Nmode = Nloop+2;

r = 0.14;
ph = pi/4;
th = 0;

H = zeros(Nloop+1,1);
loop = zeros(Nloop+1,1);
In = zeros(Nloop+1,1);
phasevar = zeros(Nloop+1,1);

c0 = [cosh(r) 0;0 cosh(r)];
s0 = [0 -exp(i*th)*sinh(r); -exp(-i*th)*sinh(r) 0];
u = [exp(i*ph) 0;0 exp(-i*ph)];
du = [i*exp(i*ph) 0;0 -i*exp(-i*ph)];

sig = eye(Nmode*2)/2;
dsig = zeros(Nmode*2);

    
for k = 0: Nloop
    
    Co = blkdiag(c0,eye(2*(k)),c0,eye(2*Nloop - 2*k));
    S1 = zeros(2*Nmode);
     
%     So = blkdiag(s0,zeros(2*Nmode-2));    
%     S1 = circshift(So,2*(k+1),2);
%     S2 = S1';
    
    S1(1,2*(k+1)+2) = -exp(i*th)*sinh(r);
    S1(2,2*(k+1)+1) = -exp(-i*th)*sinh(r);
    S1(2*(k+1)+1,2) = -exp(-i*th)*sinh(r);
    S1(2*(k+1)+2,1) = -exp(i*th)*sinh(r);
        
    S = Co+S1;
    
    U = blkdiag(u,eye(2*Nmode-2));
    dU = blkdiag(du,zeros(2*Nmode-2));   
    
    dsig = (S*dU*S*sig*S'*U'*S' + S*U*S*dsig*S'*U'*S' + S*U*S*sig*S'*dU'*S');  
    sig1 = S*sig*S';
    sig = (S*U*S*sig*S'*U'*S');  
    
    Isig = inv(sig);
    X = trace(Isig*dsig*Isig*dsig)/4

    H(k+1,1) = X/((k+1)^2); 
    %H(k+1,1) = X;

    In(k+1,1) = sig1(1,1)/2 - 1/2;
    phasevar(k+1,1) = 1/sqrt(H(k+1,1));
    loop(k+1,1) = k+1;    

%     dsig = U*dsig*U' + dU*sig*U' + U*sig*dU';
%     sig = U*sig*U';
end

hold on;
plot(loop,H,'r')
xlabel('No. of loops');
ylabel('Intensity');
