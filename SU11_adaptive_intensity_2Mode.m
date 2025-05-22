clear all;

Nloop = 30;
Nmode = Nloop+2;

r = 0.1;
%ph1 = pi/4;
th = 0;

Nph = 20;
%ph = zeros(Nph,1);
loop = zeros(Nloop+1,1);
In = zeros(Nloop+1,Nph);

c0 = [cosh(r) 0;0 cosh(r)];
s0 = [-exp(i*th)*sinh(r) 0;0 -exp(-i*th)*sinh(r)];

loop(1,1) = 0;
In(1,1) = 0;
for l = 0:Nph
    ph1 = (l+10)*pi/60;
    %ph1 = pi/4;
    
    u = [exp(i*ph1) 0;0 exp(-i*ph1)];
    du = [i*exp(i*ph1) 0;0 -i*exp(-i*ph1)];
    
    sig = eye(Nmode*2)/2;
    dsig = zeros(Nmode*2);

    for k = 0: Nloop
        Co = zeros(2*Nmode);
        Co = blkdiag(c0,eye(2*(k)),c0,eye(2*Nloop - 2*k));
     
%      S0 = zeros(2*Nmode);
%      S0 = blkdiag(s0,zeros(2*(k)),s0,zeros(2*Nloop - 2*k));
%      S1 = flip(S0);
        S1 = zeros(2*Nmode);
        S2 = zeros(2*Nmode);
    
        S1(1,2*(k+1)+2) = -exp(i*th)*sinh(r);
        S1(2,2*(k+1)+1) = -exp(-i*th)*sinh(r);
        S1(2*(k+1)+2,1) = -exp(-i*th)*sinh(r);
        S1(2*(k+1)+1,2) = -exp(i*th)*sinh(r);

        S2(1,2*(k+1)+2) = exp(i*th)*sinh(r);
        S2(2,2*(k+1)+1) = exp(-i*th)*sinh(r);
        S2(2*(k+1)+2,1) = exp(-i*th)*sinh(r);
        S2(2*(k+1)+1,2) = exp(i*th)*sinh(r);
        
        S = Co+S1;
        Sq = Co+S2;
        
        U = blkdiag(u,eye(2*Nmode-2));
        dU = blkdiag(du,zeros(2*Nmode-2));   
    
        dsig = (S*dU*S*sig*S'*U'*S' + S*U*S*dsig*S'*U'*S' + S*U*S*sig*S'*dU'*S');  
        sig1 = S*sig*S';
        sig = (S*U*S*sig*S'*U'*S');
    
    
        Isig = inv(sig);
        %H(k+1,1) = trace(Isig*dsig*Isig*dsig)/4; 

        In(k+1,l+1) = sig1(1,1) + sig1(2*(k+1)+1,2*(k+1)+1) - 1;
        %In(k+1,l) = sig(2*(k+1)+1,2*(k+1)+1)/2 - 1/2;
        %phasevar(k+1,1) = 1/sqrt(H(k+1,1));
        loop(k+1,1) = k;    
        ph(l+1,1) = l;
    end
end

Z = abs(In);
hold on;
% plot(ph,Z(10,:),'c')
% xlabel('\phi');
% ylabel('Intensity');
surf(ph,loop,Z)
%ribbon(abs(In))
