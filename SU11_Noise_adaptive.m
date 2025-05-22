clear all;
num = 100;
Nloop = 14;
Nmode = Nloop+2;

H = zeros(Nloop,num);
trans = zeros(num,1);
loop = zeros(Nloop,1);
In = zeros(Nloop,num);
Entropy = zeros(Nloop,num);

r = 0.1;
ph = pi/4;
th = 0;

c0 = [cosh(r) 0;0 cosh(r)];
s0 = [0 -exp(i*th)*sinh(r); -exp(-i*th)*sinh(r) 0];
u = [exp(i*ph) 0;0 exp(-i*ph)];
du = [i*exp(i*ph) 0;0 -i*exp(-i*ph)];

%for partial transpose
mode_a = eye(2*Nmode-2);
mode_b = [1 0;0 -1];
E = blkdiag(mode_b,mode_a);

%partially transposed symplectic matrix
o = [0 -1;1 0];
Om = o;
for j = 1:Nmode-1
    Om = blkdiag(Om,o);
end
Transposed_Om = E*Om*E;


transfer = 1;
for r1 = 0: num-1
    transfer = 1-0.001*(r1);
    sig = eye(Nmode*2)/2;
    dsig = zeros(Nmode*2);
    for k = 0: Nloop
        Co = blkdiag(c0,eye(2*(k)),c0,eye(2*Nloop - 2*k));
        S1 = zeros(2*Nmode);
    
        S1(1,2*(k+1)+2) = -exp(i*th)*sinh(r);
        S1(2,2*(k+1)+1) = -exp(-i*th)*sinh(r);
        S1(2*(k+1)+1,2) = -exp(-i*th)*sinh(r);
        S1(2*(k+1)+2,1) = -exp(i*th)*sinh(r);
        
        S = Co+S1;
    
        U = blkdiag(u,eye(2*Nmode-2));
        dU = blkdiag(du,zeros(2*Nmode-2));   
    
        dsig = (S*dU*S*sig*S'*U'*S'+S*U*S*sig*S'*dU'*S'+S*U*S*dsig*S'*U'*S');
        dsig = transfer*dsig;
     
        sig1 = S*sig*S';
        sig = (S*U*S*sig*S'*U'*S');
        sig = (transfer*(sig-eye(Nmode*2)/2)) + eye(Nmode*2)/2;      
        Isig = inv(sig);

        %partially transposed covaraince matrix
        Transposed_sig = E*sig*E;

        %symplectic transformation of covaraince matrix
        J = (Transposed_sig*Transposed_Om*Transposed_sig*Transposed_Om);
    
        [Ev,Eg] = eig(-J);

         S11=0;
        S22=0;
        for l = 1:2*Nmode
            S11 = log2(2*sqrt(Eg(l,l)));
            % if Eg(j,j)) > 0
            %     la = abs(Eg(j,j));
            %     S1 = S1+((la+1/2)*log2((la+1/2)) - (la-1/2)*log2((la-1/2)));
            % end
            if S11< 0
                S22=S22+S11;
            end
        end

        loop(k+1,1) = k+1;
        H(k+1,r1+1) = trace(Isig*dsig*Isig*dsig)/2;
        In(k+1,r1+1) = sig(1,1) -1/2;
        Entropy(k+1,r1+1) =-S22;
end
    trans(r1+1,1) =  1 - transfer;
end

hold on;
%plot(loop,H,'r')
surf(trans,loop,abs(Entropy))

