clear all;
nloop = 15;
num = 100;

r = 0.1;
ph = pi/4;

H = zeros(nloop,num);
H1 = zeros(nloop,num);
loop = zeros(nloop,1);
trans_vec = zeros(num,1);
In = zeros(nloop,num);
ent = zeros(nloop,num);

U = [exp(i*ph) 0 0 0;0 exp(-i*ph) 0 0;0 0 1 0;0 0 0 1];
S = [cosh(r) 0 0 -sinh(r);0 cosh(r) -sinh(r) 0;0 -sinh(r) cosh(r) 0;-sinh(r) 0 0 cosh(r)];
dU = [i*exp(i*ph) 0 0 0;0 -i*exp(-i*ph) 0 0;0 0 0 0;0 0 0 0];
Dis = [1; 1; 1; 1]/2;

%for partial transpose
mode_a = eye(2);
mode_b = [1 0;0 -1];
E = blkdiag(mode_b,mode_a);

%partially transposed symplectic matrix
om = [0 1;-1 0];
Om= blkdiag(om,om);
Transposed_Om = E*Om*E;

%transfer = rand(nloop,1);
transfer = 1
for j = 0: num-1
    transfer = 1-0.001*(j);
    sig = eye(4)/2;
    dsig = zeros(4);
for k = 1:nloop
    dsig = (S*dU*S*sig*S'*U'*S'+S*U*S*sig*S'*dU'*S'+S*U*S*dsig*S'*U'*S');
    dsig = transfer*dsig;

    sig1 = S*sig*S';
    sig = (S*U*S*sig*S'*U'*S');
    sig = (transfer*(sig-eye(4)/2)) + eye(4)/2;
    isig = inv(sig);

    % dsig = (dU*S*sig*S'*U'+U*S*sig*S'*dU'+U*S*dsig*S'*U');
    % dX = S*dsig*S';
    % dX = transfer*dX;
    % 
    % sig = (U*S*sig*S'*U');
    % X = S*sig*S';
    % X = (transfer*(X-eye(4)/2)) + eye(4)/2;
    % iX = inv(X);


    %partially transposed covaraince matrix
    Transposed_sig = E*sig*E;

    %symplectic transformation of covaraince matrix
    J = (Transposed_sig*Transposed_Om*Transposed_sig*Transposed_Om);
    
    [Ev,Eg] = eig(-J);

    S1=0;
    S2=0;
    for l = 1:4
        S1 = log2(2*sqrt(Eg(l,l)));
        if S1< 0
            S2=S2+S1;
        end
    end

    loop(k,1) = k;
    H(k,j+1) = trace(isig*dsig*isig*dsig)/2;
    % H1(k,j+1) = trace(iX*dX*iX*dX)/2;
    In(k,j+1) = sig1(1,1) -1/2;
    ent(k,j+1) =-S2; %divided by the number of modes = 4
end
    trans_vec(j+1,1) = 1 - transfer;
end

hold on
surf(trans_vec,loop,abs(H))
%surf(trans_vec,loop,abs(ent))
%plot(trans_vec,H,'b')
%ribbon(abs(H))
% xlabel('Loss');
% ylabel('H(\phi)');