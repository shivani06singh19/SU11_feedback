clear all
ntime = 15;
num = 20;

r = 0.1;
ph1 = pi/4;
th = 0;
time = zeros(ntime,1);
ph = zeros(num,1);
Qf = zeros(ntime,num);
QH = zeros(ntime,num);
In = zeros(ntime,num);
Entropy = zeros(ntime,num);

%for partial transpose
mode_a = eye(2);
mode_b = [1 0;0 -1];
E = blkdiag(mode_b,mode_a);

%partially transposed symplectic matrix
om = [0 1;-1 0];
Om= blkdiag(om,om);
Transposed_Om = E*Om*E;

for l = 0:num
    ph1 = l*pi/40;    

    sig = eye(4)/2;
    dsig = zeros(4);

    U = [1 0 0 0;0 1 0 0;0 0 exp(i*ph1) 0;0 0 0 exp(-i*ph1)];
    S = [cosh(r) 0 0 -exp(i*th)*sinh(r); 0 cosh(r) -exp(-i*th)*sinh(r) 0; 0 -exp(-i*th)*sinh(r) cosh(r) 0; -exp(i*th)*sinh(r) 0 0 cosh(r)];
    S1 = [cosh(-r) 0 0 -exp(i*th)*sinh(-r); 0 cosh(-r) -exp(-i*th)*sinh(-r) 0; 0 -exp(-i*th)*sinh(-r) cosh(-r) 0; -exp(i*th)*sinh(-r) 0 0 cosh(-r)];
    dU = [0 0 0 0;0 0 0 0;0 0 i*exp(i*ph1) 0;0 0 0 -i*exp(-i*ph1)];

    for k = 0:ntime  
        proje = blkdiag(zeros(2),eye(2));

        if k < 4
            dsig = (S*dU*S*sig*S'*U'*S' + S*U*S*dsig*S'*U'*S' + S*U*S*sig*S'*dU'*S');
            sig = (S*U*S*sig*S'*U'*S');
        elseif k<9
            dsig = (S1*dU*S1*sig*S1'*U'*S1' + S1*U*S1*dsig*S1'*U'*S1' + S1*U*S1*sig*S1'*dU'*S1');
            sig = (S1*U*S1*sig*S1'*U'*S1'); 
        elseif k<14
            dsig = (S*dU*S*sig*S'*U'*S' + S*U*S*dsig*S'*U'*S' + S*U*S*sig*S'*dU'*S');
            sig = (S*U*S*sig*S'*U'*S');
        elseif k<19
            dsig = (S1*dU*S1*sig*S1'*U'*S1' + S1*U*S1*dsig*S1'*U'*S1' + S1*U*S1*sig*S1'*dU'*S1');
            sig = (S1*U*S1*sig*S1'*U'*S1'); 
        elseif k<24
            dsig = (S*dU*S*sig*S'*U'*S' + S*U*S*dsig*S'*U'*S' + S*U*S*sig*S'*dU'*S');
            sig = (S*U*S*sig*S'*U'*S');
        else
            dsig = (S1*dU*S1*sig*S1'*U'*S1' + S1*U*S1*dsig*S1'*U'*S1' + S1*U*S1*sig*S1'*dU'*S1');
            sig = (S1*U*S1*sig*S1'*U'*S1'); 
        end

    % J = Om*sig;
    % 
    % Y = proje*J*proje;
    % 
    % [Ev,Eg] = eig(Y);
    % S2=0;
    % for j = 1:4
    %     if abs(Eg(j,j)) > 1
    %         la = abs(Eg(j,j));
    %         S2 = S2+((la+1)*log2((la+1)/2) - (la-1)*log2((la-1)/2))/2;
    %     end
    % end
     
        H = trace(dsig*inv(sig)*dsig*inv(sig))/4; 
     
    %partially transposed covaraince matrix
        Transposed_sig = E*sig*E;

    %symplectic transformation of covaraince matrix
        J = (Transposed_sig*Transposed_Om*Transposed_sig*Transposed_Om);

        [Ev,Eg] = eig(-J);

        S11=0;
        S22=0;
        for a = 1:4
            S11 = log2(2*sqrt(Eg(a,a)));
            if S11< 0
                S22=S22+S11;
            end
        end
         
        time(k+1,1) = k;
        Qf(k+1,l+1) = H; %trace(H);
        In(k+1,l+1) = sig(1,1) + sig(3,3) - 1;
        QH(k+1,l+1) = 1/sqrt(H);
        Entropy(k+1,l+1) = -S22;
        ph(l+1,1) = l;

    end
 end

hold on;
Z = abs(Entropy);
surf(ph,time,abs(Qf))
% ribbon(Entropy)
% plot(time,Entropy,'g')
xlabel('\phi');
ylabel('No. of loops');
zlabel('H')
