%% handbook
%the program is for the complex-value |u|^2u
%Klein-Gordon equation in the non-relativistic limit regime,
%first-order Nested Picard method
%with Fourier pseudospectral method,
%function: function [us,t]=newnpi1(flag)
%input:  the explanation of flag can be found below
%output: us is the numerical solution
%(the index of column is fe, the index of row is fN, t is the same as us),
%t is the cpu time
%%  function
%function [err,vs,t]=newnpi22(flag)
clc;
clear;

tic
%% preperation
T=10;
%the T_max
lada=1;
%the coefficient of nonlinear term
%ht=waitbar(0,'Please wait...');
%the progress bar
%% fe
for fe=1:8
    %control the value of epsilon
    epsilon=1/5^(1*(fe));
    %epsilon=1/2^(2*(fe-1));
    for  fN=1:1
        fN
        %control the value of M or N
        %M=128*2^(fN-1);
        M=6*32;
        %the spatial discretization
        
        %the temporal discretization
        %N=10*2^(fN-1);
        %N=1000000;
        N=1000000;
        tau=T/N; tt=tau;
        %% initial value
        u0=zeros(M,1);
        u1=u0;
        flag1=4;
        switch flag1
            %the value of flag1 can be 1,2,3
            % flag1=1,  [a,b]=[-32,32]   M is changed with a,b
            case 1  %来自文章：kgs regime
                a=-32;b=32;
                h=(b-a)/M;
                x=linspace(a,b-h,M);
                u0(:,1)=(1+1i)/2*sech(x.^2);
                v0(:,1)=1/2*exp(-x.^2);
                v1(:,1)=1/(sqrt(2)*epsilon^2)*exp(-x.^2);
                % flag1=2,  [a,b]=[-512,512]   M is changed with a,b
            case 2   %来自文章：kgs regime
                a=-512;b=512;
                h=(b-a)/M;
                x=linspace(a,b-h,M);
                u0(:,1)=(x.^m.*abs(x))/2.*sech(x.^2);
                v0(:,1)=(x.^m.*abs(x))/2.*exp(-x.^2);
                v1(:,1)=1/(sqrt(2)*epsilon^2)*exp(-x.^2);
                % flag1=3,  [a,b]=[-32,32]   M is changed with a,b
            case 3  %来自文章：baosu
                a=-30;b=30;
                h=(b-a)/M;
                x=linspace(a,b-h,M);
                u0(:,1)=exp(-x.^2/2+1i*x/2);
                v0(:,1)=0;
                v1(:,1)=0;
            case 4  %来自文章：baosu
                a=0;b=2*pi;
                h=(b-a)/M;
                x=linspace(a,b-h,M);
                u0(:,1)=1+1i*sin(x)./(2-cos(x));
                v0(:,1)=1/2*cos(x).^2./(2-cos(x));
                v1(:,1)=1/epsilon^2*1/2*sin(x).*cos(x)./(2-cos(x));
            otherwise
                disp('wrong number!');
        end
        
        mul=zeros(M,1); betal=mul;
        betag=mul;B=mul;A=mul;
        for l=1:M
            %notice the solver of 0->0(the first method for l: 0~M/2-1<=>1~M/2, -M/2~-1<=>M/2+1~M)
            if l>M/2
                mul(l)=2*pi*(l-M-1)/(b-a);
            else
                mul(l)=2*pi*(l-1)/(b-a);
            end
            betal(l,1)=sqrt(1+epsilon^2*mul(l)^2)/(epsilon^2);
        end
        u0=fft(u0);
        v0=fft(v0); v1=fft(v1);
        betag=betal-1/epsilon^2;
        B=1i*sin(tt*(betag))/tt;
        B1=1i*sin(tt*(-mul.^2))/tt;
        A=1i*lada./(2*epsilon^2*betal);
        vp=1/2*(v0+1i./betal.*v1);
        vm=1/2*(v0-1i./betal.*v1);
        %upload the symbol integration
        load('symboltrue.mat');
        for k=1:2
            %integral coefficient
            p1s(k)=double(p1{k}(tau,epsilon));
            q1s(k)=double(q1{k}(tau,epsilon));
        end
        for k=1:10
            p2s(k)=double(p2{k}(tau,epsilon));
        end
        for k=1:12
            q2s(k)=double(q2{k}(tau,epsilon));
        end
        %computing
        %% computing
        us=zeros(M,1);vps=us;vms=us;
         unum{1,fN,fe}=ifft(u0);
         vnum{1,fN,fe}=ifft(vp+vm); 
        for tk=1:N                                                  %notice the tk can not be used again below
            %delta^{n,1}_{+} delta_n1p
            us=ifft(u0);
            vps=ifft(vp);
            vms=ifft(vm);
            F_1p=fft(1i*us.*vps);       %F_{1,+}
            F_1m=fft(1i*us.*vms);       %F_{1,-}
            G_1p=+A.*fft(us.*conj(us)); %G_{1,+}
            G_1m=-A.*fft(us.*conj(us)); %G_{1,-}
            
            delta_n1psi=p1s(1)*F_1p+p1s(2)*F_1m;
            delta_n1p=q1s(1)*G_1p;
            delta_n1m=q1s(2)*G_1m;
            
            F_1p=ifft(F_1p);F_1m=ifft(F_1m);G_1p=ifft(G_1p);G_1m=ifft(G_1m);
            %second-order
            F_24p=1i*fft(G_1p.*us); F_24m=1i*fft(G_1m.*us);
            F_22p=-1i*fft(ifft(B.*fft(vps)).*us)+ 1i*fft(vps.*ifft(B1.*fft(us))) ;
            F_22m=1i*fft(ifft(B.*fft(vms)).*us)+ 1i*fft(vms.*ifft(B1.*fft(us)));
            F_20p=1i*fft(F_1p.*vps); F_20m=1i*fft(F_1p.*vms);
            F_21p=1i*fft(F_1m.*vps); F_21m=1i*fft(F_1m.*vms);
            F_23p=1i*B1.*fft(vps.*us); F_23m=1i*B1.*fft(vms.*us);
            
            G_20=A.*fft(conj(us).*F_1p);G_21=A.*fft(conj(us).*F_1m);
            G_22=A.*fft(us.*conj(F_1p));G_23=A.*fft(us.*conj(F_1m));
            G_24=A.*fft(us.*conj(ifft(B1.*fft(us))) )+ A.*fft(conj(us).*ifft(B1.*fft(us)) );
            G_25=A.*B.*fft(us.*conj(us) );
            
            delta_n2psi=p2s(3)*F_20p+p2s(4)*F_20m...
                +p2s(5)*F_21p+p2s(6)*F_21m...
                +p2s(7)*F_22p+p2s(8)*F_22m...
                +p2s(9)*F_23p+p2s(10)*F_23m...
                +p2s(1)*F_24p+p2s(2)*F_24m;
%             delta_n2psi=p2s(1)*F_20p+p2s(2)*F_20m...
%                 +p2s(3)*F_21p+p2s(4)*F_21m...
%                 +p2s(5)*F_22p+p2s(6)*F_22m...
%                 +p2s(7)*F_23p+p2s(8)*F_23m...
%                 +p2s(9)*F_24p+p2s(10)*F_24m;
            
            delta_n2p=q2s(1)*G_20+ q2s(3)*G_21+q2s(5)*G_22+q2s(7)*G_23+q2s(9)*G_24-q2s(11)*G_25;
            delta_n2m=-q2s(2)*G_20- q2s(4)*G_21-q2s(6)*G_22-q2s(8)*G_23-q2s(10)*G_24-q2s(12)*G_25;

            u0=exp(-1i*tau*mul.^2).*u0+ delta_n1psi+delta_n2psi;
            vp=exp(-1i*tau*betal).*vp + delta_n1p+delta_n2p;
            vm=exp(1i*tau*betal).*vm + delta_n1m+delta_n2m;
            if mod(tk,1000)==0
                unum{(tk)/1000+1,fN,fe}=ifft(u0);
                vnum{(tk)/1000+1,fN,fe}=ifft(vp+vm); 
            end
        end
        %waitbar(fN*fe/21,ht);                                       %show the progress bar
        u{fN,fe}=ifft(u0);
        vs{fN,fe}=ifft(vp+vm);                                              %save the numerical solution
        if fN>1
            err1(fN,fe)=sqrt(h)*norm(u{fN,fe}-u{fN-1,fe},'fro');
            err2(fN,fe)=sqrt(h)*norm(vs{fN,fe}-vs{fN-1,fe},'fro');
        end
 
        t(fN,fe)=toc;
    end%fN   
end%fe
%close(ht);                                                      %close the progress bar
%end%function



