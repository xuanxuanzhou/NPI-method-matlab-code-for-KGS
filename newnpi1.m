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
%function [err,vs,t]=newnpi1(flag)
clc;
clear;
tic
%% preperation
T=1;
syms s w epsilo real;
%% first-order used funciton
    p1{1}=symfun(int(exp(-1i*w/epsilo^2),w,0,s),[w,epsilo]);    %p_{1,0},p_{1,-}
    p1{2}=symfun(int(exp(1i*w/epsilo^2),w,0,s),[w,epsilo]);     %p_{1,1},p_{1,+}
    q1{1}=symfun(int(exp(-1i*(s-w)/epsilo^2),w,0,s),[w,epsilo]);%q_{1,-}
    q1{2}=symfun(int(exp(1i*(s-w)/epsilo^2),w,0,s),[w,epsilo]); %q_{1,+}
    for k=1:2
        p1{k}=subs(p1{k},s,w);
        q1{k}=subs(q1{k},s,w);
    end
lada=1;
%the coefficient of nonlinear term
%ht=waitbar(0,'Please wait...');
%the progress bar
%% fe
for fe=1:6
    %control the value of epsilon
    epsilon=1/2^(2*(fe-1));
    %epsilon=1/2^(fe-1);
    for  fN=1:7
        fN
        %control the value of M or N
        %M=2^3*2^(fN-1);
        M=6*32;
        %the spatial discretization
        
        %the temporal discretization
        N=10*2^(fN-1);
        %N=100000;
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
                m=1;
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
                v0(:,1)=1/2*cos(x).^2./(2-cos(x));
                v1(:,1)=1/epsilon^2*1/2*sin(x).*cos(x)./(2-cos(x));
                u0(:,1)=1+1i*sin(x)./(2-cos(x));
            otherwise
                disp('wrong number!');
        end
        
        mul=zeros(M,1); betal=mul;
        betag=mul;B=mul;A=mul;up=mul;um=mul;
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
        A=1i*lada./(2*epsilon^2*betal);
        vp=1/2*(v0+1i./betal.*v1);
        vm=1/2*(v0-1i./betal.*v1);
%          %upload the symbol integration
%         load('symboltrue.mat');
        for k=1:2
            %integral coefficient
            p1s(k)=double(p1{k}(tau,epsilon));
            q1s(k)=double(q1{k}(tau,epsilon));
        end
        %computing
        uz0(:,1)=ifft(u0);
        vsz0(:,1)=ifft(vp+vm);
        uz(1,fe)=uz0(M/2+1,1);
        vsz(1,fe)=vsz0(M/2+1,1);
        %% computing
        us=zeros(M,1);vps=us;vms=us;F_1p=us;F_1m=us;F_2p=us;F_2m=us;delta_n1psi=us;delta_n1p=us;delta_n1m=us;
        for tk=1:N                                                    %notice the tk can not be used again below
            %delta^{n,1}_{+} delta_n1p
            us=ifft(u0);
            vps=ifft(vp);
            vms=ifft(vm);
            F_1p=fft(1i*us.*vps);
            F_1m=fft(1i*us.*vms);
%             G_1p=+A.*fft(us.^2);
%             G_1m=-A.*fft(us.^2);
            G_1p=+A.*fft(us.*conj(us));
            G_1m=-A.*fft(us.*conj(us));
            
            delta_n1psi=p1s(1)*F_1p+p1s(2)*F_1m;
            delta_n1p=q1s(1)*G_1p;
            delta_n1m=q1s(2)*G_1m;
            
            us=exp(-1i*tau*mul.^2).*u0+ delta_n1psi;
            vps=exp(-1i*tau*betal).*vp+ delta_n1p;
            vms=exp(1i*tau*betal).*vm+ delta_n1m;
            u0=us;
            vp=vps;
            vm=vms;
            uz0(:,1)=ifft(u0);
            vsz0(:,1)=ifft(vp+vm);
            uz(tk+1,fe)=uz0(M/2+1,1);
            vsz(tk+1,fe)=vsz0(M/2+1,1);
        end
        %waitbar(fN*fe/21,ht);                                       %show the progress bar
        u{fN,fe}=ifft(u0);
        vs{fN,fe}=ifft(vp+vm);                                              %save the numerical solution
%         if fN>1
%             err1(fN,fe)=sqrt(h)*norm(u{fN,fe}-u{fN-1,fe},'fro');
%             err2(fN,fe)=sqrt(h)*norm(vs{fN,fe}-vs{fN-1,fe},'fro');
%         end
    end%fN
    t(fN,fe)=toc;                                                         %save the cpu time
end%fe
%close(ht);                                                      %close the progress bar
%end%function



