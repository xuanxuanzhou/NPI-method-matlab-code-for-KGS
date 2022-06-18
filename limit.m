%the fourier method for the limit equation
%iu_t+u_xx=0
%2iv_t-v_xx=0
%u0,v0=1/2(phi_0-iphi_1)
clc;
clear;
%% preperation
T=10;
lada=1;
%% fe
for fe=1:6
    %control the value of epsilon
    epsilon=1/10^(1*(fe));
    for  fN=1:1
        fN
        %control the value of M or N
        %M=2^3*2^(fN-1);
        M=6*32;%64*16;%
        N=1000000;
        %N=10*2^(fN-1);
        tau=T/N; tt=tau;
        %% initial value
        u0=zeros(M,1);
        v0=u0;v1=u0;
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
                v1(:,1)=1/(sqrt(2))*exp(-x.^2);
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
                u0(:,1)=1+1i*sin(x)./(2-cos(x));
                v0(:,1)=1/2*cos(x).^2./(2-cos(x));
                v1(:,1)=1/2*sin(x).*cos(x)./(2-cos(x));
            otherwise
                disp('wrong number!');
        end
        
        mul=zeros(M,1); 
        for l=1:M
            %notice the solver of 0->0(the first method for l: 0~M/2-1<=>1~M/2, -M/2~-1<=>M/2+1~M)
            if l>M/2
                mul(l)=2*pi*(l-M-1)/(b-a);
            else
                mul(l)=2*pi*(l-1)/(b-a);
            end
        end
        u0=fft(u0);
        v0=fft(v0); v1=fft(v1);
        
        %initial value
        v2=1/2*(v0-1i*v1);
        %computing
        %% computing
        %u->n  vp->u
        us=zeros(M,1);vs=us;
        ulimit{1,fN,fe}=ifft(u0);
        vlimit{1,fN,fe}=1*(exp(1i*(0)*tau/epsilon^2)*ifft(v2)+exp(-1i*(0)*tau/epsilon^2)*conj(ifft(v2)));
        for tk=1:N
            us=u0.*exp(-1i*mul.^2*tau);
            vs=v2.*exp(1i/2*mul.^2*tau);
            u0=us;
            v2=vs;
             if mod(tk,1000)==0
                ulimit{tk/1000+1,fN,fe}=ifft(u0);
                vlimit{tk/1000+1,fN,fe}=1*(exp(1i*(tk)*tau/epsilon^2)*ifft(vs)+exp(-1i*(tk)*tau/epsilon^2)*conj(ifft(vs)));
            end
        end
        %waitbar(fN*fe/21,ht);                                       %show the progress bar
        u{fN,fe}=ifft(u0);
        v{fN,fe}=1/2*(exp(1i*(tk)*tau/epsilon^2)*ifft(vs)+exp(-1i*(tk)*tau/epsilon^2)*conj(ifft(vs)));
        if fN>1
            err1(fN,fe)=sqrt(h)*norm(u{fN,fe}-u{fN-1,fe},'fro');
            err2(fN,fe)=sqrt(h)*norm(v{fN,fe}-v{fN-1,fe},'fro');
        end
    end%fN
    %t(fN,fe)=toc;                                                         %save the cpu time
end%fe
%close(ht);                                                      %close the progress bar
%end%function



