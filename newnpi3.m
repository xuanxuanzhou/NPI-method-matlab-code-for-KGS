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
T=1;
%the T_max
lada=1;
%the coefficient of nonlinear term
%ht=waitbar(0,'Please wait...');
%the progress bar
%% fe
for fe=6:6
    %control the value of epsilon
    epsilon=1/2^(2*(fe-1));
    %epsilon=1/10^(fe-1);
    for  fN=1:1
        fN
        %control the value of M or N
        %M=2^3*2^(fN-1);
        M=6*16;
        %the spatial discretization
        
        %the temporal discretization
        %N=10*2^(fN-1);
        N=10000;
        tau=T/N; 
        tt=tau;
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
                v0(:,1)=1/2*cos(x).^2./(2-cos(x));
                v1(:,1)=1/epsilon^2*1/2*sin(x).*cos(x)./(2-cos(x));
                u0(:,1)=1+1i*sin(x)./(2-cos(x));
            otherwise
                disp('wrong number!');
        end
%         tt=h/(2*pi);
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
        %why the B1 is so important
        B1=1i*sin(tt*(-mul.^2))/tt;%1i*(-mul.^2);
        A=1i*lada./(2*epsilon^2*betal);
        vp=1/2*(v0+1i./betal.*v1);
        vm=1/2*(v0-1i./betal.*v1);
        %upload the symbol integration
        load('symboltrue.mat');
        %first-order
        for k=1:2
            %integral coefficient
            p1s(k)=double(p1{k}(tau,epsilon));
            q1s(k)=double(q1{k}(tau,epsilon));
        end
        %second-order
        for k=1:10
            p2s(k)=double(p2{k}(tau,epsilon));
        end
        for k=1:12
            q2s(k)=double(q2{k}(tau,epsilon));
        end
        %third-order
        for k=1:18
            p3s(k)=double(p3{k}(tau,epsilon));
        end
        for k=1:5
            for kk=1:4
                p33s(k,kk)=double(p33{k,kk}(tau,epsilon));
            end
        end
        for k=1:4
            p37s(k)=double(p37{k}(tau,epsilon));
        end
        for k=1:6
            for kk=1:2
                p38s(k,kk)=double(p38{k,kk}(tau,epsilon));
            end
        end
        for k=1:5
            for kk=1:8
                q33s(k,kk)=double(q33{k,kk}(tau,epsilon));
            end
        end
        for k=1:22
            q3s(k)=double(q3{k}(tau,epsilon));
        end
        for k=1:8
            q39s(k)=double(q39{k}(tau,epsilon));
        end
        %computing
        %% computing
        us=zeros(M,1);vps=us;vms=us;
        for tk=1:N                                                    %notice the tk can not be used again below
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
            %             G_21p=i*fft(F_2p.*us); G_21m=i*fft(F_2m.*us);
            %             G_22p=-i*fft(ifft(B.*fft(vps)).*us); G_22m=i*fft(ifft(B.*fft(vms)).*us);
            %             G_23p=i*fft(F_1p.*vps); G_23m=i*fft(F_1p.*vms);
            %             G_24p=i*fft(F_1m.*vps); G_24m=i*fft(F_1m.*vms);
            %             G_25p=i*fft(vps.*ifft(B1.*fft(us))); G_25m=i*fft(vms.*ifft(B1.*fft(us)));
            %             G_26p=i*B1.*fft(vps.*us); G_26m=i*B1.*fft(vms.*us);
            %
            %             F_21p=A.*fft(conj(us).*F_1p);F_21m=A.*fft(conj(us).*F_1m);
            %             F_22p=A.*fft(us.*conj(F_1p));F_22m=A.*fft(us.*conj(F_1m));
            %             F_23=A.*fft(us.*conj(ifft(B1.*fft(us))) );
            %             F_24=A.*fft(conj(us).*ifft(B1.*fft(us)) );
            %             F_25=A.*B.*fft(us.*conj(us) );
            %
            %             delta_n2psi=q1s(1)*G_21p+q1s(2)*G_21m...
            %                 +q1s(3)*G_23m+ q1s(4)*G_23p+ q1s(5)*G_24m+ q1s(6)*G_24p...
            %                              +ms(1)*(G_22p+G_25p)+ms(2)*(G_22m+G_25m)...
            %                              +ms(3)*G_26p+ms(4)*G_26m ;
            %    delta_n2p=q2s(1)*F_21p+ q2s(3)*F_21m+q3s(1)*F_22p+q3s(3)*F_22m+q4s(1)*(F_23+F_24)-q5s(1)*F_25;
            %             delta_n2m=-q2s(2)*F_21p- q2s(4)*F_21m-q3s(2)*F_22p-q3s(4)*F_22m-q4s(2)*(F_23+F_24)-q5s(2)*F_25;
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
            
            delta_n2psi=p2s(1)*F_20p+p2s(2)*F_20m...
                +p2s(3)*F_21p+p2s(4)*F_21m...
                +p2s(5)*F_22p+p2s(6)*F_22m...
                +p2s(7)*F_23p+p2s(8)*F_23m...
                +p2s(9)*F_24p+p2s(10)*F_24m;
            
            delta_n2p=q2s(1)*G_20+ q2s(3)*G_21+q2s(5)*G_22+q2s(7)*G_23+q2s(9)*G_24-q2s(11)*G_25;
            delta_n2m=-q2s(2)*G_20- q2s(4)*G_21-q2s(6)*G_22-q2s(8)*G_23-q2s(10)*G_24-q2s(12)*G_25;
            
            %third-order(p2的1，2顺序已经改好了)
            F_2{1,1}=ifft(F_20m);F_2{1,2}=ifft(F_20p);
            F_2{2,1}=ifft(F_21m);F_2{2,2}=ifft(F_21p);
            F_2{3,1}=ifft(F_22m);F_2{3,2}=ifft(F_22p);
            F_2{4,1}=ifft(F_23m);F_2{4,2}=ifft(F_23p);
            F_2{5,1}=ifft(F_24m);F_2{5,2}=ifft(F_24p);
            
            G_2{1}=ifft(G_20); G_2{2}=ifft(G_21); G_2{3}=ifft(G_22);
            G_2{4}=ifft(G_23); G_2{5}=ifft(G_24); G_2{6}=ifft(G_25);
            %空间转换
            F_30p=1i*fft(ifft(B1.*u0).*ifft(B.*vm))+ 1i/2*fft(vms.*ifft(B1.^2.*u0))+ 1i/2*fft(us.*ifft(B.^2.*vm));
            F_30m=-1i*fft(ifft(B1.*u0).*ifft(B.*vp))+ 1i/2*fft(vps.*ifft(B1.^2.*u0))+ 1i/2*fft(us.*ifft(B.^2.*vp));
            F_31p=1i/2*B1.^2.*fft(us.*vms);F_31m=1i/2*B1.^2.*fft(us.*vps);
            F_32p=1i*B1.*fft(ifft(B1.*u0).*vms+ us.*ifft(B.*vm));
            F_32m=1i*B1.*fft(ifft(B1.*u0).*vps- us.*ifft(B.*vp));
            F_35p=1i*B1.*fft(G_1m.*us); F_35m=1i*B1.*fft(G_1p.*us);
            F_36p=1i*fft(G_1m.*ifft(B1.*u0)); F_36m=1i*fft(G_1p.*ifft(B1.*u0));
            
            F_330p=1i*B1.*fft(F_1p.*vms);F_330m=1i*B1.*fft(F_1p.*vps);
            F_331p=1i*B1.*fft(F_1m.*vms);F_331m=1i*B1.*fft(F_1m.*vps);
            F_332p=1i*fft(F_1p.*ifft(B.*vm));F_332m=-1i*fft(F_1p.*ifft(B.*vp));
            F_333p=1i*fft(F_1m.*ifft(B.*vm));F_333m=-1i*fft(F_1m.*ifft(B.*vp));
            for k=1:5%将F_2写成F_2(k,1->m,2->+)%再看
                F_3{k,1}=1i*fft(F_2{k,2}.*vps);
                F_3{k,2}=1i*fft(F_2{k,2}.*vms);
                F_3{k,3}=1i*fft(F_2{k,1}.*vps);
                F_3{k,4}=1i*fft(F_2{k,1}.*vms);
            end
            F_37{1}=1i*fft(F_1p.*G_1m);
            F_37{2}=1i*fft(F_1p.*G_1p);
            F_37{3}=1i*fft(F_1m.*G_1m);
            F_37{4}=1i*fft(F_1m.*G_1p);
            for k=1:5%将G_2改写
                F_38{k,1}=1i*fft(us.*(G_2{k}));F_38{k,2}=1i*fft(us.*(-G_2{k}));
            end
            F_38{6,1}=1i*fft(us.*(-G_2{6}));F_38{6,2}=1i*fft(us.*(-G_2{6}));
            delta_n3psi=p3s(1)*F_30m+p3s(2)*F_30p+p3s(3)*F_31m+p3s(4)*F_31p+...
                p3s(5)*F_32m+p3s(6)*F_32p+p3s(15)*F_35m+p3s(16)*F_35p+...
                p3s(17)*F_36m+p3s(18)*F_36p+...
                p3s(7)*F_330m+p3s(8)*F_330p+p3s(9)*F_331m+p3s(10)*F_331p+...
                p3s(11)*F_332m+p3s(12)*F_332p+p3s(13)*F_333m+p3s(14)*F_333p+...
                p37s(1)*F_37{1}+p37s(2)*F_37{2}+p37s(3)*F_37{3}+p37s(4)*F_37{4};
            for k=1:6
                delta_n3psi=delta_n3psi+p38s(k,1)*F_38{k,1}+p38s(k,2)*F_38{k,2};
            end
            for k=1:5
                delta_n3psi=delta_n3psi+p33s(k,1)*F_3{k,1}+p33s(k,2)*F_3{k,2}+...
                    p33s(k,3)*F_3{k,3}+p33s(k,4)*F_3{k,4};
            end
            %%%%%%%%delta_n3pm%%%%%%%
            G_30p=A.*fft(1/2*us.*conj(ifft(B1.^2.*u0))+ 1/2*conj(us).*ifft(B1.^2.*u0)+ ifft(B1.*u0).*conj( ifft(B1.*u0)) );
            G_30m=-A.*fft(1/2*us.*conj(ifft(B1.^2.*u0))+ 1/2*conj(us).*ifft(B1.^2.*u0)+ ifft(B1.*u0).*conj( ifft(B1.*u0)) );
            G_31p=1/2*A.*B.^2.*fft(us.*conj(us) ); G_31m=-1/2*A.*B.^2.*fft(us.*conj(us) );
            G_32p=-A.*B.*fft( us.*conj(ifft(B1.*u0)) +conj(us).*ifft(B1.*u0) );
            G_32m=-A.*B.*fft( us.*conj(ifft(B1.*u0)) +conj(us).*ifft(B1.*u0) );
            
            G_330p=-A.*B.*fft(conj(us).*F_1p);
            G_330m=-A.*B.*fft(conj(us).*F_1p);
            G_331p=-A.*B.*fft(conj(us).*F_1m);
            G_331m=-A.*B.*fft(conj(us).*F_1m);
            G_340p=-A.*B.*fft(us.*conj(F_1p));
            G_340m=-A.*B.*fft(us.*conj(F_1p));
            G_341p=-A.*B.*fft(us.*conj(F_1m));
            G_341m=-A.*B.*fft(us.*conj(F_1m));
            
            G_370p=A.*fft(conj(ifft(B1.*u0)).*F_1p);
            G_370m=-A.*fft(conj(ifft(B1.*u0)).*F_1p);
            G_371p=A.*fft(conj(ifft(B1.*u0)).*F_1m);
            G_371m=-A.*fft(conj(ifft(B1.*u0)).*F_1m);
            G_380p=A.*fft(ifft(B1.*u0).*conj(F_1p));
            G_380m=-A.*fft(ifft(B1.*u0).*conj(F_1p));
            G_381p=A.*fft(ifft(B1.*u0).*conj(F_1m));
            G_381m=-A.*fft(ifft(B1.*u0).*conj(F_1m));
            
            for k=1:5
                G35{k,1}=A.*fft( F_2{k,2}.*conj(us)  );
                G35{k,2}=-A.*fft( F_2{k,2}.*conj(us)  );
                G35{k,3}=A.*fft( F_2{k,1}.*conj(us)  );
                G35{k,4}=-A.*fft( F_2{k,1}.*conj(us)  );
                
                G36{k,1}=A.*fft( conj(F_2{k,2}).*us  );
                G36{k,2}=-A.*fft( conj(F_2{k,2}).*us  );
                G36{k,3}=A.*fft( conj(F_2{k,1}).*us  );
                G36{k,4}=-A.*fft( conj(F_2{k,1}).*us  );
            end
            %G39可简化
            G39{1}=A.*fft(conj(F_1p).*F_1p);
            G39{2}=-A.*fft(conj(F_1p).*F_1p);
            G39{3}=A.*fft(conj(F_1m).*F_1p);
            G39{4}=-A.*fft(conj(F_1m).*F_1p);
            G39{5}=A.*fft(conj(F_1p).*F_1m);
            G39{6}=-A.*fft(conj(F_1p).*F_1m);
            G39{7}=A.*fft(conj(F_1m).*F_1m);
            G39{8}=-A.*fft(conj(F_1m).*F_1m);
            
            delta_n3p=q3s(1)*G_30p+q3s(3)*G_31p+...
                q3s(5)*G_32p+...
                q3s(7)*G_330p+ q3s(9)*G_331p+...
                q3s(11)*G_340p+ q3s(13)*G_341p+...
                q3s(15)*G_370p+q3s(17)*G_371p+...
                q3s(19)*G_380p+q3s(21)*G_381p+...
                q39s(1)*G39{1}+q39s(3)*G39{3}+...
                q39s(5)*G39{5}+q39s(7)*G39{7};
            for k=1:5
                delta_n3p=delta_n3p+q33s(k,1)*G35{k,1}+ q33s(k,3)*G35{k,3}+...
                    q33s(k,5)*G36{k,1}+ q33s(k,7)*G36{k,3};
            end

            delta_n3m=q3s(2)*G_30m+q3s(4)*G_31m+...
                q3s(6)*G_32m+...
                q3s(8)*G_330m+  q3s(10)*G_331m+...
                q3s(12)*G_340m+ q3s(14)*G_341m+...
                q3s(16)*G_370m+q3s(18)*G_371m+...
                q3s(20)*G_380m+q3s(22)*G_381m+...
                q39s(2)*G39{2}+q39s(4)*G39{4}+...
                q39s(6)*G39{6}+q39s(8)*G39{8};
            for k=1:5
                delta_n3m=delta_n3m+q33s(k,2)*G35{k,2}+ q33s(k,4)*G35{k,4}+...
                    q33s(k,6)*G36{k,2}+ q33s(k,8)*G36{k,4};
            end
%                           max(abs(delta_n3psi))/tau^3
%                          max(abs(delta_n3p))/tau^3
%                         max(abs(delta_n3m))/tau^3
            u0=exp(-1i*tau*mul.^2).*u0+ delta_n1psi+delta_n2psi+delta_n3psi;
            vp=exp(-1i*tau*betal).*vp + delta_n1p+delta_n2p+delta_n3p;
            vm=exp(1i*tau*betal).*vm + delta_n1m+delta_n2m+delta_n3m;      
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



