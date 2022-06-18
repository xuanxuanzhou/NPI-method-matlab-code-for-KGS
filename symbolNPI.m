%% 符号计算
%tips 1: 符号计算中可以嵌套共轭运算
%tips 2: 符号计算可以单独成文件先算，保存，载入直接赋值用symbolNPI_complex.m
function symbolNPI(savefile,x)
%x=1->first-order NPI symbol integration
%x=2->first-order NPI symbol integration
%x=3->first-order NPI symbol integration
syms s w epsilo real;
%先判断x是否为整数
flagg=0;%用于指示是否已保存符号积分
if x==1 | x==2 | x==3
    %% first-order used funciton
    p1{1}=symfun(int(exp(-1i*w/epsilo^2),w,0,s),[w,epsilo]);    %p_{1,0},p_{1,-}
    p1{2}=symfun(int(exp(1i*w/epsilo^2),w,0,s),[w,epsilo]);     %p_{1,1},p_{1,+}
    q1{1}=symfun(int(exp(-1i*(s-w)/epsilo^2),w,0,s),[w,epsilo]);%q_{1,-}
    q1{2}=symfun(int(exp(1i*(s-w)/epsilo^2),w,0,s),[w,epsilo]); %q_{1,+}
    for k=1:2
        p1{k}=subs(p1{k},s,w);
        q1{k}=subs(q1{k},s,w);
    end
    if x==2 | x==3
        %% second-order used funciton
        %         p2{3}=symfun(int(exp(1i*w/epsilo^2)*p1{1},w,0,s),[w,epsilo]); %p_{2,0,+}
        %         p2{4}=symfun(int(exp(-1i*w/epsilo^2)*p1{1},w,0,s),[w,epsilo]);%p_{2,0,-}
        %         p2{5}=symfun(int(exp(1i*w/epsilo^2)*p1{2},w,0,s),[w,epsilo]); %p_{2,1,+}
        %         p2{6}=symfun(int(exp(-1i*w/epsilo^2)*p1{2},w,0,s),[w,epsilo]);%p_{2,1,-}
        
        p2{1}=symfun(int(exp(-1i*w/epsilo^2)*p1{1},w,0,s),[w,epsilo]); %p_{2,0,-}
        p2{2}=symfun(int(exp(1i*w/epsilo^2)*p1{1},w,0,s),[w,epsilo]);%p_{2,0,+}
        p2{3}=symfun(int(exp(-1i*w/epsilo^2)*p1{2},w,0,s),[w,epsilo]); %p_{2,1,-}
        p2{4}=symfun(int(exp(1i*w/epsilo^2)*p1{2},w,0,s),[w,epsilo]);%p_{2,1,+}
        
        p2{5}=symfun(int(exp(-1i*w/epsilo^2)*w,w,0,s),[w,epsilo]);    %p_{2,2,-}
        p2{6}=symfun(int(exp(1i*w/epsilo^2)*w,w,0,s),[w,epsilo]);     %p_{2,2,+}
        p2{7}=symfun(int(exp(-1i*w/epsilo^2)*(s-w),w,0,s),[w,epsilo]);%p_{2,3,-}
        p2{8}=symfun(int(exp(1i*w/epsilo^2)*(s-w),w,0,s),[w,epsilo]); %p_{2,3,+}
        p2{9}=symfun(int(q1{1},w,0,s),[w,epsilo]);%1->m,2-p,           p_{2,4,-}
        p2{10}=symfun(int(q1{2},w,0,s),[w,epsilo]);                    %p_{2,4,+}
        for k=1:10
            p2{k}=subs(p2{k},s,w);
        end
        q2{1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*p1{1},w,0,s),[w,epsilo]);     %q_{2,+,0}
        q2{2}=symfun(int(exp(1i*(s-w)/epsilo^2)*p1{1},w,0,s),[w,epsilo]);      %q_{2,-,0}
        q2{3}=symfun(int(exp(-1i*(s-w)/epsilo^2)*p1{2},w,0,s),[w,epsilo]);     %q_{2,+,1}
        q2{4}=symfun(int(exp(1i*(s-w)/epsilo^2)*p1{2},w,0,s),[w,epsilo]);      %q_{2,-,1}
        
        q2{5}=symfun(int(exp(-1i*(s-w)/epsilo^2)*conj(p1{1}),w,0,s),[w,epsilo]); %q_{\bar{2},+,0}
        q2{6}=symfun(int(exp(1i*(s-w)/epsilo^2)*conj(p1{1}),w,0,s),[w,epsilo]);  %q_{\bar{2},-,0}
        q2{7}=symfun(int(exp(-1i*(s-w)/epsilo^2)*conj(p1{2}),w,0,s),[w,epsilo]); %q_{\bar{2},+,1}
        q2{8}=symfun(int(exp(1i*(s-w)/epsilo^2)*conj(p1{2}),w,0,s),[w,epsilo]);  %q_{\bar{2},-,1}
        
        q2{9}=symfun(int(exp(-1i*(s-w)/epsilo^2)*w,w,0,s),[w,epsilo]);    %q_{2,+,4}
        q2{10}=symfun(int(exp(1i*(s-w)/epsilo^2)*w,w,0,s),[w,epsilo]);     %q_{2,-,4}
        q2{11}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(s-w),w,0,s),[w,epsilo]);%q_{2,+,5}
        q2{12}=symfun(int(exp(1i*(s-w)/epsilo^2)*(s-w),w,0,s),[w,epsilo]); %q_{2,-,5}
        for k=1:12
            q2{k}=subs(q2{k},s,w);
        end
        if x==3
            %%  third-order used function
            p3{1}=symfun(int(exp(-1i*w/epsilo^2)*w^2,w,0,s),[w,epsilo]);     %p_{3,0,-}
            p3{2}=symfun(int(exp(1i*w/epsilo^2)*w^2,w,0,s),[w,epsilo]);     %p_{3,0,+}
            p3{3}=symfun(int(exp(-1i*w/epsilo^2)*(s-w)^2,w,0,s),[w,epsilo]);     %p_{3,1,-}
            p3{4}=symfun(int(exp(1i*w/epsilo^2)*(s-w)^2,w,0,s),[w,epsilo]);     %p_{3,1,+}
            p3{5}=symfun(int(exp(-1i*w/epsilo^2)*(s-w)*w,w,0,s),[w,epsilo]);     %p_{3,2,-}
            p3{6}=symfun(int(exp(1i*w/epsilo^2)*(s-w)*w,w,0,s),[w,epsilo]);     %p_{3,2,+}
            
            p3{7}=symfun(int(exp(-1i*w/epsilo^2)*(s-w)*p1{1},w,0,s),[w,epsilo]);     %p_{3,3,0,-}
            p3{8}=symfun(int(exp(1i*w/epsilo^2)*(s-w)*p1{1},w,0,s),[w,epsilo]);     %p_{3,3,0,+}
            p3{9}=symfun(int(exp(-1i*w/epsilo^2)*(s-w)*p1{2},w,0,s),[w,epsilo]);     %p_{3,3,1,-}
            p3{10}=symfun(int(exp(1i*w/epsilo^2)*(s-w)*p1{2},w,0,s),[w,epsilo]);     %p_{3,3,1,+}
            p3{11}=symfun(int(exp(-1i*w/epsilo^2)*w*p1{1},w,0,s),[w,epsilo]);     %p_{3,3,2,-}
            p3{12}=symfun(int(exp(1i*w/epsilo^2)*w*p1{1},w,0,s),[w,epsilo]);     %p_{3,3,2,+}
            p3{13}=symfun(int(exp(-1i*w/epsilo^2)*w*p1{2},w,0,s),[w,epsilo]);     %p_{3,3,3,-}
            p3{14}=symfun(int(exp(1i*w/epsilo^2)*w*p1{2},w,0,s),[w,epsilo]);     %p_{3,3,3,+}
            
            p3{15}=symfun(int(q1{1}*(s-w),w,0,s),[w,epsilo]);     %p_{3,5,-}
            p3{16}=symfun(int(q1{2}*(s-w),w,0,s),[w,epsilo]);     %p_{3,5,+}
            p3{17}=symfun(int(q1{1}*w,w,0,s),[w,epsilo]);     %p_{3,6,-}
            p3{18}=symfun(int(q1{2}*w,w,0,s),[w,epsilo]);     %p_{3,6,+}
            for k=1:18
                p3{k}=subs(p3{k},s,w);
            end
            for k=1:2:9
                p33{(k-1)/2+1,1}=symfun(int(exp(-1i*w/epsilo^2)*p2{k},w,0,s),[w,epsilo]);     %p_{3,k,-,-}
                p33{(k-1)/2+1,2}=symfun(int(exp(1i*w/epsilo^2)*p2{k},w,0,s),[w,epsilo]);     %p_{3,k,-,+}
                p33{(k-1)/2+1,3}=symfun(int(exp(-1i*w/epsilo^2)*p2{k+1},w,0,s),[w,epsilo]);     %p_{3,k,+,-}
                p33{(k-1)/2+1,4}=symfun(int(exp(1i*w/epsilo^2)*p2{k+1},w,0,s),[w,epsilo]);     %p_{3,k,+,+}
                p33{(k-1)/2+1,1}=subs(p33{(k-1)/2+1,1},s,w);
                p33{(k-1)/2+1,2}=subs(p33{(k-1)/2+1,2},s,w);
                p33{(k-1)/2+1,3}=subs(p33{(k-1)/2+1,3},s,w);
                p33{(k-1)/2+1,4}=subs(p33{(k-1)/2+1,4},s,w);
            end
            for k=1:2:11
                p38{(k-1)/2+1,1}=symfun(int(q2{k},w,0,s),[w,epsilo]);     %p_{3,8,k,+}
                p38{(k-1)/2+1,2}=symfun(int(q2{k+1},w,0,s),[w,epsilo]);     %p_{3,8,k,-}
                p38{(k-1)/2+1,1}=subs(p38{(k-1)/2+1,1},s,w);
                p38{(k-1)/2+1,2}=subs(p38{(k-1)/2+1,2},s,w);
            end
            p37{1}=symfun(int(p1{1}*q1{2},w,0,s),[w,epsilo]);     %p_{3,7,0,+}
            p37{2}=symfun(int(p1{1}*q1{1},w,0,s),[w,epsilo]);     %p_{3,7,0,-}
            p37{3}=symfun(int(p1{2}*q1{2},w,0,s),[w,epsilo]);     %p_{3,7,1,+}
            p37{4}=symfun(int(p1{2}*q1{1},w,0,s),[w,epsilo]);     %p_{3,7,1,-}
            for k=1:4
                p37{k}=subs(p37{k},s,w);
            end
            %%%%%%%%%%%%%%%%%%%
            q3{1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*w^2,w,0,s),[w,epsilo]);     %q_{3,0,+}
            q3{2}=symfun(int(exp(1i*(s-w)/epsilo^2)*w^2,w,0,s),[w,epsilo]);     %q_{3,0,-}
            q3{3}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(s-w)^2,w,0,s),[w,epsilo]);     %q_{3,1,+}
            q3{4}=symfun(int(exp(1i*(s-w)/epsilo^2)*(s-w)^2,w,0,s),[w,epsilo]);     %q_{3,1,-}
            q3{5}=symfun(int(exp(-1i*(s-w)/epsilo^2)*w*(s-w),w,0,s),[w,epsilo]);     %q_{3,2,+}
            q3{6}=symfun(int(exp(1i*(s-w)/epsilo^2)*w*(s-w),w,0,s),[w,epsilo]);     %q_{3,2,-}
            
            q3{7}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(s-w)*p1{1},w,0,s),[w,epsilo]);     %q_{3,3,0,+}
            q3{8}=symfun(int(exp(1i*(s-w)/epsilo^2)*(s-w)*p1{1},w,0,s),[w,epsilo]);     %q_{3,3,0,-}
            q3{9}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(s-w)*p1{2},w,0,s),[w,epsilo]);     %q_{3,3,1,+}
            q3{10}=symfun(int(exp(1i*(s-w)/epsilo^2)*(s-w)*p1{2},w,0,s),[w,epsilo]);     %q_{3,3,1,-}
            
            q3{11}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(s-w)*conj(p1{1}),w,0,s),[w,epsilo]);     %q_{3,4,0,+}
            q3{12}=symfun(int(exp(1i*(s-w)/epsilo^2)*(s-w)*conj(p1{1}),w,0,s),[w,epsilo]);     %q_{3,4,0,-}
            q3{13}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(s-w)*conj(p1{2}),w,0,s),[w,epsilo]);     %q_{3,4,1,+}
            q3{14}=symfun(int(exp(1i*(s-w)/epsilo^2)*(s-w)*conj(p1{2}),w,0,s),[w,epsilo]);     %q_{3,4,1,-}
            
            q3{15}=symfun(int(exp(-1i*(s-w)/epsilo^2)*w*p1{1},w,0,s),[w,epsilo]);     %q_{3,7,0,+}
            q3{16}=symfun(int(exp(1i*(s-w)/epsilo^2)*w*p1{1},w,0,s),[w,epsilo]);     %q_{3,7,0,-}
            q3{17}=symfun(int(exp(-1i*(s-w)/epsilo^2)*w*p1{2},w,0,s),[w,epsilo]);     %q_{3,7,1,+}
            q3{18}=symfun(int(exp(1i*(s-w)/epsilo^2)*w*p1{2},w,0,s),[w,epsilo]);     %q_{3,7,1,-}
            
            q3{19}=symfun(int(exp(-1i*(s-w)/epsilo^2)*w*conj(p1{1}),w,0,s),[w,epsilo]);     %q_{3,8,0,+}
            q3{20}=symfun(int(exp(1i*(s-w)/epsilo^2)*w*conj(p1{1}),w,0,s),[w,epsilo]);     %q_{3,8,0,-}
            q3{21}=symfun(int(exp(-1i*(s-w)/epsilo^2)*w*conj(p1{2}),w,0,s),[w,epsilo]);     %q_{3,8,1,+}
            q3{22}=symfun(int(exp(1i*(s-w)/epsilo^2)*w*conj(p1{2}),w,0,s),[w,epsilo]);     %q_{3,8,1,-}
            for k=1:22
                q3{k}=subs(q3{k},s,w);
            end
            for k=1:2:9
                q33{(k-1)/2+1,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*p2{k},w,0,s),[w,epsilo]);     %p_{3,5,k,-,+}
                q33{(k-1)/2+1,2}=symfun(int(exp(1i*(s-w)/epsilo^2)*p2{k},w,0,s),[w,epsilo]);     %p_{3,5,k,-,-}
                q33{(k-1)/2+1,3}=symfun(int(exp(-1i*(s-w)/epsilo^2)*p2{k+1},w,0,s),[w,epsilo]);     %p_{3,5,k,+,+}
                q33{(k-1)/2+1,4}=symfun(int(exp(1i*(s-w)/epsilo^2)*p2{k+1},w,0,s),[w,epsilo]);     %p_{3,5,k,+,-}
                q33{(k-1)/2+1,5}=symfun(int(exp(-1i*(s-w)/epsilo^2)*conj(p2{k}),w,0,s),[w,epsilo]);     %p_{3,6,k,-,+}
                q33{(k-1)/2+1,6}=symfun(int(exp(1i*(s-w)/epsilo^2)*conj(p2{k}),w,0,s),[w,epsilo]);     %p_{3,6,k,-,-}
                q33{(k-1)/2+1,7}=symfun(int(exp(-1i*(s-w)/epsilo^2)*conj(p2{k+1}),w,0,s),[w,epsilo]);     %p_{3,6,k,+,+}
                q33{(k-1)/2+1,8}=symfun(int(exp(1i*(s-w)/epsilo^2)*conj(p2{k+1}),w,0,s),[w,epsilo]);     %p_{3,6,k,+,-}
                q33{(k-1)/2+1,1}=subs(q33{(k-1)/2+1,1},s,w);
                q33{(k-1)/2+1,2}=subs(q33{(k-1)/2+1,2},s,w);
                q33{(k-1)/2+1,3}=subs(q33{(k-1)/2+1,3},s,w);
                q33{(k-1)/2+1,4}=subs(q33{(k-1)/2+1,4},s,w);
                q33{(k-1)/2+1,5}=subs(q33{(k-1)/2+1,5},s,w);
                q33{(k-1)/2+1,6}=subs(q33{(k-1)/2+1,6},s,w);
                q33{(k-1)/2+1,7}=subs(q33{(k-1)/2+1,7},s,w);
                q33{(k-1)/2+1,8}=subs(q33{(k-1)/2+1,8},s,w);
            end
            q39{1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*p1{1}*conj(p1{1}),w,0,s),[w,epsilo]);     %q_{3,9,0,0,+}
            q39{2}=symfun(int(exp(1i*(s-w)/epsilo^2)*p1{1}*conj(p1{1}),w,0,s),[w,epsilo]);     %q_{3,9,0,0,-}
            q39{3}=symfun(int(exp(-1i*(s-w)/epsilo^2)*p1{1}*conj(p1{2}),w,0,s),[w,epsilo]);     %q_{3,9,0,1,+}
            q39{4}=symfun(int(exp(1i*(s-w)/epsilo^2)*p1{1}*conj(p1{2}),w,0,s),[w,epsilo]);     %q_{3,9,0,1,-}
            
            q39{5}=symfun(int(exp(-1i*(s-w)/epsilo^2)*p1{2}*conj(p1{1}),w,0,s),[w,epsilo]);     %q_{3,9,1,0,+}
            q39{6}=symfun(int(exp(1i*(s-w)/epsilo^2)*p1{2}*conj(p1{1}),w,0,s),[w,epsilo]);     %q_{3,9,1,0,-}
            q39{7}=symfun(int(exp(-1i*(s-w)/epsilo^2)*p1{2}*conj(p1{2}),w,0,s),[w,epsilo]);     %q_{3,9,1,1,+}
            q39{8}=symfun(int(exp(1i*(s-w)/epsilo^2)*p1{2}*conj(p1{2}),w,0,s),[w,epsilo]);     %q_{3,9,1,1,-}
            for k=1:8
                q39{k}=subs(q39{k},s,w);
            end
            flagg=flagg+1;
        else
            disp('x=2');
        end%if
        flagg=flagg+1;
    else
        disp('x=1');
    end
    flagg=flagg+1;
else
    disp('error');
end

if flagg==3
    save(savefile,'p1','p2','q1','q2','q3','p3','p33','p38','p37','q39','q33');
elseif flagg ==2
    save(savefile,'p1','q1','q2','p2');
else
    save(savefile,'p1','q1');
end%function





