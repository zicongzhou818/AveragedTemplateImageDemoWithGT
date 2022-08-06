function [T1_new,T2_new]=PJCDtestPhi_fixedboundary_v1(Avg_J_T_1ab,Avg_Curl_Phi,N,imax,tstep,tstepratio,ratiotol,alpha)

%% function fixed boundary version 1
% e.g. PJCDtestPhi_fixedboundary_v1('plusnoisePhi6Fbrain.txt',65,20000,1e-7,1e-4,0,1)
%% consider only interior nodes
% 2015-05-26
% input:
% str:groundtruth transformation file
% N:grid size
% imax: max iteration step
% tstep: optimazition parameter
% tstepratio: desired tstep value
% ratiotol: desired ssd decreased ratio
% alpha: weighting parameter
% if (nargin==6)
%     alpha=1;
% end

% given a transformation T0
% compute f0 as the prescribed Jacobian determinant
% g0 as the curl field
% Try to find T=(T1,T2), such that
% ssd = 1/2*integral [(J(T(x))-f0(x))^2 +alpha*(curl(T)-g0)^2] dA minimized
% T=x+u
% Delta u =f as the monitor function
% partial ssd/partial f = g = (g1,g2)
% where Delta gi = nabla dot ai
% a1 = -P(T_2y,-T_2x)+ Q(0,1); a2 = -P(-T_1y,T_1x) - Q(1,0); 
%P = (J(T)-f0); Q= alpha*(curlT-g0) 
% tic;
% construct a uniform square grid on 
% Physican Domain [1,N]*[1,N],actual grid size N
h=1;
% interior grid
Npts=N-2;
% read jacobian and curl file
% f0=myfuncf0(str,N,h);
% g0=myfuncg0(str,N,h);
f0=Avg_J_T_1ab;
g0=Avg_Curl_Phi;
% g0=zeros(N);
% move boundary
f0=f0(2:Npts+1,2:Npts+1);
g0=g0(2:Npts+1,2:Npts+1);

%initial value for control function
f1=zeros(Npts);
f2=zeros(Npts);

%initial value for T=x
[T1,T2]=ndgrid(1:N,1:N);

%initial value for ssd
[T1y,T1x]=gradient(T1,h);
[T2y,T2x]=gradient(T2,h);
JT=T1x.*T2y-T1y.*T2x;
curlT=T2x-T1y;
JT_interior=JT(2:Npts+1,2:Npts+1);
curlT_interior=curlT(2:Npts+1,2:Npts+1);
T1y=T1y(2:Npts+1,2:Npts+1);
T1x=T1x(2:Npts+1,2:Npts+1);
T2y=T2y(2:Npts+1,2:Npts+1);
T2x=T2x(2:Npts+1,2:Npts+1);
%only compute the interior part
P=JT_interior-f0;
Q=curlT_interior-g0;
Q=Q*sqrt(alpha);
ssd_old=sum(sum(P.^2+Q.^2));
ssd_initial=ssd_old;

%iteration parameters
better=1;
iter=1;
ratio=1;

while iter<imax && tstep>tstepratio && ratio>ratiotol
    if better
        iter=iter+1;
        %ssd decreases
        %(b1,b2)=-P(T2y,-T2x)+ Q(0,1)
        % laplacian g1 = b1x+b2y
        b1=(-P.*T2y);
        b2=(P.*T2x+Q);
        [~,b1x]=gradient(b1,h);
        [b2y,~]=gradient(b2,h);
        a1=b1x+b2y;
        
        g1=pois2fft2(a1);

        %(c1,c2)=P(-T1y,T1x)-Q(1,0)
        % laplacian g2=c1x+c2y
        c1=P.*T1y-Q;
        c2=-P.*T1x;
        [~,c1x]=gradient(c1,h);
        [c2y,~]=gradient(c2,h);
        a2=c1x+c2y;
        g2=pois2fft2(a2);
    end
    %update f1n,f2n to calculate u1,u2
    f1n=f1-tstep*g1;
    f2n=f2-tstep*g2;
    %update u1,u2
    u1=pois2fft2(f1n);
    u2=pois2fft2(f2n);
    %update T1_new,T2_new to compute ssd
    T1_new=T1;
    T2_new=T2;
    T1_new(2:Npts+1,2:Npts+1)=T1_new(2:Npts+1,2:Npts+1)+u1;
    T2_new(2:Npts+1,2:Npts+1)=T2_new(2:Npts+1,2:Npts+1)+u2;
    %update JT,curlT
    [T1y,T1x]=gradient(T1_new,h);
    [T2y,T2x]=gradient(T2_new,h);
    JT=T1x.*T2y-T1y.*T2x;
    curlT=T2x-T1y;
    curlT_interior=curlT(2:Npts+1,2:Npts+1);
    JT_interior=JT(2:Npts+1,2:Npts+1);
    T1y=T1y(2:Npts+1,2:Npts+1);
    T1x=T1x(2:Npts+1,2:Npts+1);
    T2y=T2y(2:Npts+1,2:Npts+1);
    T2x=T2x(2:Npts+1,2:Npts+1);
    P=JT_interior-f0;
    Q=curlT_interior-g0;
    Q=Q*sqrt(alpha);
    ssd_new=sum(sum(P.^2+Q.^2));
    ratio=ssd_new/ssd_initial;
    % if ssd decrease, then update f,T, increase tstep
    % if not, decrese tstep
    if (ssd_new<ssd_old)
        tstep=tstep*1.1;
        f1=f1n;
        f2=f2n;
        ssd_old=ssd_new;
        %display([' tstep:',num2str(tstep),' ssd:',num2str(ssd_old)]);
        better=1;
    else
        better=0;
        tstep=tstep*2/3;
    end 
end
T1_new=T1;T2_new=T2;
T1_new(2:Npts+1,2:Npts+1)=T1_new(2:Npts+1,2:Npts+1)+u1;
T2_new(2:Npts+1,2:Npts+1)=T2_new(2:Npts+1,2:Npts+1)+u2;

% display results
% disp(['running time:',num2str(toc)]);
% disp(['ssd initial:',num2str(ssd_initial)]);
% disp(['ssd now:',num2str(ssd_old)]);
% disp([' avg_trans_ratio: ',num2str(ssd_old/ssd_initial)]);
% disp(['iteration step:',num2str(iter)]);

% gridplot2D3(str,N,T1_new,T2_new,2);

% s=input('save reconstructed T file for analysis ? (yes or no):','s');
% if (strcmp(s,'yes'))
%     strnew=input('save as file name:','s');
%     savePhi(T1_new,T2_new,strnew);
% end

% mygroundtruth(str,N,2);
% gridplot2D(T1_new,T2_new,2);

end
