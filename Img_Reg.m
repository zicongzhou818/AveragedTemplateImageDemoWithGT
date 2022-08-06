function [pos_x,pos_y]=Img_Reg(T,R,N)
m=N;n=N;
%% Image registration 
%% step-0
% read images to T as template or moving image, R as reference or fixed image
% T=double(rgb2gray(imread(str)));T=imresize(T,[m,n]); 
% R=I_temp;
% R=double(imread('T_1.png'));R=imresize(R,[m,n]); % Reference 1

% register T to R and minimize ssd=(T(Phi)-R)^2 on Physical computation domain: [1,m]X[1,n]
% set up optimization parameters
tstep = 1e-9;% tstep: initial tstep 
tstep_ratio = 1e-32;% tstep_ratio: stop when tstep less than tstep_ratio
tstep_up=1.15;% tstep_up: magnification factor when ssd decreases
tstep_down=0.85;% tstep_down: reduction factor when doesn't ssd decrease
ite_max=1e4;% ite_max: preset maximum number of iterations
ratiotol=1e-5;% ratiotol: preset ratio to stop 
ratio=1;

%% Step-1
grid_size=[m-2,n-2];% remove zero boundaries
% set up initial values for f=(f1,f2)=0
f1=zeros(grid_size);
f2=zeros(grid_size);

% Compute initial ssd on uniform grid
ssd_old=sum(sum(((T-R).^2)));
ssd_initial=ssd_old;

% iteration flag
Better = 1;% Better: indicate whether the ssd decreases(Better=1) or not 
ii=0;% ii: total iteration steps
iii = 0;% iii: iteration steps that ssd decreases

% set up initial values for Phi(x,y)=(x+u1(x,y),y+u2(x,y))
[xI,yI] = ndgrid(1:m,1:n);% (i,j)--->(pos_x(i,j),pos_y(i,j)) on (1 to m,1 to n)
pos_x=xI;
pos_y=yI;

%% Step-9 Main loop
while  (tstep>tstep_ratio) && (ii<ite_max) && (ratio>ratiotol) 
    ii=ii+1; % record the number of iterations
    %% Step-2
    if Better  % --(PS: Better holds for the first step in "while" loop)
        iii=iii+1;
        % compute the difference from R to T (ps: not the same as from T to R)
        T_temp = interp2d2(pos_x,pos_y,T); 
        TR_diff = T_temp-R; 
        [T_x2, T_x1] = gradient(T); 
        T_x1 = interp2d2(pos_x,pos_y, T_x1); 
        T_x2 = interp2d2(pos_x,pos_y, T_x2); 
        % construct  lap(a)=[R(phi(X))-T(X)]*grad(R(phi(X)))
        lap_a1 = (TR_diff.*T_x1);
        lap_a2 = (TR_diff.*T_x2);
        % remove boundaries of a, since a=0 on boundaries 
        lap_a1_interior=lap_a1(2:m-1, 2:n-1);
        lap_a2_interior=lap_a2(2:m-1, 2:n-1);
        % w1,w2 
        w1=pois2fft2(lap_a1_interior);
        w2=pois2fft2(lap_a2_interior);
    end
    %% Step-3 Update F
    f1_new=f1-(w1)*tstep;
    f2_new=f2-(w2)*tstep;
    
    %% Step-4 Solve for U, where lap(U)=F
    intU1=pois2fft2(f1_new);
    intU2=pois2fft2(f2_new);
    % add boundaries back to get U, where U=0 on boundaries
    U1=matrixpad(intU1,0);
    U2=matrixpad(intU2,0); 

    %% Step-6  % phi(x,y)=X+U
    pos_x=xI+U1;
    pos_y=yI+U2;
    %% Step-7
    T_temp=interp2d2(pos_x,pos_y,T); 
    %% Step-8 Key part for gradient descent
    ssd=sum(sum(((T_temp-R).^2))); % use sums to approximate integrals
    ratio=ssd/ssd_initial; 
%     display([' tstep: ',num2str(tstep),' ssd: ',num2str(ssd),' ratio: ',num2str(ratio),' ii: ',num2str(ii)]);
    if ssd>ssd_old
        tstep=tstep*tstep_down;% reduction factor when ssd doesn't decrease
        Better = 0;
    else
        tstep=tstep*tstep_up;% magnification factor when ssd decreases
        ssd_old=ssd;
        Better = 1;
        f1 = f1_new;
        f2 = f2_new;
    end    
end
% display([' Img-Reg_ratio: ',num2str(ratio)]);
