clc
clear
close all
%% Uniform grid
str='Preprocessed_MRI.png';
M=257;
N=65;% N=129;
m=N;n=N;h=4;
x=1:n;
y=1:m;
[X,Y]=ndgrid(x,y);
X_new1=X;Y_new1=Y;
X_new2=X;Y_new2=Y;
X_new3=X;Y_new3=Y;
X_new4=X;Y_new4=Y;
X_new5=X;Y_new5=Y;
X_new6=X;Y_new6=Y;
%% draw white pic
Blank_white=zeros(M)+255;
%% construct transformations
b1=0.47;
b2=0.47;
for i=1:n 
    for j=1:m
          %% shift x left
          [x_new1, y_new1]=cut_off(X(i,j),Y(i,j),N,-pi/8,b1);
          X_new1(i,j)=x_new1;
          Y_new1(i,j)=y_new1;
          %% shift y up
          [x_new2, y_new2]=cut_off(X(i,j),Y(i,j),N,pi/8,b1);
          X_new2(i,j)=x_new2;
          Y_new2(i,j)=y_new2;
          %% both x and y 
          [x_new3, y_new3]=shift_right_down(X(i,j),Y(i,j),N,b2);
          X_new3(i,j)=x_new3;
          Y_new3(i,j)=y_new3;
          %% both x and y 
          [x_new4, y_new4]=shift_Left_up(X(i,j),Y(i,j),N,b2);
          X_new4(i,j)=x_new4;
          Y_new4(i,j)=y_new4;
          %% both x and y 
          [x_new5, y_new5]=shift_right_up(X(i,j),Y(i,j),N,b2);
          [x_new5, y_new5]=cut_off(x_new5, y_new5,N,-pi/12,b1);
          X_new5(i,j)=x_new5;
          Y_new5(i,j)=y_new5;
          %% both x and y 
          [x_new6, y_new6]=shift_Left_down(X(i,j),Y(i,j),N,b2);
          [x_new6, y_new6]=cut_off(x_new6, y_new6,N,pi/12,b1);
          X_new6(i,j)=x_new6;
          Y_new6(i,j)=y_new6;
    end
end

%% loading Template images
I_0=double(rgb2gray(imread(str))); I_0=imresize(I_0,[m,n]);
I_1=interp2d2(X_new1,Y_new1, I_0);
I_2=interp2d2(X_new2,Y_new2, I_0);
I_3=interp2d2(X_new3,Y_new3, I_0);
I_4=interp2d2(X_new4,Y_new4, I_0);
I_5=interp2d2(X_new5,Y_new5, I_0);
I_6=interp2d2(X_new6,Y_new6, I_0);
%% testing avg imgs
h3=2;
figure(1)
imshow(imresize(I_1,[M,M]),[],'border', 'tight');
figure(2)
imshow(Blank_white,[],'border', 'tight'); hold on
for i = 1:h3:N
    plot(Y_new1(i,1:h3:end), X_new1(i,1:h3:end),'k-');hold on;
    plot(Y_new1(1:h3:end,i), X_new1(1:h3:end,i),'k-');hold on;
end
axis([1,N,1,N]);

figure(3)
imshow(imresize(I_2,[M,M]),[],'border', 'tight');
figure(4)
imshow(Blank_white,[],'border', 'tight'); hold on
for i = 1:h3:N
    plot(Y_new2(i,1:h3:end), X_new2(i,1:h3:end),'k-');hold on;
    plot(Y_new2(1:h3:end,i), X_new2(1:h3:end,i),'k-');hold on;
end
axis([1,N,1,N]);

figure(5)
imshow(imresize(I_3,[M,M]),[],'border', 'tight');
figure(6)
imshow(Blank_white,[],'border', 'tight'); hold on
for i = 1:h3:N
    plot(Y_new3(i,1:h3:end), X_new3(i,1:h3:end),'k-');hold on;
    plot(Y_new3(1:h3:end,i), X_new3(1:h3:end,i),'k-');hold on;
end
axis([1,N,1,N]);

figure(7)
imshow(imresize(I_4,[M,M]),[],'border', 'tight');
figure(8)
imshow(Blank_white,[],'border', 'tight'); hold on
for i = 1:h3:N
    plot(Y_new4(i,1:h3:end), X_new4(i,1:h3:end),'k-');hold on;
    plot(Y_new4(1:h3:end,i), X_new4(1:h3:end,i),'k-');hold on;
end
axis([1,N,1,N]);

figure(9)
imshow(imresize(I_5,[M,M]),[],'border', 'tight');
figure(10)
imshow(Blank_white,[],'border', 'tight'); hold on
for i = 1:h3:N
    plot(Y_new5(i,1:h3:end), X_new5(i,1:h3:end),'k-');hold on;
    plot(Y_new5(1:h3:end,i), X_new5(1:h3:end,i),'k-');hold on;
end
axis([1,N,1,N]);

figure(11)
imshow(imresize(I_6,[M,M]),[],'border', 'tight');
figure(12)
imshow(Blank_white,[],'border', 'tight'); hold on
for i = 1:h3:N
    plot(Y_new6(i,1:h3:end), X_new6(i,1:h3:end),'k-');hold on;
    plot(Y_new6(1:h3:end,i), X_new6(1:h3:end,i),'k-');hold on;
end
axis([1,N,1,N]);
%% initialization
k=1;kmax=4;
% ssd_initial=sum(sum((I_1-I_0).^2))+sum(sum((I_2-I_0).^2))+sum(sum((I_3-I_0).^2))+sum(sum((I_4-I_0).^2))+sum(sum((I_5-I_0).^2))++sum(sum((I_6-I_0).^2))
ratio=1;
ratioi=ones(kmax,1);
ssd=zeros(kmax,1);
ratiotol=1e-4;
I_1_temp=I_1;
I_2_temp=I_2;
I_3_temp=I_3;
I_4_temp=I_4;
I_5_temp=I_5;
I_6_temp=I_6;

ssd_initial_VS_GT(1)=sum(sum((I_1_temp-I_0).^2));
ssd_initial_VS_GT(2)=sum(sum((I_2_temp-I_0).^2));
ssd_initial_VS_GT(3)=sum(sum((I_3_temp-I_0).^2));
ssd_initial_VS_GT(4)=sum(sum((I_4_temp-I_0).^2));
ssd_initial_VS_GT(5)=sum(sum((I_5_temp-I_0).^2));
ssd_initial_VS_GT(6)=sum(sum((I_6_temp-I_0).^2));



%% main iteration

while (ratio>=ratiotol) && (k<=kmax)
%     tic; % % !!! Example !!!!
%     [Avg_diffeo_x,Avg_diffeo_y,I_0_avg]=avg_diffeo_gen_imgs_and_I0(I_0,I_1_temp,I_2_temp,I_3_temp,I_4_temp,I_5_temp,I_6_temp,N);
%     toc; 

    [I_1_avg,Avg1_diffeo_x,Avg1_diffeo_y]=avg_diffeo_gen_on_6_imgs(I_1_temp,I_2_temp,I_3_temp,I_4_temp,I_5_temp,I_6_temp,N);
    [I_2_avg,Avg2_diffeo_x,Avg2_diffeo_y]=avg_diffeo_gen_on_6_imgs(I_2_temp,I_1_temp,I_3_temp,I_4_temp,I_5_temp,I_6_temp,N);
    [I_3_avg,Avg3_diffeo_x,Avg3_diffeo_y]=avg_diffeo_gen_on_6_imgs(I_3_temp,I_1_temp,I_2_temp,I_4_temp,I_5_temp,I_6_temp,N);
    [I_4_avg,Avg4_diffeo_x,Avg4_diffeo_y]=avg_diffeo_gen_on_6_imgs(I_4_temp,I_1_temp,I_2_temp,I_3_temp,I_5_temp,I_6_temp,N);
    [I_5_avg,Avg5_diffeo_x,Avg5_diffeo_y]=avg_diffeo_gen_on_6_imgs(I_5_temp,I_1_temp,I_2_temp,I_3_temp,I_4_temp,I_6_temp,N);
    [I_6_avg,Avg6_diffeo_x,Avg6_diffeo_y]=avg_diffeo_gen_on_6_imgs(I_6_temp,I_1_temp,I_2_temp,I_3_temp,I_4_temp,I_5_temp,N);


    ssd_VS_GT(k,1)=sum(sum((I_1_avg-I_0).^2));
    ssd_VS_GT(k,2)=sum(sum((I_2_avg-I_0).^2));
    ssd_VS_GT(k,3)=sum(sum((I_3_avg-I_0).^2));
    ssd_VS_GT(k,4)=sum(sum((I_4_avg-I_0).^2));
    ssd_VS_GT(k,5)=sum(sum((I_5_avg-I_0).^2));
    ssd_VS_GT(k,6)=sum(sum((I_6_avg-I_0).^2))
    
    ratio_vs_gt(k,1)=ssd_VS_GT(k,1)/ssd_initial_VS_GT(1);
    ratio_vs_gt(k,2)=ssd_VS_GT(k,2)/ssd_initial_VS_GT(2);
    ratio_vs_gt(k,3)=ssd_VS_GT(k,3)/ssd_initial_VS_GT(3);
    ratio_vs_gt(k,4)=ssd_VS_GT(k,4)/ssd_initial_VS_GT(4);
    ratio_vs_gt(k,5)=ssd_VS_GT(k,5)/ssd_initial_VS_GT(5);
    ratio_vs_gt(k,6)=ssd_VS_GT(k,6)/ssd_initial_VS_GT(6);
%     ratioi=ssd_current(k)/ssd_initial
%     display(['ratio: ',num2str(ratio_all),' k: ',num2str(k)]);
    I_1_temp=I_1_avg;
    I_2_temp=I_2_avg;
    I_3_temp=I_3_avg;
    I_4_temp=I_4_avg;
    I_5_temp=I_5_avg;
    I_6_temp=I_6_avg;
    
    I_1_Avg(:,:,k)=I_1_avg;
    I_2_Avg(:,:,k)=I_2_avg;
    I_3_Avg(:,:,k)=I_3_avg;
    I_4_Avg(:,:,k)=I_4_avg;
    I_5_Avg(:,:,k)=I_5_avg;
    I_6_Avg(:,:,k)=I_6_avg;    
    
    Avg1X(:,:,k)=Avg1_diffeo_x;
    Avg2X(:,:,k)=Avg2_diffeo_x;
    Avg3X(:,:,k)=Avg3_diffeo_x;
    Avg4X(:,:,k)=Avg4_diffeo_x;
    Avg5X(:,:,k)=Avg5_diffeo_x;
    Avg6X(:,:,k)=Avg6_diffeo_x;   
    
    Avg1Y(:,:,k)=Avg1_diffeo_y;
    Avg2Y(:,:,k)=Avg2_diffeo_y;
    Avg3Y(:,:,k)=Avg3_diffeo_y;
    Avg4Y(:,:,k)=Avg4_diffeo_y;
    Avg5Y(:,:,k)=Avg5_diffeo_y;
    Avg6Y(:,:,k)=Avg6_diffeo_y;  
    
    k=k+1;
end


kk=3;
I_1_temp=squeeze(I_1_Avg(:,:,kk));
I_2_temp=squeeze(I_2_Avg(:,:,kk));
I_3_temp=squeeze(I_3_Avg(:,:,kk));
I_4_temp=squeeze(I_4_Avg(:,:,kk));
I_5_temp=squeeze(I_5_Avg(:,:,kk));
I_6_temp=squeeze(I_6_Avg(:,:,kk));

Avg1x=squeeze(Avg1X(:,:,kk));
Avg2x=squeeze(Avg2X(:,:,kk));
Avg3x=squeeze(Avg3X(:,:,kk));
Avg4x=squeeze(Avg4X(:,:,kk));
Avg5x=squeeze(Avg5X(:,:,kk));
Avg6x=squeeze(Avg6X(:,:,kk));

Avg1y=squeeze(Avg1Y(:,:,kk));
Avg2y=squeeze(Avg2Y(:,:,kk));
Avg3y=squeeze(Avg3Y(:,:,kk));
Avg4y=squeeze(Avg4Y(:,:,kk));
Avg5y=squeeze(Avg5Y(:,:,kk));
Avg6y=squeeze(Avg6Y(:,:,kk));
    
%% display 1
figure(111)
imshow(imresize(I_0,[M,M]),[],'border', 'tight');


figure(21)
imshow(imresize(I_1_temp,[M,M]),[],'border', 'tight');
figure(31)
imshow(Blank_white,[],'border', 'tight'); hold on
for i = 1:h3:N
    plot(Avg1y(i,1:h3:end), Avg1x(i,1:h3:end),'r-');hold on;
    plot(Avg1y(1:h3:end,i), Avg1x(1:h3:end,i),'r-');hold on;
end
axis([1,N,1,N]);

figure(22)
imshow(imresize(I_2_temp,[M,M]),[],'border', 'tight');
figure(32)
imshow(Blank_white,[],'border', 'tight'); hold on
for i = 1:h3:N
    plot(Avg2y(i,1:h3:end), Avg2x(i,1:h3:end),'r-');hold on;
    plot(Avg2y(1:h3:end,i), Avg2x(1:h3:end,i),'r-');hold on;
end
axis([1,N,1,N]);

figure(23)
imshow(imresize(I_3_temp,[M,M]),[],'border', 'tight');
figure(33)
imshow(Blank_white,[],'border', 'tight'); hold on
for i = 1:h3:N
    plot(Avg3y(i,1:h3:end), Avg3x(i,1:h3:end),'r-');hold on;
    plot(Avg3y(1:h3:end,i), Avg3x(1:h3:end,i),'r-');hold on;
end
axis([1,N,1,N]);

figure(24)
imshow(imresize(I_4_temp,[M,M]),[],'border', 'tight');
figure(34)
imshow(Blank_white,[],'border', 'tight'); hold on
for i = 1:h3:N
    plot(Avg4y(i,1:h3:end), Avg4x(i,1:h3:end),'r-');hold on;
    plot(Avg4y(1:h3:end,i), Avg4x(1:h3:end,i),'r-');hold on;
end
axis([1,N,1,N]);

figure(25)
imshow(imresize(I_5_temp,[M,M]),[],'border', 'tight');
figure(35)
imshow(Blank_white,[],'border', 'tight'); hold on
for i = 1:h3:N
    plot(Avg5y(i,1:h3:end), Avg5x(i,1:h3:end),'r-');hold on;
    plot(Avg5y(1:h3:end,i), Avg5x(1:h3:end,i),'r-');hold on;
end
axis([1,N,1,N]);

figure(26)
imshow(imresize(I_6_temp,[M,M]),[],'border', 'tight');
figure(36)
imshow(Blank_white,[],'border', 'tight'); hold on
for i = 1:h3:N
    plot(Avg6y(i,1:h3:end), Avg6x(i,1:h3:end),'r-');hold on;
    plot(Avg6y(1:h3:end,i), Avg6x(1:h3:end,i),'r-');hold on;
end
axis([1,N,1,N]);
%%
for i=1:kmax %stanard deviation od ssd
    std_v(i)=sum(ssd_VS_GT(i,1:6).^2)^(.5);
end
ssd_VS_GT(:,7)=std_v

for i=1:kmax %stanard deviation od ssd_ratio
    std_v_ratio(i)=sum(ratio_vs_gt(i,1:6).^2)^(.5);
end
ratio_vs_gt(:,7)=std_v_ratio


%% results
% Elapsed time is 94.174952 seconds.
% Elapsed time is 94.870533 seconds.
% Elapsed time is 93.522089 seconds.
% Elapsed time is 94.396645 seconds.
% Elapsed time is 98.700774 seconds.
% Elapsed time is 97.458275 seconds.
% ssd_VS_GT =
%    1.0e+06 *
%   Columns 1 through 3
%    1.227102679654169   1.357450150332655   2.280336032668330
%   Columns 4 through 6
%    1.321297557157216   2.890097417028171   2.905271278362469
% Elapsed time is 94.854763 seconds.
% Elapsed time is 101.548022 seconds.
% Elapsed time is 99.534126 seconds.
% Elapsed time is 99.975274 seconds.
% Elapsed time is 103.015870 seconds.
% Elapsed time is 100.922464 seconds.
% ssd_VS_GT =
%    1.0e+06 *
%   Columns 1 through 3
%    1.227102679654169   1.357450150332655   2.280336032668330
%    0.651471658149111   1.014394523951919   0.492021574621959
%   Columns 4 through 6
%    1.321297557157216   2.890097417028171   2.905271278362469
%    0.467789791836571   0.810710605141036   0.539820219536505
% Elapsed time is 100.702041 seconds.
% Elapsed time is 92.886979 seconds.
% Elapsed time is 102.551625 seconds.
% Elapsed time is 98.962992 seconds.
% Elapsed time is 102.632450 seconds.
% Elapsed time is 102.405514 seconds.
% ssd_VS_GT =
%    1.0e+06 *
%   Columns 1 through 3
%    1.227102679654169   1.357450150332655   2.280336032668330
%    0.651471658149111   1.014394523951919   0.492021574621959
%    0.506279308732094   0.517915196874645   0.505585687881774
%   Columns 4 through 6
%    1.321297557157216   2.890097417028171   2.905271278362469
%    0.467789791836571   0.810710605141036   0.539820219536505
%    0.481602262115011   0.599851944021664   0.546320927051006
% Elapsed time is 98.420260 seconds.
% Elapsed time is 99.787985 seconds.
% Elapsed time is 100.355420 seconds.
% Elapsed time is 102.369615 seconds.
% Elapsed time is 47.615442 seconds.
% Elapsed time is 100.261055 seconds.
% ssd_VS_GT =
%    1.0e+06 *
%   Columns 1 through 3
%    1.227102679654169   1.357450150332655   2.280336032668330
%    0.651471658149111   1.014394523951919   0.492021574621959
%    0.506279308732094   0.517915196874645   0.505585687881774
%    0.512137292095735   0.520902872177089   0.508610237470093
%   Columns 4 through 6
%    1.321297557157216   2.890097417028171   2.905271278362469
%    0.467789791836571   0.810710605141036   0.539820219536505
%    0.481602262115011   0.599851944021664   0.546320927051006
%    0.488787791325666   0.601814855073190   0.546986512342256
% ssd_VS_GT =
%    1.0e+06 *
%   Columns 1 through 3
%    1.227102679654169   1.357450150332655   2.280336032668330
%    0.651471658149111   1.014394523951919   0.492021574621959
%    0.506279308732094   0.517915196874645   0.505585687881774
%    0.512137292095735   0.520902872177089   0.508610237470093
%   Columns 4 through 6
%    1.321297557157216   2.890097417028171   2.905271278362469
%    0.467789791836571   0.810710605141036   0.539820219536505
%    0.481602262115011   0.599851944021664   0.546320927051006
%    0.488787791325666   0.601814855073190   0.546986512342256
%   Column 7
%    5.204562982064735
%    1.692034764874455
%    1.292440105517301
%    1.300998713018968
% ratio_vs_gt =
%   Columns 1 through 3
%    0.240169215712451   0.266562086605592   0.537035481297772
%    0.127506393548629   0.199196354193667   0.115874607667713
%    0.099089266550943   0.101702854819250   0.119069053568889
%    0.100235794297504   0.102289543739299   0.119781356677116
%   Columns 4 through 6
%    0.297105591736206   0.606792231829544   0.606781709666410
%    0.105186725093766   0.170213258059396   0.112744389194347
%    0.108292582766162   0.125942294448495   0.114102097319254
%    0.109908313376202   0.126354418680301   0.114241108427773
%   Column 7
%    1.114359140851585
%    0.349445403850051
%    0.273765641685640
%    0.275597978853802
