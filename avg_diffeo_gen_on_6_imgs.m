function [I_avg_Template,Avg_diffeo_x,Avg_diffeo_y]=avg_diffeo_gen_on_6_imgs(Template,I_a,I_b,I_c,I_d,I_e,N)
%% register a to b and  register b to a
[pos_xT2a, pos_yT2a]=Img_Reg(Template, I_a,N);
[pos_xT2b, pos_yT2b]=Img_Reg(Template, I_b,N);
[pos_xT2c, pos_yT2c]=Img_Reg(Template, I_c,N);
[pos_xT2d, pos_yT2d]=Img_Reg(Template, I_d,N);
[pos_xT2e, pos_yT2e]=Img_Reg(Template, I_e,N);
%% average jacobian determinant and average curl-vector
[J_T2a, Curl_T2a]=compute_JD_and_Curl(pos_xT2a, pos_yT2a,N,1);
[J_T2b, Curl_T2b]=compute_JD_and_Curl(pos_xT2b, pos_yT2b,N,1);
[J_T2c, Curl_T2c]=compute_JD_and_Curl(pos_xT2c, pos_yT2c,N,1);
[J_T2d, Curl_T2d]=compute_JD_and_Curl(pos_xT2d, pos_yT2d,N,1);
[J_T2e, Curl_T2e]=compute_JD_and_Curl(pos_xT2e, pos_yT2e,N,1);
Avg_J=(J_T2a+J_T2b+J_T2c+J_T2d+J_T2e+1)./6;
Avg_Curl=(Curl_T2a+Curl_T2b+Curl_T2c+Curl_T2d+Curl_T2e)./6;
%% reconstruct the averaged diffeomorphism a2b and b2a
[Avg_diffeo_x,Avg_diffeo_y]=PJCDtestPhi_fixedboundary_v1(Avg_J, Avg_Curl, N, 4000, 1e-4, 1e-16,1e-6,1);
I_avg_Template=interp2d2(Avg_diffeo_x,Avg_diffeo_y, Template);
end




