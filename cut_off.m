function [x_new, y_new]=cut_off(x,y,N,theta,cute)
cut_off_bound=floor(cute*N);
r=sqrt((x-(N+1)*0.5)^2+(y-(N+1)*0.5)^2);
    if r>cut_off_bound
         x_new=x;
         y_new=y;
    else
        a=(cos((pi*r)/(cut_off_bound))+1)*0.5;
        u=[x-(N+1)*0.5; y-(N+1)*0.5];
        Alpha=a*theta;
        u=[cos(Alpha) -sin(Alpha); sin(Alpha) cos(Alpha)]*u;
        x_new=u(1)+(N+1)*0.5;
        y_new=u(2)+(N+1)*0.5;
    end
end