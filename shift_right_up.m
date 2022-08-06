function [x_new, y_new]=shift_right_up(x,y,N,b)
cut_off_bound=floor(0.4*N);
r=sqrt((x-(N+1)*0.5)^2+(y-(N+1)*0.5)^2);
    if r>cut_off_bound
         x_new=x;
         y_new=y;
    else
        a=(cos((pi*r)/(cut_off_bound))+1)*0.5*b;
        x_temp=(x)/N; 
        x_temp=(a*-(x_temp^2-x_temp)*0.5+x_temp)*N; 
        x_new=x_temp; 
        y_temp=(y)/N; 
        y_temp=(a*(y_temp^2-y_temp)*0.5+y_temp)*N; 
        y_new=y_temp; 
    end
end