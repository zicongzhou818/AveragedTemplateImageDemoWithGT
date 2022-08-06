function [Ri] = interp2d2(pos_x_total,pos_y_total,R)
%linear interpolation
%(i,j)--->(pos_x(i,j),pos_y(i,j))
[m,n] = size(R);
Ri=zeros(m,n);
for i=1:m
    for j=1:n
        x=pos_x_total(i,j);
        y=pos_y_total(i,j);
        if (x>= m)
            x=m;
            ceil_x=m;
            floor_x=ceil_x-1;
        else
            if (x<=1)
                x=1;
            end
            floor_x=floor(x);
            ceil_x=floor_x+1;
        end
        if (y>= n)
            y=n;
            ceil_y=n;
            floor_y=ceil_y-1;
        else
            if (y<=1)
                y=1;
            end
            floor_y=floor(y);
            ceil_y=floor_y+1;
        end
        
      Ri(i,j)=R(ceil_x,ceil_y)*(x-floor_x)*(y-floor_y)+R(floor_x,floor_y)*(ceil_x-x)*(ceil_y-y)+R(floor_x,ceil_y)*(ceil_x-x)*(y-floor_y)+R(ceil_x,floor_y)*(x-floor_x)*(ceil_y-y);
    end
end

end