function [y_p] = def(y,v,a,c,d,nz,pd_1,pd_2)
y_1=zeros(1,pd_2-pd_1+1);
y_n=zeros(1,pd_2-pd_1+1);
y=[y_1;y;y_n];
y_p=zeros(((c-a)/d)+1,pd_2-pd_1+1);
y_p(1:end,:)=(y(3:end,:)-...
    y(1:end-2,:))./(2*d);
end

