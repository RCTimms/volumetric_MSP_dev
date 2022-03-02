function [X,X_x,X_y,X_z]=get_timeseries(Dout,is_vector);
X=Dout.inv{1}.inverse.M*Dout.inv{1}.inverse.Y*Dout.inv{1}.inverse.T';

if is_vector==1;
    X_x=X(1:3:end,:);
    X_y=X(2:3:end,:);
    X_z=X(3:3:end,:);
else
    X_x=[];
    X_y=[];
    X_z=[];
end
end