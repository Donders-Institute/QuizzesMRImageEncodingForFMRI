function T2map=T2starMap(Magnitude,te)
% function T2map=T2starMap(Magnitude,te)
%This function takes for imput:
% Magnitude:a dataset where in the 1st and 2nd
%dimensions you have different pixels and in the 3rd dimension you have the
%different echo times
% values of the different echo times
    temp=0;
    for k=1:length(te)-1
        temp=temp+0.5*(Magnitude(:,:,k)+Magnitude(:,:,k+1))*(te(k+1)-te(k));
    end;
    % very fast estimation
    T2map=temp./(Magnitude(:,:,1)-Magnitude(:,:,end));
