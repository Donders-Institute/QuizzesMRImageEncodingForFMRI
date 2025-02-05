function [ FlipAngleDegrees ] = ErnstAngle(TR,T1)
%ERNSTANGLE computes the optimum flip angle for a given TR and T1
%    [ FlipAngleDegrees ] = ErnstAngle(TR,T1)


FlipAngleDegrees=acosd(exp(-TR/T1));


end

