
function simApproachSteadyState(T1,TR);
%%
% function simApproachSteadyState(T1,TR);
% T1 should be a T1 value
% TR should be the value of the repetition time used
%
% this function plots the signal of a GRE sequence for that T1 value run with repetition time 
% TR as a function time... to show the rate at which Steady state is
% achieved



%% checking that all the input variables exist
if ~exist('T1','var') || ~exist('TR','var')
    disp('ERROR: T1 and  TR must be specified')
    return
end

if ~exist('tColours','var')
    tColours = {[66 122 223]/255,[108 158 80]/255,[223 76 76]/255};
end

%% compute signal dependence as a function of T1
[ FlipAngleDegrees ] = round(ErnstAngle(TR,T1));

theta = [FlipAngleDegrees/2 FlipAngleDegrees 90];

time=0:TR:5*T1;

sig = (1-exp(-TR/T1)).*sin(theta*pi/180)./(1-cos(theta*pi/180)*exp(-TR/T1));

sigtime = zeros(length(time),length(theta));
Mz = ones(1,length(theta));

        for k=1:length(time)
        sigtime(k,:)=sind(theta).*Mz;
        % effect of flip angle
        Mz=cosd(theta).*Mz;
        %effect of decay of Mz
        Mz=bsxfun(@times,Mz,exp(-(TR)/T1));
        % effect of recovery
        Mz=bsxfun(@plus,Mz,1-exp(-(TR)/T1));
        end;



%% plot 
figure
set(gcf,'Color',[1 1 1]);

for k=1:length(theta)
plot(time,sigtime(:,k),'linewidth',2,'color',tColours{k}); 
hold on 
end
legend(['\theta= ' num2str(theta(1)) '°'],['\theta= ' num2str(theta(2)) 'Ernst'],['\theta= ' num2str(theta(3)) '°']);

for k=1:length(theta)
plot(linspace(time(1),time(end),100),repmat(sig(k),[100,1]),'-.','linewidth',2,'color',tColours{k}); 
end;

xlabel('time (ms)')
ylabel('Relative signal')
grid on
title(['{\bf Tissue:} T1 = ' num2str(T1) 'ms   {\bf Scan:} TR = ' num2str(TR) 'ms'])

fontScale(1.2)
xlim([0 time(end)])

