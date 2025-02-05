
function MaxSignal = simSignalvFlip(T1,TR);
%%
% function simSignalvFlip(T1,TR);
% T1 should be a T1 value
% TR should be the value of the repetition time used
%
% this function plots the signal of a GRE sequence for that T1 value run with repetition time 
% TR as a function of the flip angle
% additionally it calculates the optimum flip angle also known as the ERNST
% angle



%% checking that all the input variables exist
if ~exist('T1','var') || ~exist('TR','var')
    disp('ERROR: T1 and  TR must be specified')
    return
end
%% compute signal dependence as a function of T1
theta = linspace(0,180,1000);
sig = (1-exp(-TR/T1)).*sin(theta*pi/180)./(1-cos(theta*pi/180)*exp(-TR/T1));
[MaxSignal position]=max(sig);
OptimumTheta=theta(position);

%% plot 
figure
set(gcf,'Color',[1 1 1]);

plot(theta,sig,'linewidth',2); 
xlabel('Flip Angle (\circ)')
ylabel('Relative signal')
grid on
title(['{\bf Tissue:} T1 = ' num2str(T1) 'ms   {\bf Scan:} TR = ' num2str(TR) 'ms'])
text(45 ,MaxSignal/2,['Optimum Flip=',num2str(round(OptimumTheta)),'degrees'])
fontScale(1.2)
xlim([0 180])

