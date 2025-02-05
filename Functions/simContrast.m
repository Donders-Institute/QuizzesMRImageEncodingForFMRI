function simContrast(useInversion,TR,TE,flipAngle,fieldStrength,varargin)
% function simContrast(useInversion,TR,TE,flipAngle,fieldStrength,varargin)
% this function simulates the contrast observed in the brain when using a
% GRE sequence with a predefined set of parameters 
% useInversion defines if the sequence is an inversion recovery sequence with one excitation at a time TI or
% standard gradient echo
% TR TE are the repetition times in ms
% flip angle is the flip angle used in degrees
% fieldStrength can be 1.5 3 or 7 in which case the respective relaxation
% times ares used
% function written by Daniel Gallichan Daniel.Gallichan@epfl.ch
% function modified by Jose Marques Jose.Marques@epfl.ch



%% checking that all the parameters exist
if ~exist('TR','var') || ~exist('flipAngle','var') || ~exist('TE','var')
    disp('ERROR: TR, flipAngle and TE must all be defined')
    return
end
if nargin>=6
    TI=varargin{1};
end;
    % 1- WM
% 2- GM
% 3- CSF
if ~exist('fieldStrength','var'),    fieldStrength=1.5; end
if fieldStrength==3
    T1vals = [830 1330 4000];
    T2starvals = [53 66 100];
else
    if fieldStrength==1.5;
    T1vals = [600 900 3500]; % T1 in ms
    T2starvals = [66 84 100];  % T2 in ms    
    else
    fieldStrength=7;
    T1vals = [1100 1800 4000]; % T1 in ms
    T2starvals = [26 33 100];  % T2 in ms    
        
    end
end

if ~exist('baseSNR','var'), baseSNR = 200; end;

% tColours = {[0 0 1],[0 .1 0],[1 0 0]};
tColours = {[66 122 223]/255,[108 158 80]/255,[223 76 76]/255};

load segBrain

if ~exist('useInversion','var'), useInversion = 0; end
% useInversion = 1;
if (useInversion) && ~exist('TI','var')
    disp(['ERROR: if inversion is chosen, TI must be set you have to add one more argument to the function']);
    return
end
% TI = 600;
%% checking it the parameters don't have comflicts

flip = flipAngle;

maxTE = min(TR,200);

if TE>maxTE
    disp(['ERROR: TE may not be greater than TR or 200ms (whichever is smaller)'])
    return
end

if useInversion && (TI>=TR)
    disp(['ERROR: TI may not be greater than (or equal to) TR'])
    return
end

% if ~exist('hSimContrast','var')
    hSimContrast = figure;
    % set(gcf,'Position',[    50   164   592   945])
     set(gcf,'Position',[    50   512   937   597])
%     set(gcf,'Position',[581   103   975   712])
    set(gcf,'Color',[1 1 1]);

% axHeight = .22; topY = .73; midY = .4; botY = .05; thinWidth = .08;
 axLeft = .12; mainWidth = .6876; loff = .1;
 axHeight = .35; topY = .58; midY = .1; thinWidth = .05; axLeft = .09; mainWidth=.5;loff=.02;
 subT1 = subplot('Position',[    axLeft    topY    mainWidth    axHeight]);
hold on; grid on; xlabel('Time since last excitation pulse (ms)')
if useInversion
    title(['{\bf T1-weighting:} TR = ' num2str(TR) 'ms, Flip Angle = ' num2str(flip) ' degrees, TI = ' num2str(TI) 'ms'])
else
    title(['{\bf T1-weighting:} TR = ' num2str(TR) 'ms, Flip Angle = ' num2str(flip) ' degrees'])
end
ylabel('Relative M_z')
    
subT1w = axes('Position',[    mainWidth+axLeft+loff   topY    thinWidth    axHeight]);
set(gca,'xtick',[],'ytick',[],'box','off')
hold on; axis([0 2 0 1]);
title('T1-weighting')

mainWidthT2 = .3;
axLeftT2 = .29;
subT2 = axes('Position',[   axLeftT2    midY   mainWidthT2    axHeight]);
hold on; grid on; title(['{\bf T2-weighting:} TE = ' num2str(TE) 'ms'])
xlabel('Time since excitation pulse (ms)')
ylabel('Relative M_{xy}')
axis([0 maxTE 0 1.05])
subT2w = subplot('Position',[    mainWidthT2+axLeftT2+loff    midY    thinWidth    axHeight]);
hold on; axis([0 2 0 1]);
set(gca,'xtick',[],'ytick',[],'box','off')
title('T2-weighting')

subText = axes('Position',[    0.0164    0.9397    0.1010    0.0457]);
axis off
axis([0 1 0 1])
text(0,1,['{\bfField Strength:}' num2str(fieldStrength) 'T'])

% 
subT1T2w = subplot('Position',[ 0.7256    0.6388    0.1987    0.2381]);


hold on; axis([0 2 0 1.1]); title({'\bf Combined contrast','(normalised)'})
set(gca,'xtick',[.5 1 1.5],'ytick',[],'xticklabel',{'WM','GM','CSF'})
ylabel('Relative Signal')


subIm = axes('Position',[        0.7028    0.0573    0.2517    0.4459]);

subImSeg = axes('Position',[       0.0653    0.0999    0.1535    0.3640]);
imSeg = zeros(size(segBrain,1),size(segBrain,2),3);
for iC = 1:3
    [inX inY] = find(segBrain==iC);
    for iX = 1:length(inX)
        imSeg(inX(iX),inY(iX),1) = tColours{iC}(1);
        imSeg(inX(iX),inY(iX),2) = tColours{iC}(2);
        imSeg(inX(iX),inY(iX),3) = tColours{iC}(3);
    end
end
image(imSeg); axis equal tight off

lw = 10;

if useInversion
    M_minus = (1 - exp(-(TR-TI)./T1vals) + (1 - exp(-TI./T1vals)).*cos(flip*pi/180).*exp(-(TR-TI)./T1vals))./(1 + exp(-TR./T1vals)*cos(flip*pi/180));
    MTI_minus = ((1 - exp(-TI./T1vals)) - M_minus.*exp(-TI./T1vals));
    sigT1 = abs(MTI_minus*sin(flip*pi/180));
    t = linspace(0,TR,500);
    iTI = round(interp1(t,1:length(t),TI));
    M = zeros(length(t),length(T1vals));
    for iT1 = 1:length(T1vals)
        M(1:iTI,iT1) = (1 - exp(-t(1:iTI)/T1vals(iT1))) - M_minus(iT1)*exp(-t(1:iTI)/T1vals(iT1));
        M(iTI+1:end,iT1) = (1 - exp(-(t(iTI+1:end)-TI)/T1vals(iT1))) + MTI_minus(iT1)*cos(flip*pi/180)*exp(-(t(iTI+1:end)-TI)/T1vals(iT1));
    end
 
else
    M_minus = (1-exp(-TR./T1vals))./(1-cos(flip*pi/180)*exp(-TR./T1vals));
    sigT1 = M_minus*sin(flip*pi/180);
    t = linspace(0,TR,500);
    M = zeros(length(t),length(T1vals));
    for iT1 = 1:length(T1vals)
        M(:,iT1) = 1 - exp(-t/T1vals(iT1)) + M_minus(iT1)*cos(flip*pi/180)*exp(-t/T1vals(iT1));
    end


end
tt = [-t(125:-1:1) t  t(1:125)+TR];
for iLine = 1:length(T1vals)
    subplot(subT1)
    plot(tt,[M(end-124:end,iLine); M(:,iLine); M(1:125,iLine)],'linewidth',2,'color',tColours{iLine});
    axis([tt(1) tt(end) min( min(M(:)) ,0) 1.05*max(M(:))])
        
    subplot(subT1w)
    line([0 0]+iLine*.5,[0 sigT1(iLine)],'linewidth',lw,'color',tColours{iLine})
%     subplot(subT1w2)
%     line([0 0]+iLine*.5,[0 sigT1(iLine)],'linewidth',lw,'color',tColours{iLine})
end
subplot(subT1)

if ~useInversion
    line([0 0], [-1 1], 'linewidth',2,'color','k','linestyle','--')
    line([TR TR], [-1 1], 'linewidth',2,'color','k','linestyle','--')
else
    line([0 0], [-1 1], 'linewidth',4,'color','k','linestyle',':')
    line([TR TR], [-1 1], 'linewidth',4,'color','k','linestyle',':')
    line([TI TI], [-1 1], 'linewidth',2,'color','k','linestyle',':')
    xlabel('Time since last inversion pulse (ms)')
end


t = linspace(0,maxTE,100);
sigT2 = exp(-TE./T2starvals);
sigT1T2 = sigT1.*sigT2;
subplot(subT2)
line([TE TE],[0 1],'linewidth',2,'color','k','linestyle','--')
for iLine = 1:length(T2starvals)
    subplot(subT2)
    plot(t,exp(-t/T2starvals(iLine)),'linewidth',2,'color',tColours{iLine}); 
    
    subplot(subT2w)
    line([0 0]+iLine*.5,[0 sigT2(iLine)],'linewidth',lw,'color',tColours{iLine})
       
    subplot(subT1T2w)
    line([0 0]+iLine*.5,[0 sigT1T2(iLine)/max(sigT1T2)],'linewidth',20,'color',tColours{iLine})
end

fontScale(1.2)
bIm = zeros(size(segBrain));
for iLine = 1:length(T1vals)
    bIm(segBrain==iLine) = sigT1T2(iLine)/max(sigT1T2);
end
subplot(subIm)
noiseIm = randn(size(segBrain))/(max(sigT1T2)*baseSNR);
% noiseIm = randn(size(segBrain))/baseSNR;
% imagesc(abs(bIm+noiseIm),[0 1]); 
imagesc(abs(bIm+noiseIm)); 
axis equal tight off; colormap(gray)
% title({['rel. GM/WM contrast: ' abs(num2str(abs(sigT1T2(2)-sigT1T2(1))/sigT1T2(1),3))],...
title({    ['{\bf Total acq. time:} ' num2str(round(prod(size(bIm))*TR/1000/6)/10) ' mins'],...
    ['{\bf WM SNR:} ' num2str(sigT1T2(1)*baseSNR,3) ... num2str(max(sigT1T2)*200,3) ...
    ' {\bf GM/WM CNR}: ' num2str(abs((sigT1T2(2)-sigT1T2(1))*baseSNR),3)] })



