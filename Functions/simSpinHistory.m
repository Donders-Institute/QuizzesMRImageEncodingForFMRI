function [ExampleSignal]=simSpinHistory(zshiftPercentage,T1,TR,NSlices)
%%
% function [OptimumTheta MaxContrast]=simSpinHistory(zshiftPercentage,T1,TR,NSlices)
%
%
% zshiftPercentage   through slice shift
% T1                tissue T1
% TR                repetition time
% NSlices           NSlices
%
%
% simulates the Spin History effects of movement in the slice direction
% just before and after a given slice excitation
% For simplicity movement always hapens halfway through a 20 volume
% timecourse
%
%
% the ouptput from the function are the signals in two slices, one prior
% and one post movement


if ~exist('T1','var') || ~exist('TR','var')|| ~exist('zshiftPercentage','var')
    disp('ERROR: you seem to be missing one of these variables: zshiftPercentage,T1,TR,NSlices')
    return
end
% if ~or (strcmp(SliceMode,'interleaved'),strcmp(SliceMode,'sequential'));
%     disp('ERROR: SliceMode is either interleaved or sequential')
%     return
% end

if zshiftPercentage>100 || zshiftPercentage<-100
    disp('ERROR: for the purpose of this exercise you subject can only have a movement smaller between -100 100')
    return
end

if ~exist('tColours','var')
    tColours = {[66 122 223]/255,[108 158 80]/255,[223 76 76]/255};
end


%% doing all the calculations
% zshiftPercentage=50,T1=1.300,TR=3.000,NSlices=64;
SliceModeType{1}='interleaved';
SliceModeType{2}='sequential';
figure
set(gcf,'Position',[  69         200        1181         745]);
set(gcf,'Color',[1 1 1]);
counter=0;
for plusminus = [-1 1]
    zshiftPercentage = plusminus * abs(zshiftPercentage);
    for SliceModecounter = [1,2]
        counter=counter+1;
        SliceMode=SliceModeType{SliceModecounter};
        %initialize
        NVol=20;
        
        % SliceMode='interleaved';
        % SliceMode='sequential'
        if strcmp(SliceMode,'interleaved');
            order=(cat(2,mod(NSlices,2)+2:2:NSlices,mod(NSlices,2)+1:2:NSlices));
            [a order]=sort(order);
            
        elseif strcmp(SliceMode,'sequential');
            order=1:NSlices;
        end
        
        
        timing=order/NSlices*TR;
        
        timecourse = 1:NVol;
        
        NSlices_of_interest=4;
        timing=timing(1:NSlices_of_interest);
        order=order(1:NSlices_of_interest);
        
        % signal after the excitation
        % size timepoints x number of slices, subcells
        
        sig = ones(length(timecourse),NSlices_of_interest,100);
        
        % Mz after the excitationand recovery before the next volume
        % size timepoints x number of slices, subcells
        
        Mz = ones(length(timecourse),NSlices_of_interest,100);
        
        % finding Ernst angle
        [ FlipAngleDegrees ] = ErnstAngle(TR,T1);
        
        
        sig(1,:,:)=sind(FlipAngleDegrees);
        
        % effect of flip angle
        Mz(1,:,:)=cosd(FlipAngleDegrees);
        %effect of decay of Mz
        Mz(1,:,:)=bsxfun(@times,Mz(1,:,:),exp(-(TR-timing)/T1));
        % effect of recovery
        Mz(1,:,:)=bsxfun(@plus,Mz(1,:,:),1-exp(-(TR-timing)/T1));
        
        
        
        
        %% Driving into steady state
        
        
        for k=2:NVol/2
            %effect of decay of Mz
            Mz(k,:,:)=bsxfun(@times,Mz(k-1,:,:),exp(-(timing)/T1));
            % effect of recovery
            Mz(k,:,:)=bsxfun(@plus,Mz(k,:,:),1-exp(-(timing)/T1));
            
            sig(k,:,:)=sind(FlipAngleDegrees)*Mz(k,:,:);
            
            % effect of flip angle
            Mz(k,:,:)=cosd(FlipAngleDegrees).*Mz(k,:,:);
            %effect of decay of Mz
            Mz(k,:,:)=bsxfun(@times,Mz(k,:,:),exp(-(TR-timing)/T1));
            % effect of recovery
            Mz(k,:,:)=bsxfun(@plus,Mz(k,:,:),1-exp(-(TR-timing)/T1));
            
        end
        %%
        %MovementVolume
        
        
        midpoint=timing(order==2);
        k=k+1;
        
        temp=Mz(k,:,:);
        
        
        slicesbefore=find( or (order==1 , order==2) );
        slicesafter=find( and (order~=1 , order~=2) );
        
        
        %evolution until midpoint / first set of slices
        %effect of decay of Mz
        Mz(k,slicesbefore,:)=bsxfun(@times,Mz(k-1,slicesbefore,:),exp(-(timing(slicesbefore))/T1));
        % effect of recovery
        Mz(k,slicesbefore,:)=bsxfun(@plus,Mz(k,slicesbefore,:),1-exp(-(timing(slicesbefore))/T1));
        %signal
        sig(k,slicesbefore,:)=sind(FlipAngleDegrees)*Mz(k,slicesbefore,:);
        % effect of flip angle
        Mz(k,slicesbefore,:)=cosd(FlipAngleDegrees).*Mz(k,slicesbefore,:);
        %effect of decay of Mz
        Mz(k,slicesbefore,:)=bsxfun(@times,Mz(k,slicesbefore,:),exp(-(midpoint-timing(slicesbefore))/T1));
        % effect of recovery
        Mz(k,slicesbefore,:)=bsxfun(@plus,Mz(k,slicesbefore,:),1-exp(-(midpoint-timing(slicesbefore))/T1));
        
        
        %evolution until midpoint / second set of slices
        %effect of decay of Mz
        Mz(k,slicesafter,:)=bsxfun(@times,Mz(k-1,slicesafter,:),exp(-(midpoint)/T1));
        % effect of recovery
        Mz(k,slicesafter,:)=bsxfun(@plus,Mz(k,slicesafter,:),1-exp(-(midpoint)/T1));
        
        
        %%
        % introduce shift
        % figureJ(1);
        
        temp= squeeze(Mz(k,:,:))';
        % subplot(121)
        % imagesc(temp)
        temp=temp(:);
        temp=reshape(circshift(temp,[zshiftPercentage 0]),[100 NSlices_of_interest])';
        % subplot(122)
        % imagesc(temp)
        
        Mz(k,:,:)=temp;
        
        %%  go on with signal evolution
        
        %
        
        
        %evolution from midpoint / second set of slices
        %effect of decay of Mz
        Mz(k,slicesafter,:)=bsxfun(@times,Mz(k,slicesafter,:),exp(-(timing(slicesafter)-midpoint)/T1));
        % effect of recovery
        Mz(k,slicesafter,:)=bsxfun(@plus,Mz(k,slicesafter,:),1-exp(-(timing(slicesafter)-midpoint)/T1));
        %signal
        sig(k,slicesafter,:)=sind(FlipAngleDegrees)*Mz(k,slicesafter,:);
        % effect of flip angle
        Mz(k,slicesafter,:)=cosd(FlipAngleDegrees).*Mz(k,slicesafter,:);
        %effect of decay of Mz
        Mz(k,slicesafter,:)=bsxfun(@times,Mz(k,slicesafter,:),exp(-(TR-timing(slicesafter))/T1));
        % effect of recovery
        Mz(k,slicesafter,:)=bsxfun(@plus,Mz(k,slicesafter,:),1-exp(-(TR-timing(slicesafter))/T1));
        
        
        %evolution from midpoint / first set of slices
        %effect of decay of Mz
        Mz(k,slicesbefore,:)=bsxfun(@times,Mz(k-1,slicesbefore,:),exp(-(TR-midpoint)/T1));
        % effect of recovery
        Mz(k,slicesbefore,:)=bsxfun(@plus,Mz(k,slicesbefore,:),1-exp(-(TR-midpoint)/T1));
        %%
        movementpoint=k;
        for k=movementpoint+1:NVol
            %effect of decay of Mz
            Mz(k,:,:)=bsxfun(@times,Mz(k-1,:,:),exp(-(timing)/T1));
            % effect of recovery
            Mz(k,:,:)=bsxfun(@plus,Mz(k,:,:),1-exp(-(timing)/T1));
            
            sig(k,:,:)=sind(FlipAngleDegrees)*Mz(k,:,:);
            
            % effect of flip angle
            Mz(k,:,:)=cosd(FlipAngleDegrees).*Mz(k,:,:);
            %effect of decay of Mz
            Mz(k,:,:)=bsxfun(@times,Mz(k,:,:),exp(-(TR-timing)/T1));
            % effect of recovery
            Mz(k,:,:)=bsxfun(@plus,Mz(k,:,:),1-exp(-(TR-timing)/T1));
            
        end
        
        %% plotting the results
        
        subplot(3,2,2+counter)
        slice1=squeeze(sum(sig(:,slicesbefore(or (slicesbefore==2 , slicesbefore==3)),:),3)/100);
        slice2=squeeze(sum(sig(:,slicesafter(or (slicesafter==2 , slicesafter==3)),:),3)/100);
        % keyboard
        plot(1:NVol,slice1,'linewidth',2,'color',tColours{1}); hold all
        plot(1:NVol,slice2,'linewidth',2,'color',tColours{2});
        xlabel('Volume Number ');
        ylabel('Relative signal');
        grid on
        title(['through slice movement of ',num2str(zshiftPercentage),'%, Slice mode ' ,SliceMode ,', TR = ' num2str(TR) 'ms']);
if counter ==2
        legend(['last slice not to experience movement'], ['first slice to experience movement'],'Location','NorthOutside');
end;
        xlim([1 NVol])
        ylim([0 1])
        
        steadystate = slice1(movementpoint-1);
        variation1 = round(max(max(slice1(movementpoint:end))-steadystate,steadystate-min(slice1(movementpoint:end)))/steadystate*100);
        steadystate = slice2(movementpoint-1);
        variation2 = round(max(max(slice2(movementpoint:end))-steadystate,steadystate-min(slice2(movementpoint:end)))/steadystate*100);
        
        
        text(NVol/5 ,0.35,['maximum Signal Variation = ',num2str(variation1),'%']);
        text(NVol/5 ,0.15,['maximum Signal Variation = ',num2str(variation2),'%']);
        
    ExampleSignal{counter}=cat(2,slice1,slice2);
    end;
end;

        subplot(3,2,1)

        plot(1:NVol,cat(1,zeros(movementpoint,1),abs(zshiftPercentage)*ones(NVol-movementpoint,1)),'linewidth',2,'color',tColours{2});
        xlabel('Volume Number ');
        ylabel('Movement z');

fontScale(1.2)
