function SeqParamOut = writeEPIRS_HandsOn(SeqParam,sys,Show)
% Help function for writeEPIRS_HandsOn
% 
% This function generates an experimental high-performance EPI sequence
% using split gradients to overlap blips with the readout gradients combined with ramp-sampling.
% based on the Pulseq EPI example written by MAxim Zaitsev
%
% Usage:
%   SeqParamOut = writeEPIRS_HandsOn(SeqParam, sys, Show)
%
% Inputs:
%   SeqParam - Structure containing sequence parameters with the following possible fields:
%       Nx              - Number of readout points (default: 64)
%       Ny              - Number of phase encoding steps (default: Nx)
%       Nz              - Number of slices (default: 4)
%       readoutTime     - Readout time in seconds (default: 4.2e-4)
%       TR              - Repetition time in seconds (default: 2000e-3)
%       TE              - Echo time in seconds (default: 30e-3)
%       ME              - Number of echoes (default: 1)
%       FOV             - Field of view in meters (default: 256e-3)
%       sliceThickness  - Slice thickness in meters (default: 4e-3)
%       sliceGap        - Gap between slices in meters (default: 0)
%       flipAngle       - Flip angle in degrees (default: 90)
%       timeBwProduct   - Time-bandwidth product (default: 4)
%       iPat            - Parallel imaging factor (default: 1)
%
%   sys - Structure containing system limits obtained from mr.opts from Pulseq.
%
%   Show - Optional structure specifying which plots to output with the following fields:
%       SequenceDiagram - Plot the sequence diagram (default: 1)
%       SliceProfile    - Plot the slice profile (default: 1)
%       KspaceTimeCourse- Plot the k-space time course (default: 1)
%       KspaceCoverage  - Plot the k-space coverage (default: 1)
%       Sound  - plays the sequence sound (default: 0)
%
% Outputs:
%   SeqParamOut - Updated SeqParam structure with any missing fields filled with default values.

% Check if SeqParam fields exist, otherwise use default values
if ~isfield(SeqParam, 'Nx')
    warning('SeqParam.Nx is missing. Using default value 64.');
    SeqParam.Nx = 64;
end
if ~isfield(SeqParam, 'Ny')
    warning('SeqParam.Ny is missing. Using default value = SeqParam.Nx.');
    SeqParam.Ny = SeqParam.Nx;
end
if ~isfield(SeqParam, 'Nz')
    warning('SeqParam.Nz is missing. Using default value 4.');
    SeqParam.Nz = 4;
end
if ~isfield(SeqParam, 'readoutTime')
    warning('SeqParam.readoutTime is missing. Using default value 4.2e-4.');
    SeqParam.readoutTime = 4.2e-4;
end
if ~isfield(SeqParam, 'TR')
    warning('SeqParam.TR is missing. Using default value 2000e-3.');
    SeqParam.TR = 2000e-3;
end
if ~isfield(SeqParam, 'TE')
    warning('SeqParam.TE is missing. Using default value 30e-3.');
    SeqParam.TE = 30e-3;
end
if ~isfield(SeqParam, 'ME')
    warning('SeqParam.ME is missing. Using default value 1 - single Echo.');
    SeqParam.ME = 1;
end
if length(SeqParam.TE) ~= SeqParam.ME
    warning('SeqParam.TE length was not consistent with number of echos - ');
    SeqParam.TE = repmat(SeqParam.TE(1),SeqParam.ME,1) + (0:(SeqParam.ME-1))*1e-3; % default TE values only with 1 ms increments are temporarily given 
end
if ~isfield(SeqParam, 'FOV')
    warning('SeqParam.FOV is missing. Using default value 256e-3.');
    SeqParam.FOV = 256e-3;
end
if ~isfield(SeqParam, 'sliceThickness')
    warning('SeqParam.sliceThickness is missing. Using default value 4e-3.');
    SeqParam.sliceThickness = 4e-3;
end
if ~isfield(SeqParam, 'sliceGap')
    warning('SeqParam.sliceGap is missing. Using default value 0.');
    SeqParam.sliceGap = 0;
end
if ~isfield(SeqParam, 'flipAngle')
    warning('SeqParam.flipAngle is missing. Using default value 90.');
    SeqParam.flipAngle = 90;
end
if ~isfield(SeqParam, 'timeBwProduct')
    warning('SeqParam.timeBwProduct is missing. Using default value 4.');
    SeqParam.timeBwProduct = 4;
    Show.SliceProfile = 0;
end
if ~isfield(SeqParam, 'iPat')
    warning('SeqParam.iPat is missing. Using default value 1 - no parallel imaging.');
    SeqParam.iPat = 1;
end



% Assign SeqParam values to local variables
Nx = SeqParam.Nx;
Ny = SeqParam.Ny;
Nslices = SeqParam.Nz;
readoutTime = SeqParam.readoutTime;
TR = SeqParam.TR;
TE = SeqParam.TE;
ME = SeqParam.ME;
fov = SeqParam.FOV;
thickness = SeqParam.sliceThickness;
sliceGap = SeqParam.sliceGap;
flipAngle = SeqParam.flipAngle;
timeBwProduct = SeqParam.timeBwProduct;	
iPat = SeqParam.iPat;


% Define figures to be shown
if nargin < 3
    Show.SequenceDiagram = 1;
    Show.SliceProfile = 1;
    Show.KspaceTimeCourse = 1;
    Show.KspaceCoverage = 1;
    Show.Sound = 0;
    else
    if ~isfield(Show, 'SequenceDiagram')
        Show.SequenceDiagram = 1;
    end
    if ~isfield(Show, 'SliceProfile')
        Show.SliceProfile = 1;
    end
    if ~isfield(Show, 'KspaceTimeCourse')
        Show.KspaceTimeCourse = 1;
    end
    if ~isfield(Show, 'KspaceCoverage')
        Show.KspaceCoverage = 1;
    end
    if ~isfield(Show, 'Sound')
        Show.Sound = 0;
    end
end


seq=mr.Sequence(sys);      % Create a new sequence object

pe_enable=1;               % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
ro_os=1;                   % oversampling factor (in contrast to the product sequence we don't really need it)
partFourierFactor=1;       % partial Fourier factor: 1: full sampling 0: start with ky=0

Ny_eff=round(Ny/iPat /2)*2; % number of lines acquired in one excitation accounting for iPat
fov_y_eff=fov/iPat; % reduced FOV in y-direction

% Create fat-sat pulse 
sat_ppm=-3.45; % fat saturation pulse parameters
sat_freq=sat_ppm*1e-6*sys.B0*sys.gamma; % saturation pulse frequency offset in Hz
rf_fs = mr.makeGaussPulse(110*pi/180,'system',sys,'Duration',8e-3,'dwell',10e-6,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
rf_fs.phaseOffset = -2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase    
rf_fs.name = 'fat-sat'; % useful for debugging, can be seen in seq.plot
gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',0.1/1e-4); % spoil up to 0.1mm

% Create excitation slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(flipAngle,'system',sys,'Duration',4e-3,...
    'SliceThickness',thickness,'apodization',0.42,'timeBwProduct',timeBwProduct,'use','excitation');
rf.name = 'rf90'; % useful for debugging, can be seen in seq.plot

% define the output trigger to play out with every slice excitatuion
trig=mr.makeDigitalOutputPulse('osc0','duration', 100e-6); % possible channels: 'osc0','osc1','ext1'

% Define other gradients and ADC events
deltak=1/fov;
deltak_y=1/fov_y_eff;
kWidth = Nx*deltak;

% Phase blip in shortest possible time
blip_dur = ceil(2*sqrt(deltak_y/sys.maxSlew)/10e-6/2)*10e-6*2; % we round-up the duration to 2x the gradient raster time
% the split code below fails if this really makes a trapezoid instead of a triangle...
gy = mr.makeTrapezoid('y',sys,'Area',-deltak_y,'Duration',blip_dur); % we use negative blips to save one k-space line on our way towards the k-space center


% readout gradient is a truncated trapezoid with dead times at the beginnig
% and at the end each equal to a half of blip_dur
% the area between the blips should be defined by kWidth
% we do a two-step calculation: we first increase the area assuming maximum
% slewrate and then scale down the amlitude to fix the area 

extra_area=blip_dur/2*blip_dur/2*sys.maxSlew; % check unit!;
gx = mr.makeTrapezoid('x',sys,'Area',kWidth+extra_area,'duration',readoutTime+blip_dur);
actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
gx.amplitude=gx.amplitude/actual_area*kWidth;
gx.area = gx.amplitude*(gx.flatTime + gx.riseTime/2 + gx.fallTime/2);
gx.flatArea = gx.amplitude*gx.flatTime;
gx.name='Gro'; % useful for debugging, can be seen in seq.plot

% calculate ADC
% we use ramp sampling, so we have to calculate the dwell time and the
% number of samples, which are will be qite different from Nx and
% readoutTime/Nx, respectively. 
adcDwellNyquist=deltak/gx.amplitude/ro_os;
% round-down dwell time to 100 ns
adcDwell=floor(adcDwellNyquist*1e7)*1e-7;
adcSamples=floor(readoutTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
% MZ: no idea, whether ceil,round or floor is better for the adcSamples...
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',blip_dur/2);
% realign the ADC with respect to the gradient
time_to_center=adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)*1e6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us 
% this rounding actually makes the sampling points on odd and even readouts
% to appear misalligned. However, on the real hardware this misalignment is
% much stronger anyways due to the grdient delays

% FOV positioning requires alignment to grad. raster... -> TODO

% split the blip into two halves and produce a combined synthetic gradient
gy_parts = mr.splitGradientAt(gy, blip_dur/2, sys);
[gy_blipup, gy_blipdown,~]=mr.align('right',gy_parts(1),'left',gy_parts(2),gx);
gy_blipdownup=mr.addGradients({gy_blipdown, gy_blipup}, sys);

% pe_enable support
gy_blipup.waveform=gy_blipup.waveform*pe_enable;
gy_blipdown.waveform=gy_blipdown.waveform*pe_enable;
gy_blipdownup.waveform=gy_blipdownup.waveform*pe_enable;

% phase encoding and partial Fourier

Ny_pre=round(partFourierFactor*Ny_eff   /2-1); % PE steps prior to ky=0, excluding the central line
Ny_post=round(Ny_eff/2+1); % PE lines after the k-space center including the central line
Ny_meas=Ny_pre+Ny_post;

% Pre-phasing gradients
gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2);
gyPre = mr.makeTrapezoid('y',sys,'Area',Ny_pre*deltak_y);
[gxPre,gyPre,gzReph]=mr.align('right',gxPre,'left',gyPre,gzReph);
% relax the PE prepahser to reduce stimulation
gyPre = mr.makeTrapezoid('y',sys,'Area',gyPre.area,'Duration',mr.calcDuration(gxPre,gyPre,gzReph));
gyPre.amplitude=gyPre.amplitude*pe_enable;

% Pre-phasing gradients in case of Multi-echo acquisition
if mod(Ny, 2) == 0
    gxPre_ME = mr.makeTrapezoid('x',sys,'Area',0);
%    gxPre_ME = gx;
%    gxPre_ME.amplitude = 0*gx.amplitude
else
    gxPre_ME = mr.makeTrapezoid('x',sys,'Area',-gx.area);
%    gxPre_ME = gx;
%    gxPre_ME.amplitude = -gx.amplitude
end
gyPre_ME = mr.makeTrapezoid('y',sys,'Area',(Ny_pre+Ny_post-1)*deltak_y);
[gxPre_ME,gyPre_ME]=mr.align('right',gxPre_ME,'left',gyPre_ME);
% relax the PE prepahser to reduce stimulation
gyPre_ME = mr.makeTrapezoid('y',sys,'Area',gyPre_ME.area,'Duration',mr.calcDuration(gxPre_ME,gyPre_ME));


% slice positions
slicePositions=(thickness+sliceGap)*((0:(Nslices-1)) - (Nslices-1)/2);
slicePositions=slicePositions([1:2:Nslices 2:2:Nslices]); % reorder slices for an interleaved acquisition (optional)

% Define sequence blocks
%seq.addBlock(mr.makeDelay(1)); % older scanners like Trio may need this
                                % dummy delay to keep up with timing

minTR = 0;                           
minTR_prep = mr.calcDuration(gz_fs) + mr.calcDuration(gz) + mr.calcDuration(gxPre,gyPre,gzReph)     ;                          
minTE = 0;                           
minTE = mr.calcDuration(rf) / 2 + gz.fallTime + mr.calcDuration(gxPre,gyPre,gzReph)   ;                             

for i=1:Ny_meas
        if i==1
            minTR = minTR + mr.calcDuration(gx,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
            minTE = minTE + mr.calcDuration(gx,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
        elseif i==Ny_meas
            minTR = minTR + mr.calcDuration(gx,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
        else
            minTR = minTR + mr.calcDuration(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
            if i< Ny_pre %
                minTE = minTE + mr.calcDuration(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
            elseif i == Ny_pre
                minTE = minTE + 0.5 * mr.calcDuration(gx,gy_blipdownup,adc); % Read the midle line, just take the half of it 
            end; 
        end 
end;

minDTE = minTR + mr.calcDuration(gxPre_ME,gyPre_ME);
minTR = minTR_prep + ME * (minTR)  + (ME-1) * mr.calcDuration(gxPre_ME,gyPre_ME);



TEfill = TE(1) - minTE; 
SeqParamOut = SeqParam;

if TEfill < seq.gradRasterTime
    TEfill = seq.gradRasterTime;
    display(['Your provided TE was too short, increased to ', num2str(minTE)]); % TE too short, increasing to the next gradient raster time
    SeqParamOut.TE = minTE; % update the TE in the output structure
end;
TEfill = ceil(TEfill / seq.gradRasterTime) * seq.gradRasterTime;
TRfill = TR/Nslices - minTR;
if TRfill < seq.gradRasterTime
    TRfill = seq.gradRasterTime;
    display(['% your provided TR was too short, increased to ', num2str(minTR * Nslices)]); % TR too short, increasing to the next gradient raster time
    SeqParamOut.TR = minTR * Nslices; % update the TR in the output structure
end;

SeqParamOut.TE = max(TE(1),minTE) +minDTE * [ (0:(ME-1))]


TRfill = ceil(TRfill / seq.gradRasterTime) * seq.gradRasterTime;

for s=1:Nslices
    seq.addBlock(rf_fs,gz_fs);
    rf.freqOffset=gz.amplitude*slicePositions(s);
    rf.phaseOffset=-2*pi*rf.freqOffset*mr.calcRfCenter(rf); % compensate for the slice-offset induced phase
    seq.addBlock(rf,gz,trig);
    seq.addBlock(mr.makeDelay(TEfill));
    seq.addBlock(gxPre,gyPre,gzReph);
    for echo=1:ME
    for i=1:Ny_meas
        if i==1
            seq.addBlock(gx,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
        elseif i==Ny_meas
            seq.addBlock(gx,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
        else
            seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
        end 
        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
    end
    if echo < ME
        seq.addBlock(gxPre_ME,gyPre_ME);
    % seq.addBlock(mr.makeDelay(minDTE)); this adds extra complexity to the sequence, so we will add it only later
    end
    end
    seq.addBlock(mr.makeDelay(TRfill));
end

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% do some visualizations
if Show.SequenceDiagram == 1;
    % seq.plot('stacked', 1);             % Plot all sequence waveforms in the new 'stacked' mode
    %seq.paperPlot('blockRange',[1 41]);
    seq.paperPlot('blockRange',[1 5+ME*Ny_meas])

end;
%seq.plot('timeDisp','us','showBlocks',1,'timeRange',[0 25e-3]); %detailed view

rf.freqOffset=0;
rf.phaseOffset=0;
[rf_bw,rf_f0,rf_spectrum,rf_w]=mr.calcRfBandwidth(rf);


%% trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();


if or(Show.KspaceTimeCourse == 1, Show.KspaceCoverage == 1);
    % plot k-spaces
    figure; 
end
if Show.KspaceTimeCourse == 1;
    nexttile;
    plot(t_ktraj, ktraj(1:2,:)'); % plot the entire k-space trajectory
    hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
    title('k-space vector components as functions of time');
    legend('k_x','k_y','ADC samples');
    xlabel('time, s');
    ylabel('k-space trajectory, 1/m');
end;

if Show.KspaceCoverage == 1;
    nexttile;
    plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
    hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
    axis('equal'); % enforce aspect ratio for the correct trajectory display
    axis ('tight')
    title('k-space trajectory (k_x/k_y)');
end;

if Show.SliceProfile == 1;
    figure;
    nexttile;
    plot(rf_w,abs(rf_spectrum)); hold on;
    title('Excitation pulse profile (low-angle approximation)');

  %  plot(rf_w,circshift(abs(rf_spectrum),0, round ( (thickness+sliceGap)/thickness * rf_bw / mean(diff(rf_w)) )),'g'); % shift the spectrum to the the next slice
  %  plot(rf_w,circshift(abs(rf_spectrum),0, -round ((thickness+sliceGap)/thickness * rf_bw / mean(diff(rf_w)) )),'r'); % shift the spectrum to the the next slice
    plot(rf_w -  ( (thickness+sliceGap)/thickness * rf_bw ),abs(rf_spectrum), 'g');
    plot(rf_w +  ( (thickness+sliceGap)/thickness * rf_bw ),abs(rf_spectrum), 'g');
    legend('slice -1','slice 0','slice 1');
    xlabel('Frequency, Hz');
    xlim(3*[-rf_bw rf_bw]);
    %keyboard
    nexttile;
    plot(t_slicepos, slicepos, '*');
    title('slice position (vector components) as a function or time');
    %axis off;


end;


%% prepare the sequence output for the scanner
seq.setDefinition('Name', 'epi'); 
seq.setDefinition('FOV', [fov fov max(slicePositions)-min(slicePositions)+thickness]);
seq.setDefinition('ReceiverGainHigh',1);
% the following definitions only have effect in conjunction with LABELs 
%seq.setDefinition('SlicePositions', slicePositions);
%seq.setDefinition('SliceThickness', thickness);
%seq.setDefinition('SliceGap', sliceGap);

seq.write('epi_rs.seq'); 

% seq.install('siemens');

if Show.Sound == 1;
    seq.sound(); % simulate the seq's tone
end;


