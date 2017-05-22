classdef TransceiverCharacterizationTests < matlab.unittest.TestCase
    % Test tranceivers (ad93XX) with test equipment
    %
    % -Harness will only configure SampleRate, CenterFrequency, and num 
    % samples to TX/RX automatically
    %   - SDRs will have RFBandwidth and FIRDataFilename configured as well
    %   - The gains will only adjust the gain of the SDR (TX attenuation)
    % -Any remaining configurations applied to the sink/source objects will
    % be maintained
    % 
    % When using an SDR an automatic calibration will be performed for
    % frequency compensation by tuning the XO_CORRECTION.  The
    % TX_ATTENUATION will be tuned as well when an SDR is acting as a
    % transmitter to prevent a PXA from saturating.
    %
    %
    % Example call:
    %  test = TransceiverCharacterizationTests;
    %  test.run()
    
    % User configured
    properties
        EVMTestsToAverage = 1 % Variance is normally small so leave this small
        CRCFailureTries = 5; % If we get LTE frame failures retry measurement
        EnableVisuals = false;
        FrequencyOffsetForceRecal = 500; % (Hz) Measured offset when a recalibation will be initiated
    end
    % Automatically configured
    properties (Access = private)
        RMConfig = 'R.7'
        FIRDataFilename = 'LTE10_MHz.ftr'
        SampleRate = 15.36e6
        RFBandwidth = 9e6
        XOOffset = 292;
        FreqOffset = 0;
        CellID = 17;
        % Initialization
        InstrumentsCleared = false;
        PlutoSetup = false;
        cfgACLR = []
    end
    
    methods(TestClassSetup)
        % Check for LST
        function checkForLicenses(testCase)
            % Check for LST presence
            testCase.verifyTrue(~isempty(ver('lte')),'LTE System Toolbox not installed/availabe')
        end
        % Clear any instruments connected
        function clearVSAInstruments(testCase)
            instr = instrfind;
            if ~isempty(instr)
                fclose(instr); 
            end
            testCase.InstrumentsCleared = true;
        end
        % Setup environment for SDR
        function setupPlutoEnv(testCase)
           d1 = sdrdev('Pluto'); setupSession(d1);
           testCase.PlutoSetup = true;
        end
        
    end
    
    methods(TestClassTeardown)
        function disconnectVSAInstruments(~)
            instr = instrfind;
            if ~isempty(instr)
                fclose(instr);
                instrreset;
            end
        end
    end
    
    methods (Static)
    end
    
    methods % Non-Static Test Scaffolding
        
        % Generate Power Plot
        function genPowerResultsPlot(testCase,sigpowerdB,gain,sigpowerdBVar,frequencies)
            %plot(frequencies,evm);
            freqs = repmat(frequencies.',1,size(sigpowerdB,2));
            if testCase.EVMTestsToAverage==1
                plot(freqs,sigpowerdB)
            else
                errorbar(freqs,sigpowerdB,sigpowerdBVar)
            end
            xlabel('Frequency (Hz)');
            ylabel('Inband Power (dBm)');
            title('Inband Power vs Frequency for LTE-10 QAM-64');
            ylim([min(min(sigpowerdB))-5 max(max(sigpowerdB))+5]);
            xlim([min(frequencies)-100 max(frequencies)+100]);
            str = {};
            for s = 1:length(gain(1,:))
                str = {str{:},['Gain: ',num2str(gain(1,s))]};  %#ok<CCAT>
            end
            legend(str);
            drawnow;
        end
        
        % Generate EVM Plot
        function genEVMResultsPlot(testCase,evm,gain,rssi,evmVar,rssiVar,frequencies)
            
            %plot(frequencies,evm);
            freqs = repmat(frequencies.',1,size(evm,2));
            if testCase.EVMTestsToAverage==1
                plot(freqs,evm)
            else
                errorbar(freqs,evm,evmVar)
            end
            xlabel('Frequency (Hz)');
            ylabel('EVM (dB)');
            title('EVM vs Frequency for LTE-10 QAM-64');
            ylim([min(min(evm))-5 max(max(evm))+5]);
            xlim([min(frequencies)-100 max(frequencies)+100]);
            str = {};
            for s = 1:length(gain(1,:))
                str = {str{:},['Gain: ',num2str(gain(1,s))]};  %#ok<CCAT>
            end
            legend(str);
            drawnow;
        end
        
        % Set class parameters based on LTE configuration set
        function LTEConfigurationLibrary(obj,LTEMode)
            if strcmp(LTEMode,'LTE1.4') == 1
                obj.RMConfig = 'R.4';
                obj.SampleRate = 1.92e6;
                obj.RFBandwidth = 1.08e6;
                obj.FIRDataFilename = 'LTE1p4_MHz.ftr';
            elseif strcmp(LTEMode,'LTE3') == 1
                obj.RMConfig = 'R.5';
                obj.SampleRate = 3.84e6;
                obj.RFBandwidth = 2.7e6;
                obj.FIRDataFilename = 'LTE3_MHz.ftr';
            elseif strcmp(LTEMode,'LTE5') == 1
                obj.RMConfig = 'R.6';
                obj.SampleRate = 7.68e6;
                obj.RFBandwidth = 4.5e6;
                obj.FIRDataFilename = 'LTE5_MHz.ftr';
            elseif strcmp(LTEMode,'LTE10') == 1
                obj.RMConfig = 'R.7';
                obj.SampleRate = 15.36e6;
                obj.RFBandwidth = 9e6;
                obj.FIRDataFilename = 'LTE10_MHz.ftr';
            else
                error('Please input LTE1.4, LTE3, LTE5 or LTE10.');
            end 
        end
        
        % Measure ACLR (ACPR)
        function aclr = ACLRMeasurement(testCase,RxIQ)
           
            % Measure ACLR for both UTRA and E-UTRA scenarios
            aclr = hACLRMeasurementEUTRA(testCase.cfgACLR,RxIQ);
            aclr = hACLRMeasurementUTRA(aclr,RxIQ,testCase.cfgACLR.nRC,testCase.cfgACLR.Rc,testCase.cfgACLR.UTRAbw);
            
        end
        
        % Measure in-band power
        function EUTRAdB = MeasureInBandPower(testCase,resampled)
            
            Ndft = 12288;
            resampled = resampled(1:Ndft);
            %window(@blackmanharris,Ndft)
            win = repmat(window(@blackmanharris,Ndft),1,1);
            win = win * sqrt(size(win,1)/sum(win(:,1).^2));
            demod = fftshift(fft(resampled.*win,Ndft),1)/Ndft;
            
            % Calculate energy in relevant DFT bins
            m = abs(demod.^2);
            
            % Calculate power
            p = sum(m(:))/size(m,2);
            EUTRAdB = 10*log10(p)+30;
            
%             %%
%             aclr = testCase.cfgACLR.aclr;
%             
%             gcdBWConfig = gcd(aclr.SamplingRate,aclr.BandwidthConfig);
%             Ndft = aclr.SamplingRate/gcdBWConfig;
%             
%             % Also ensure that Ndft is even
%             multEvenNdft = 1+mod(Ndft,2);
%             Ndft = Ndft*multEvenNdft;
%             
%             % Calculate NbinsConfig, the number of DFT bins which spans channel
%             % transmission bandwidth configuration of interest.
%             NbinsConfig = aclr.BandwidthConfig*multEvenNdft/gcdBWConfig;
%             multEvenNbinsConfig = 1+mod(NbinsConfig,2);
%             
%             % Also ensure that NbinsConfig is even
%             NbinsConfig = NbinsConfig*multEvenNbinsConfig;
%             Ndft = Ndft*multEvenNbinsConfig;
%             
%             % Increase the DFT size, if necessary, so that at least 15kHz frequency
%             % resolution is obtained.
%             delta_f = aclr.SamplingRate/Ndft;
%             multFreqRes15kHz = ceil(delta_f/15e3);
%             Ndft = Ndft*multFreqRes15kHz;
%             NbinsConfig = NbinsConfig*multFreqRes15kHz;
%             
%             % Calculate NbinsChannel, the number of DFT bins between adjacent
%             % E-UTRA channel center frequencies, increasing the number of DFT bins
%             % if necessary to make NbinsChannel an integer.
%             gcdBW = gcd(aclr.SamplingRate,aclr.Bandwidth);
%             NdftTemp = aclr.SamplingRate/gcdBW;
%             multNbinsChannelInt = NdftTemp/gcd(Ndft,NdftTemp);
%             Ndft = Ndft*multNbinsChannelInt;
%             
%             NbinsConfig = NbinsConfig*multNbinsChannelInt;
%             NbinsChannel = aclr.Bandwidth*Ndft/aclr.SamplingRate;
%             
%             % OFDM demodulate the measurement signal using the DFT.
%             Ndfts = floor(size(resampled,1)/Ndft);
%             dftinput = reshape(resampled(1:(Ndfts*Ndft)),Ndft,Ndfts);
%             win = repmat(window(@blackmanharris,Ndft),1,Ndfts);
%             win = win * sqrt(size(win,1)/sum(win(:,1).^2));
%             demod = fftshift(fft(dftinput.*win,Ndft),1)/Ndft;
%             
%             % For each channel index in the array [-2 -1 0 1 2] record the channel
%             % center frequency and calculate the channel power by calculating the
%             % energy in relevant DFT bins.
%             EUTRAdB = zeros(1,5);
%             for i = -2:2
%                 
%                 % Calculate energy in relevant DFT bins
%                 m = abs(demod((Ndft/2)+1+(i*NbinsChannel)+...
%                     (-(NbinsConfig/2):NbinsConfig/2),:).^2);
%                 
%                 % Calculate power
%                 p = sum(m(:))/size(m,2);
%                 EUTRAdB(i+3) = 10*log10(p)+30;
%                 
%                 % Record center frequency
%                 aclr.EUTRACenterFreq(i+3) = i*NbinsChannel/Ndft*aclr.SamplingRate;
%                 
%             end
% 
%             ibPowerDBM = EUTRAdB(3);
        end
        
        % Recover LTE signal
        function [finalEVM_db,CRCError,doRecal] = LTEReceiver(testCase,Rx)
            
            rxWaveform = reshape(Rx,[length(Rx),1]);
            
            % Derived parameters
            samplesPerFrame = 10e-3*testCase.SampleRate; % LTE frames period is 10 ms
            
            % *LTE Setup*
            % The parameters for decoding the MIB are contained in the structure |enb|.
            % FDD duxplexing mode and a normal cyclic prefix length are assumed. Four
            % cell-specific reference ports (CellRefP) are assumed for the MIB decode.
            % The number of actual CellRefP is provided by the MIB.
            
            enb.DuplexMode = 'FDD';
            enb.CyclicPrefix = 'Normal';
            enb.CellRefP = 4;
            
            % The sampling rate of the signal controls the captured bandwidth. The
            % number of number of RBs captured is obtained from a lookup table using
            % the chosen sampling rate.
            % Bandwidth: {1.4 MHz, 3 MHz, 5 MHz, 10 MHz}
            SampleRateLUT = [1.92e6 3.84e6 7.68e6 15.36e6];
            NDLRBLUT = [6 15 25 50];
            enb.NDLRB = NDLRBLUT(SampleRateLUT==testCase.SampleRate);
            
            % Channel estimation configuration using cell-specific reference signals. A
            % conservative 9-by-9 averaging window is used to minimize the effect of
            % noise.
            cec.FreqWindow = 9;               % Frequency averaging window in Resource Elements (REs)
            cec.TimeWindow = 9;               % Time averaging window in REs
            cec.InterpType = 'Cubic';         % Cubic interpolation
            cec.PilotAverage = 'UserDefined'; % Pilot averaging method
            cec.InterpWindow = 'Centred';     % Interpolation windowing method
            cec.InterpWinSize = 3;            % Interpolate up to 3 subframes simultaneously
            
            % *Signal Processing*
            %
            % For each
            % captured frame the MIB is decoded and if successful the CFI and the PDCCH
            % for each subframe are decoded and channel estimate and equalized PDCCH
            % symbols are shown.
            enbDefault = enb;
            
            % Set default LTE parameters
            enb = enbDefault;
            
            % Perform frequency offset correction
            frequencyOffset = lteFrequencyOffset(enb,rxWaveform);
            doRecal = abs(frequencyOffset)>testCase.FrequencyOffsetForceRecal;
            rxWaveformFreqCorrected = lteFrequencyCorrect(enb,rxWaveform,frequencyOffset);
            fprintf('Corrected a frequency offset of %g Hz.\n',frequencyOffset)
            
            % Perform the blind cell search to obtain cell identity
            cellSearch.SSSDetection = 'PostFFT'; cellSearch.MaxCellCount = 1;
            [NCellID,frameOffset] = lteCellSearch(enb,rxWaveform,cellSearch);
            if testCase.CellID~=NCellID
                finalEVM_db = NaN;
                CRCError = true;
                fprintf('Bad cell identity of %i.\n',NCellID);
                return
            else
                fprintf('Detected a cell identity of %i.\n',NCellID);
            end
            
            % Sync the captured samples to the start of an LTE frame, and trim off
            % any samples that are part of an incomplete frame.
            fprintf('frameoffset %i\n', frameOffset);
            rxWaveform = rxWaveformFreqCorrected(frameOffset+1:end);
            tailSamples = mod(length(rxWaveform),samplesPerFrame);
            rxWaveform = rxWaveform(1:end-tailSamples);
            
            % EVM Calculation
            rmc = lteRMCDL(testCase.RMConfig);
            rmc.NCellID = 17;
            rmc.NFrame = 700;
            rmc.TotSubframes = 8*10; % 10 subframes per frame
            %rmc.OCNG = 'On'; % Add noise to unallocated PDSCH resource elements
            rmc.OCNGPDSCHEnable = 'On';
            rmc.OCNGPDCCHEnable = 'On';
            %rmc.PDSCH.Modulation = '16QAM'; % default is 64QAM
            
            % Generate RMC waveform
            trData = [1;0;0;1]; % Transport data
            [~,~,rmc] = lteRMCDLTool(rmc,trData);
            [EVM,CRCError] = PDSCHEVM(rmc,cec,rxWaveform);
            finalEVM_db= 10*log10((EVM.RMS)^2);
            fprintf('Averaged overall EVM: %0.3f dB\n', finalEVM_db);
            
        end
        
        % Generate LTE waveform with specific configuration
        function [eNodeBOutput,SampleRate,cellID] = generateLTETestWaveform(testCase)
            % Generate the LTE signal
            txsim.RC = testCase.RMConfig; % Base RMC configuration
            txsim.NCellID = 17;       % Cell identity
            txsim.NFrame = 700;       % Initial frame number
            txsim.TotFrames = 1;      % Number of frames to generate
            txsim.RunTime = 20;       % Time period to loop waveform in seconds
            txsim.DesiredCenterFrequency = 1e9; % Center frequency in Hz (NOT USED)
            cellID = txsim.NCellID;
            
            % Generate RMC configuration and customize parameters
            rmc = lteRMCDL(txsim.RC);
            rmc.NCellID = txsim.NCellID;
            rmc.NFrame = txsim.NFrame;
            rmc.TotSubframes = txsim.TotFrames*10; % 10 subframes per frame
            %rmc.OCNG = 'On'; % Add noise to unallocated PDSCH resource elements
            rmc.OCNGPDSCHEnable = 'On';
            rmc.OCNGPDCCHEnable = 'On';
            %rmc.PDSCH.Modulation = '16QAM'; % default is 64QAM
            
            % Generate RMC waveform
            trData = [1;0;0;1]; % Transport data
            [eNodeBOutput,~,rmc] = lteRMCDLTool(rmc,trData);
            %txsim.SamplingRate = rmc.SamplingRate;
            SampleRate = rmc.SamplingRate;

        end
        
        % Run ACLR tests over gain and frequency
        function [sigpowerdB, sigpowerdBVar] = testACLROverGainAndFrequency(testCase, SDRgains, frequencies, source, sink)
            
            [sigpowerdB, sigpowerdBVar, gain] = deal(NaN(length(frequencies),length(SDRgains)));
            
            % Define filename for log file
            d = dbstack(); callname = d(2).name; callname = strsplit(callname,'.');
            filename = [callname{end},' ',datestr(now,'dd-mmm-yyyy HH_MM_SS'),'.mat'];
            
            % Generate Waveform
            txWaveform = testCase.generateLTETestWaveform();
            
            % ACLR Configuration
            % Get necessary sampling rates for LTE config
            rmc = lteRMCDL(testCase.RMConfig);
            info = lteOFDMInfo(rmc);
            rmc.SamplingRate = info.SamplingRate;
            % Get necessary sampling rate for ACLR measurement
            [aclr, nRC, Rc, UTRAbw] = hACLRParameters(rmc);
            testCase.cfgACLR = struct;
            testCase.cfgACLR.aclr = aclr;
            testCase.cfgACLR.nRC = nRC;
            testCase.cfgACLR.Rc = Rc;
            testCase.cfgACLR.UTRAbw = UTRAbw;
            testCase.cfgACLR.ACLRSampleRate = aclr.OSR*rmc.SamplingRate;
            
            % Cycle over gains
            for gainIndx = 1:length(SDRgains)
                fprintf(['-------------------------\n',...
                    'Testing at gain: %f (dB) [%d of %d]\n'],...
                    SDRgains(gainIndx),gainIndx,length(SDRgains));
                
                % Set devices for LTE configuration and gain
                testCase.Configure(source,sink,SDRgains(gainIndx));
                % CalibrateXO and set RX/TX gain so we are not saturating
                testCase.Calibrate(source,sink,max(frequencies));
                
                % Set receiver to larger bandwidth for ACPR calculation
                %sink.SampleRate = testCase.cfgACLR.ACLRSampleRate;
                
                % Setup actual devices with TX/RX in frequency independent config
                source.Setup(length(txWaveform));
                sink.Setup(length(txWaveform));
                
                % Cycle across frequencies
                for freqIndx = 1:length(frequencies)
                    fprintf(['-------------------------\n',...
                        'Testing at frequency: %f (Hz) [%d of %d]\n'],...
                        frequencies(freqIndx),freqIndx,length(frequencies));
                    
                    % Change transmitter center frequency
                    source.CenterFrequency = frequencies(freqIndx);
                    sink.CenterFrequency = frequencies(freqIndx);
                    
                    % Transmit signal again with new settings
                    source.EnableTX(txWaveform);
                    
                    % Receive and evaluate EVM
                    [sigpowerdB(freqIndx,gainIndx),...
                        sigpowerdBVar(freqIndx,gainIndx)] =...
                        testACLRForConfig(testCase, sink);
                    gain(freqIndx,gainIndx) = SDRgains(gainIndx);
                    
                    % Save (overwrite) new data to file
                    testCase.SaveDataPoints(filename,...
                        'freqIndx',freqIndx,'gainIndx',gainIndx,...
                        'sigpowerdB',sigpowerdB,'sigpowerdBVar',sigpowerdBVar,...
                        'gain',gain,'frequencies',frequencies);
                    % Plot
                    testCase.genPowerResultsPlot(sigpowerdB,gain,sigpowerdBVar,frequencies);
                    
%                     % Force recalibrate
%                     if doRecal
%                         % Clean up devices
%                         source.Release();
%                         sink.Release();
%                         % Set devices for LTE configuration and gain
%                         testCase.Configure(source,sink,SDRgains(gainIndx));
%                         % CalibrateXO and set RX/TX gain so we are not saturating
%                         testCase.Calibrate(source,sink,max(frequencies));
%                         % Setup actual devices with TX/RX in frequency independent config
%                         source.Setup(length(txWaveform));
%                         sink.Setup(length(txWaveform)*5);
%                     end
                end
                % Clean up devices
                source.Release();
                sink.Release();
            end
            % Done
            fprintf('Log data written to: %s\n',filename);
        end

        % Run individual ACLR test
        function [sigpowerdB,sigpowerdBVar] = testACLRForConfig(testCase, sink)
            
            [sigpowerdB] = deal(NaN(testCase.EVMTestsToAverage,1));
            % Enable visuals
            if testCase.EnableVisuals
                SpecAnalyzer = dsp.SpectrumAnalyzer(...
                    'SampleRate',testCase.SampleRate);
            end
            % Run receiver
            for m = 1:testCase.EVMTestsToAverage
                fprintf('-- EVM Test [%d of %d]\n',m,testCase.EVMTestsToAverage);
                % Receive Data
                RxIQ = sink.EnableRX();
                if testCase.EnableVisuals
                    SpecAnalyzer(RxIQ);
                    SpecAnalyzer.release(); % Force scope to autoscale
                    pause(1);
                end
                % Measure inband power
                sigpowerdB(m) = testCase.MeasureInBandPower(RxIQ);
            end
            % Get variance
            sigpowerdBVar = var(sigpowerdB);
            % Get mean
            sigpowerdB = mean(sigpowerdB);
        end
        
        % Run individual EVM test
        function [evm, rssi, gain, evmVar, rssiVar, gainVar, doRecal] = testEVMForConfig(testCase, sink)
            
            [evm, rssi, gain] = deal(NaN(testCase.EVMTestsToAverage,1));
            
            % Enable visuals
            if testCase.EnableVisuals
                SpecAnalyzer = dsp.SpectrumAnalyzer(...
                    'SampleRate',testCase.SampleRate);
            end
            % Run receiver
            for m = 1:testCase.EVMTestsToAverage
                fprintf('-- EVM Test [%d of %d]\n',m,testCase.EVMTestsToAverage);
                % Receive Data
                Tries = testCase.CRCFailureTries;
                while Tries
                    RxIQ = sink.EnableRX();
                    if testCase.EnableVisuals
                        SpecAnalyzer(RxIQ);
                        SpecAnalyzer.release(); % Force scope to autoscale
                        pause(1);
                    end
                    % Measure EVM in dB
                    try % Receiver can fail in rare occassions
                        [evm(m),CRCError,doRecal] = testCase.LTEReceiver(RxIQ);
                    catch
                       evm(m) = NaN;
                       CRCError = true;
                    end
                    if ~CRCError
                        Tries = 0; % done
                    else
                        Tries = Tries - 1; % try again
                    end
                end
            end
            % Get variance
            evmVar = var(evm);
            rssiVar = var(rssi);
            gainVar = var(gain);
            % Get mean
            evm = mean(evm);
            rssi = mean(rssi);
            gain = mean(gain);
        end
                
        % Run EVM tests over gain and frequency
        function [evm, rssi, gain, evmVar, rssiVar] = testEVMOverGainAndFrequency(testCase, SDRgains, frequencies, source, sink)
            
            [evm, rssi, gain, evmVar, rssiVar] = deal(NaN(length(frequencies),length(SDRgains)));
            
            % Define filename for log file
            d = dbstack(); callname = d(2).name; callname = strsplit(callname,'.');
            filename = [callname{end},' ',datestr(now,'dd-mmm-yyyy HH_MM_SS'),'.mat'];
            
            % Generate Waveform
            txWaveform = testCase.generateLTETestWaveform();
            
            % Cycle over gains
            for gainIndx = 1:length(SDRgains)
                fprintf(['-------------------------\n',...
                    'Testing at gain: %f (dB) [%d of %d]\n'],...
                    SDRgains(gainIndx),gainIndx,length(SDRgains));
                
                % Set devices for LTE configuration and gain
                testCase.Configure(source,sink,SDRgains(gainIndx));
                % CalibrateXO and set RX/TX gain so we are not saturating
                testCase.Calibrate(source,sink,max(frequencies));
                % Setup actual devices with TX/RX in frequency independent config
                source.Setup(length(txWaveform));
                sink.Setup(length(txWaveform)*5);
                
                % Cycle across frequencies
                for freqIndx = 1:length(frequencies)
                    fprintf(['-------------------------\n',...
                        'Testing at frequency: %f (Hz) [%d of %d]\n'],...
                        frequencies(freqIndx),freqIndx,length(frequencies));
                    
                    % Change transmitter center frequency
                    source.CenterFrequency = frequencies(freqIndx);
                    sink.CenterFrequency = frequencies(freqIndx);
                    
                    % Transmit signal again with new settings
                    source.EnableTX(txWaveform);
                    
                    % Receive and evaluate EVM
                    [evm(freqIndx,gainIndx), rssi(freqIndx,gainIndx), ~,...
                        evmVar(freqIndx,gainIndx), rssiVar(freqIndx,gainIndx),...
                        ~,doRecal] =...
                        testEVMForConfig(testCase, sink);
                    gain(freqIndx,gainIndx) = SDRgains(gainIndx);
                    
                    % Save (overwrite) new data to file
                    testCase.SaveDataPoints(filename,...
                        'freqIndx',freqIndx,'gainIndx',gainIndx,...
                        'evm',evm,'evmVar',evmVar,...
                        'gain',gain,'frequencies',frequencies);
                    % Plot
                    testCase.genEVMResultsPlot(evm,gain,rssi,evmVar,rssiVar, frequencies);
                    
                    % Force recalibrate
                    if doRecal
                        % Clean up devices
                        source.Release();
                        sink.Release();
                        % Set devices for LTE configuration and gain
                        testCase.Configure(source,sink,SDRgains(gainIndx));
                        % CalibrateXO and set RX/TX gain so we are not saturating
                        testCase.Calibrate(source,sink,max(frequencies));
                        % Setup actual devices with TX/RX in frequency independent config
                        source.Setup(length(txWaveform));
                        sink.Setup(length(txWaveform)*5);
                    end
                end
                % Clean up devices
                source.Release();
                sink.Release();
            end
            % Done
            fprintf('Log data written to: %s\n',filename);
        end

        % Calibrate XO and TX/RX Gains
        function Calibrate(~,source,sink,CenterFrequency)
            IsSDRSource = isequal(class(source),'iioSDR');
            IsSDRSink = isequal(class(sink),'iioSDR');
            % Calibrate if we have an SDR and we are not transmitting to
            % the same SDR
            if IsSDRSource && ~IsSDRSink
                source.CalibrateRXTXGain(sink, CenterFrequency);
                source.CalibrateXO(sink, CenterFrequency);
            elseif ~IsSDRSource && IsSDRSink
                sink.CalibrateRXTXGain(sink, CenterFrequency);
                sink.CalibrateXO(source, CenterFrequency);
            elseif (IsSDRSource && IsSDRSink) && (source.PI ~= sink.IP)
                source.CalibrateXO(sink, CenterFrequency);
            end
        end
        
        % Configure Sources/Sinks and Calibrate SDR
        function Configure(testCase,source,sink,varargin)
            % Make we have a radio case
            IsSDRSource = isequal(class(source),'iioSDR');
            IsSDRSink = isequal(class(sink),'iioSDR');
            % Set static configurations
            source.SampleRate = testCase.SampleRate;
            sink.SampleRate = testCase.SampleRate;
            if IsSDRSource || IsSDRSink
                % Set static configurations
                if IsSDRSource % TX
                    source.RFBandwidth = testCase.RFBandwidth;
                    source.FIRDataFilename = testCase.FIRDataFilename;
                    if ~isempty(varargin)
                        source.TXAttenuation = varargin{1};
                        %source.RXGain = varargin{2};
                    end
                end % RX
                if IsSDRSink
                   sink.RFBandwidth = testCase.RFBandwidth;
                   sink.FIRDataFilename = testCase.FIRDataFilename;
                   if ~isempty(varargin)
                       %sink.TXAttenuation = varargin{1};
                       sink.RXGain = varargin{1};
                   end
                end
            end
               
        end
        
        % Save logs
        function SaveDataPoints(~,filename,varargin)
            % Build struct
            variablesToSave = cell2struct({varargin{2:2:end}},{varargin{1:2:end}},2); %#ok<NASGU>
            save(filename,'-struct', 'variablesToSave');
        end
    end
    
    
    methods (Test)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Tests
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % In development
        function testEVMR7_SDRtoSDR(testCase)
            tic
            % Set up configuration based on LTE Configuration
            testCase.LTEConfigurationLibrary('LTE10');
            % Set up source
            source = iioSDR('IP','192.168.2.1','Mode','TX','dev_name','pluto');
            % Set up sink
            sink = iioSDR('IP','192.168.2.1','Mode','RX','dev_name','pluto');
            % Calculate frequencies to evaluate
            startFreq = 1e9;
            freqStep = 1e9;
            endFreq = 5e9;
            frequencies = startFreq:freqStep:endFreq;
            % Enable Visuals
            testCase.EnableVisuals = true;
            % Set averaging
            testCase.EVMTestsToAverage = 1;
            % Run
            [evm, rssi, gain, evmVar, rssiVar] = testCase.testEVMOverFrequency(frequencies, source, sink);
            % Create Plot
            testCase.genEVMResultsPlot(evm,gain,rssi,evmVar,rssiVar, frequencies);
            toc
        end

        % Well tested
        function testEVMR7_SDRtoPXA(testCase)
            tic
            % Set up parameters based on LTE Configuration
            testCase.LTEConfigurationLibrary('LTE10');
            % Set up source
            %source = iioSDR('IP','192.168.2.1','Mode','TX','dev_name','pluto');
            source = iioSDR('IP','10.66.98.185','Mode','TX','dev_name','ad9361');
            % Set up sink
            sink = N9030A('IP','10.66.100.143');
            % Calculate frequencies to evaluate
            startFreq = 70e6+1e6;
            freqStep = 1e7;
            endFreq = 6e9-1e6;
            frequencies = startFreq:freqStep:endFreq;
            % Pick gains to consider
            TXgains = [-10 -10 -10 -10];
            % Enable Visuals
            %testCase.EnableVisuals = true;
            % Run            
            [evm, rssi, gain, evmVar, rssiVar] = ...
                testCase.testEVMOverGainAndFrequency(TXgains, frequencies, source, sink);
            % Create Plot
            testCase.genEVMResultsPlot...
                (evm,gain,rssi,evmVar,rssiVar, frequencies);
            toc
        end
        
        % In development
        function testACLRR7_SDRtoPXA(testCase)
            tic
            % Set up parameters based on LTE Configuration
            testCase.LTEConfigurationLibrary('LTE10');
            % Set up source
            source = iioSDR('IP','192.168.2.1','Mode','TX','dev_name','pluto');
            % Set up sink
            sink = N9030A('IP','10.66.100.143');
            % Calculate frequencies to evaluate
            startFreq = 70e6+1e6;
            freqStep = 1e7;
            endFreq = 6e9-1e6;
            frequencies = startFreq:freqStep:endFreq;
            % Pick gains to consider
            TXgains = [-10 -20 -30 -40];
            % Enable Visuals
            testCase.EnableVisuals = true;
            % Run
            testCase.testACLROverGainAndFrequency(TXgains, frequencies, source, sink);
            toc
        end
        
        % Stable
        function testEVMR7_SigGenToSDR(testCase)
            tic
            % Set up source
            source = N5182B('IP','10.66.98.158','OutputPower',-30);
            % Set up sink
            sink = iioSDR('IP','10.66.98.185',...
                          'Mode','RX',...
                          'dev_name','ad9361',...
                          'AGCMode','slow_attack');
            % Calculate frequencies to evaluate
            startFreq = 70e6;
            freqStep = 100e6;
            endFreq = 5e9;
            frequencies = startFreq:freqStep:endFreq;
            % Increase EVM test to average
            testCase.EVMTestsToAverage = 5;
            % Enable Visuals
            testCase.EnableVisuals = true;
            % SDR RX gains (ONLY USED IF AGCMode is manual)
            RXgains = [25];
            % Run
            [evm, rssi, gain, evmVar, rssiVar] = ...
                testCase.testEVMOverGainAndFrequency(RXgains, frequencies, source, sink);
            % Create Plot
            testCase.genEVMResultsPlot(evm,gain,rssi,evmVar,rssiVar, frequencies);
            toc
        end
        
        % Well tested
        function testEVMR7_SigGenToPXA(testCase)
            tic
            % NOTE: to correct any measure offset updated the
            % 'FrequencyBias' parameter of the sink object
            
            % Set up source
            source = N5182B('IP','10.66.98.158','OutputPower',-30);
            % Set up sink
            sink = N9030A('IP','10.66.100.143','Attenuation',0,'FrequencyBias',115);
            % Calculate frequencies to evaluate
            startFreq = 1e9;
            freqStep = 1e9;
            endFreq = 5e9;
            frequencies = startFreq:freqStep:endFreq;
            % Enable Visuals
            testCase.EnableVisuals = true;
            % Pick gains to consider
            % THESE MEAN NOTHING WHEN USING PXA AND MXG TOGETHER
            SDRgains = [0 0 0 0];
            % Enable Visuals
            %testCase.EnableVisuals = true;
            % Run            
            [evm, rssi, gain, evmVar, rssiVar] = ...
                testCase.testEVMOverGainAndFrequency(SDRgains, frequencies, source, sink);
            % Save data to file
            testCase.SaveDataPoints({'evm', 'frequencies','evmVar'});
            % Create Plot
            testCase.genEVMResultsPlot(evm,gain,rssi,evmVar,rssiVar, frequencies);
            toc
        end
        
    end
    
    
    
end