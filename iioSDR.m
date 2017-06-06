classdef iioSDR < matlab.System
    %iioSDR Software-Defined Radio interface
    
    properties
        Mode = 'TX'; % | TX or RX
        inputs = 0; % 2 For 1 complex input | 4 for 2 complex inputs
        outputs = 0; % 2 For 1 complex output | 4 for 2 complex outputs
        CenterFrequency = 1e9; % Hz
        SampleRate = 15.36e6; % Hz
        RFBandwidth = 9e6; % Hz
        IP = '192.168.2.1';
        txLoOffset = 0; % Hz
        FrameLength = 2^17; % Samples
        %XOOffset = 0; % Hz for LO
        dev_name = 'pluto'; % Config filename without .cfg extension
        FIRDataFilename = 'LTE10_MHz.ftr'; % Filter filename
        AGCMode = 'manual';
        RXGain = 0; % dB
        TXAttenuation = -30; % dB
        ScaleInput = true; % Scale to best use int16 range
        RXLoops = 10;
        Channel = 1;
        XOCorrectionFrequencyTolerance = 100; % Hz
        XOCorrectionTries = 15; % Tries (only really needed on RX tests)
    end
    
    properties (Access = private)
       radioDev % Handle to radio object
       inputStruct
       XOCorrectionValue = 0;
       DCXOCorrectionValue = 0;
       RadioSetup = false;
    end
    
    methods
        % Constructor
        function obj = iioSDR(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:})
%             % Add listeners
%             obj.PropListener();
            
        end
        
        function Setup(obj,varargin)
            if nargin>1
                obj.FrameLength = varargin{1};
            end
            %obj.CenterFrequency = CenterFrequency;
            % Connect to radio
            buildADIRadioSystemObject(obj);
            obj.RadioSetup = true;
        end
        
        function buildADIRadioSystemObject(obj)
            % Construct radio object
            obj.radioDev = obj.setupRadio();
            % Build example input struct
            obj.SetupInputStruct();
        end
        
        function SetupInputStruct(obj)
            % Build input struct
            input = cell(1, obj.radioDev.in_ch_no + length(obj.radioDev.iio_dev_cfg.cfg_ch));
            % Set the attributes of AD936X
            if obj.inputs>0
                input{1} = zeros(obj.FrameLength,1);
                input{2} = zeros(obj.FrameLength,1);
            end
            input{obj.radioDev.getInChannel('RX_LO_FREQ')} = obj.CenterFrequency;
            input{obj.radioDev.getInChannel('RX_SAMPLING_FREQ')} = obj.SampleRate;
            input{obj.radioDev.getInChannel('RX_RF_BANDWIDTH')} = obj.RFBandwidth;
            if strcmp(obj.dev_name,'pluto')
            input{obj.radioDev.getInChannel('RX_GAIN_MODE')} = obj.AGCMode;
            input{obj.radioDev.getInChannel('RX_GAIN')} = obj.RXGain;
            else
            input{obj.radioDev.getInChannel('RX1_GAIN_MODE')} = obj.AGCMode;
            input{obj.radioDev.getInChannel('RX1_GAIN')} = obj.RXGain;
            %input{obj.radioDev.getInChannel('RX2_GAIN_MODE')} = obj.AGCMode;
            %input{obj.radioDev.getInChannel('RX2_GAIN')} = obj.RXGain;
            end
            input{obj.radioDev.getInChannel('TX_LO_FREQ')} = obj.CenterFrequency + obj.txLoOffset;
            input{obj.radioDev.getInChannel('TX_SAMPLING_FREQ')} = obj.SampleRate;
            input{obj.radioDev.getInChannel('TX_RF_BANDWIDTH')} = obj.RFBandwidth;
            input{obj.radioDev.getInChannel('TX_RF_ATTENUATION')} = obj.TXAttenuation;
            if strcmp(obj.dev_name,'fmcomms3')
                input{obj.radioDev.getInChannel('DCXO_FINE_TUNE')} = 5765+obj.DCXOCorrectionValue;            
            else
                input{obj.radioDev.getInChannel('XO_CORRECTION')} = 40e6+obj.XOCorrectionValue;
            end
            % Save to class
            obj.inputStruct = input;
        end
        
        function UpdateInputStructField(obj,fieldname,value)
            obj.inputStruct{obj.radioDev.getInChannel(fieldname)} = value;
        end
        
        function XOCorrectionValue = getXOCorrectionValue(obj)
            XOCorrectionValue = obj.XOCorrectionValue;
        end
        
        % Setup radio object
        function radio = setupRadio(obj)
            % System Object Configuration
            r = iio_sys_obj_matlab; % MATLAB libiio Constructor
            r.ip_address = obj.IP;
            r.dev_name = obj.dev_name; % config file name (without .cfg)
            switch obj.Mode
                case 'TX'
                    i=2;o=0;
                case 'RX'
                    i=0;o=2;
                otherwise
                    i=2;o=2;
            end
            obj.inputs = i; obj.outputs = o;
            r.in_ch_no = i;
            r.out_ch_no = o;
            r.in_ch_size = obj.FrameLength*(i/2);
            r.out_ch_size = obj.FrameLength*(o/2);
            % Go out to device
            radio = r.setupImpl();
            % Configure FIR
            if ~isempty(radio)
                radio.writeFirData(obj.FIRDataFilename);
            end
        end
        
        % Generic abstract call to start transmitting data
        function EnableTX(obj,txWaveform)
            Transmit(obj,txWaveform);
        end
        % Generic abstract call to start receiving data
        function RxIQ = EnableRX(obj)
            %obj.CenterFrequency = CenterFrequency;
            %SetupInputStruct(obj)
            
            % Update struct
            obj.SetupInputStruct();
            RxIQ = obj.Receive();
        end
        
        % Transmit
        function Transmit(obj,txWaveform)
            % Check input vector size
            if length(txWaveform)~=obj.radioDev.in_ch_size
                error('Input vector incorrect size');
            elseif obj.radioDev.in_ch_no < 1
                error('No input channels configured');
            else
                if obj.ScaleInput
                    txWaveform = obj.ScaleInputData(txWaveform);
                end
                % Update struct
                obj.SetupInputStruct();
                % Send out radio
                obj.inputStruct{1} = real(txWaveform);
                obj.inputStruct{2} = imag(txWaveform);
                stepImpl(obj.radioDev, obj.inputStruct);
            end
        end
        
        % Receive
        function [RxIQ, rssi, gain] = Receive(obj)
            % Check ports
            if obj.radioDev.out_ch_no < 1
                error('No input channels configured');
            else
                for d=1:obj.RXLoops
                    d = stepImpl(obj.radioDev, obj.inputStruct); %#ok<FXSET>
                end
                % Extract signal
                if obj.Channel == 1
                    RxIQ = complex(d{1},d{2});
                elseif channel == 2
                    RxIQ = complex(d{3},d{4});
                else
                    error('Please select a channel: 1 or 2.');
                end
                % Get measurements
                rssi = -1;
                gain = -1;
                % Remove cruff
                RxIQ = RxIQ(round(length(RxIQ)/4):end);
            end
        end
        
        % Perform gain calibration for PXA
        function CalibrateRXTXGain(obj, OtherDevice, CenterFrequency)
            % Generate Tone
            Fs = obj.SampleRate;
            Fc = 1e6; % Frequency of tone
            fftSize = 2^20;
            Samples = fftSize;
            t = 1/Fs:1/Fs:Samples/Fs;
            txWaveform = complex(sin(2*pi*Fc*t+pi/2)*1024,sin(2*pi*Fc*t)*1024);
            
            %sa = dsp.SpectrumAnalyzer;
            %sa.SampleRate = Fs;
            
            % Determine device set up
            SDRisTX = strcmp(obj.Mode,'TX');%(obj.inputs~=0) && (obj.outputs==0);
            
            orig_state = warning;
            
            % Set up transmitter and receiver
            if SDRisTX
                % Set up this object as transmitter
                obj.CenterFrequency = CenterFrequency;
                obj.Setup(length(txWaveform));
                obj.EnableTX(txWaveform);
                % Set up other object as receiver
                OtherDevice.CenterFrequency = CenterFrequency;
                OtherDevice.Setup(length(txWaveform));
                sortedAttenuations = sort(OtherDevice.PossibleAttenuations);
                % Disable warnings
                %warning('off','all')
                for k=1:length(OtherDevice.PossibleAttenuations)
                    % Update attenuation
                    OtherDevice.Attenuation = sortedAttenuations(k);
                    %OtherDevice.EnableRX();
                    [~,errors] = OtherDevice.EnableRX();
                    pause(1); % Wait for measurement to overload
                    aol = OtherDevice.checkADCOverloadCondition();
                    if isempty(errors) && ~aol
                        fprintf('Best RX attenuation %d\n',sortedAttenuations(k));
                        % Cleanup
                        OtherDevice.Release();
                        obj.Release();
                        warning(orig_state);
                        return
                    end
                end
            else
                % NO CALIBRATION IS DONE HERE, TX (MXG) REMAINS FIXED
                %error('SDR RX gain calibration not set up');
                
                % Set up other object as transmitter
                %OtherDevice.CenterFrequency = CenterFrequency;
                %OtherDevice.Setup(length(txWaveform));
                %OtherDevice.EnableTX(txWaveform);
                % Set up this object as receiver
                %obj.CenterFrequency = CenterFrequency;
                %obj.Setup(length(txWaveform)*4); 
                return
            end
            
            OtherDevice.Release(); % Cleanup
            obj.Release();
            warning(orig_state);
            error('Unable to set RX/TX gain combination to appropriate value');
            
        end
        
        % Perform calibration for XO
        function CalibrateXO(obj, OtherDevice, CenterFrequency)
            
            % Generate Tone
            Fs = obj.SampleRate;
            Fc = 1e6; % Frequency of tone
            fftSize = 2^19;
            Samples = fftSize;
            t = 1/Fs:1/Fs:Samples/Fs;
            txWaveform = complex(sin(2*pi*Fc*t+pi/2)*1024,sin(2*pi*Fc*t)*1024);
            
            % Set up XO info
            XOOffset = 0;
            F_osc = 40e6; % XO Freq
            MULT_lo = CenterFrequency/F_osc;
            %MULT_fs = obj.SampleRate/F_osc;
            %binRez = F_osc*MULT_fs/fftSize;
            
%             if MULT_lo>obj.XOCorrectionFrequencyTolerance
%                 warning('Tolerance may not be attainable');
%             end
            
            %sa = dsp.SpectrumAnalyzer;
            %sa.SampleRate = Fs;
            
            % Determine device set up
            SDRisTX = strcmp(obj.Mode,'TX');%(obj.inputs~=0) && (obj.outputs==0);
            
            % Set up transmitter and receiver
            if SDRisTX
                % Set up this object as transmitter
                obj.CenterFrequency = CenterFrequency;
                obj.Setup(length(txWaveform));
                obj.EnableTX(txWaveform);
                % Set up other object as receiver
                OtherDevice.CenterFrequency = CenterFrequency;
                OtherDevice.Setup(length(txWaveform));
            else
                % Set up other object as transmitter
                OtherDevice.CenterFrequency = CenterFrequency;
                OtherDevice.Setup(length(txWaveform));
                OtherDevice.EnableTX(txWaveform);
                % Set up this object as receiver
                obj.CenterFrequency = CenterFrequency;
                obj.Setup(length(txWaveform)*4);                
            end

            DCXOOffset = 0;
            
            %% Update XO to force tone to specific frequency
            fprintf(['-------------------------\n',...
                'Starting SDR calibration\n']);
            for k=1:obj.XOCorrectionTries
                
                if SDRisTX
                    obj.EnableTX(txWaveform); % Force update to XO
                    RxIQ = OtherDevice.EnableRX();
                else
                    RxIQ = obj.EnableRX();
                end
                % Make sure we have enough data
                RxIQ = [RxIQ;zeros(fftSize-length(RxIQ),1)]; %#ok<AGROW>
                Rx = RxIQ(1:fftSize);
                %sa(Rx);
                % Estimate frequency of signal
                freq = abs(fftshift(fft(Rx)));
                [~,ind] = max(freq);
                ind = ind - fftSize/2;
                %disp( ['Perceived Offset: ',num2str(ind*Fs/fftSize)] );
                
                % Estimate current offset with correct clock
                %currentTrueOffsetHz = ind*binRez;
                %disp( ['Estimated Offset: ',num2str(currentTrueOffsetHz),' Hz'] );
                
                % Correct LO to force tone to be a specific offset
                DesiredFreq = Fc;
                if SDRisTX
                    errorHz = DesiredFreq - ind*Fs/fftSize;
                else
                    errorHz = -DesiredFreq + ind*Fs/fftSize;
                end
                disp( ['Offset Error: ',num2str(errorHz),' Hz'] );
                update = errorHz/MULT_lo; % Error with respect to original clock
                % Check if we are close enough
                if abs(errorHz)<=obj.XOCorrectionFrequencyTolerance
                    disp('Frequency offset within tolerance, calibration complete');
                    break;
                end
                % From new XO calculate actual center frequency and sample rate
                %f_s = XOOffset*MULT_fs;
                %binRez = f_s/fftSize;
                
                if strcmp(obj.dev_name,'fmcomms3')
                    DCXOOffset = DCXOOffset - sign(update);
                    disp(['DCXO Update : ',num2str(DCXOOffset)]);
                    obj.DCXOCorrectionValue = DCXOOffset;
                else
                    XOOffset = XOOffset - round(update);
                    disp(['XO Update Prescale: ',num2str(round(update))]);
                    % Update XO parameter in SDR
                    obj.XOCorrectionValue = XOOffset;
                    %obj.UpdateInputStructField('XO_CORRECTION',40e6 + XOOffset);
                end
            end
            % Save correction value
            obj.XOCorrectionValue = XOOffset;
            % Clean up
            OtherDevice.Release();
            obj.Release();
            
        end
        
        % Cleanup
        function Release(obj)
            % Release any resources used by the system object.
            obj.radioDev.releaseImpl();
%             obj.radioDev.iio_dev_cfg = {};
%             delete(obj.radioDev.libiio_data_in_dev);
%             delete(obj.radioDev.libiio_data_out_dev);
%             delete(obj.radioDev.libiio_ctrl_dev);
            obj.RadioSetup = false;
        end
        
%         function PropListener(obj)
%             addlistener(obj,'CenterFrequency','PostSet',@obj.SetupInputStruct);
%             addlistener(obj,'SampleRate','PostSet',@obj.SetupInputStruct);
%             addlistener(obj,'RFBandwidth','PostSet',@obj.SetupInputStruct);
%             addlistener(obj,'AGCMode','PostSet',@obj.SetupInputStruct);
%             addlistener(obj,'RXGain','PostSet',@obj.SetupInputStruct);
%             addlistener(obj,'TXAttenuation','PostSet',@obj.SetupInputStruct);
%             addlistener(obj,'FrameLength','PostSet',@obj.SetupInputStruct);
%         end
        
    end
    
    methods (Access=protected)
                
        function processTunedPropertiesImpl(obj)
            % Perform calculations if tunable properties change while
            % system is running
            obj.SetupInputStruct();
        end
        
        
    end
    
    methods (Static)
        % Scale Input
        function out = ScaleInputData(in)
            % Scale the signal for better power output and cast to int16. This is the
            % native format for the SDR hardware. Since we are transmitting the same
            % signal in a loop, we can do the cast once to save processing time.
            powerScaleFactor = 0.75;
            in = in.*(1/max(abs(in))*powerScaleFactor);
            out = int16(in*2^15);
        end 
    end
    
end

