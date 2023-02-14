[sample_voise,Fs] = audioread('dove.mp3');
%sound(sample_voise,2*Fs);
%time domain
time_values = linspace(0,length(sample_voise)/Fs,length(sample_voise));
%plot(time_values,sample_voise);
%xlabel('Time (ms)');
%ylabel('Amplitude');
%title('Audio Sample in Time Domain');

%frequency domain
Y=fftshift(fft(sample_voise));
magnitudeY = abs(Y);        % Magnitude of the FFT
phaseY = unwrap(angle(Y));  % Phase of the FFT
dB_mag=mag2db(magnitudeY);
Frequency_time_values = linspace(-Fs/2,Fs/2,length(magnitudeY));

%Magnitude plot 
% plot(Frequency_time_values, magnitudeY);
% title('Magnitude response of signal');
% ylabel('Magnitude');
 
%Phase plot 
% plot(Frequency_time_values, phaseY);
% title('Phase Response of the Audio Signal');
% xlabel('Frequency in kHz')
% ylabel('radians');


%Channels responses
% passing the sound signal over delta function 
delta_func= 1;
delta_func_sample=linspace(-Fs/2,Fs/2,length(delta_func));
delta_sound_1= conv(sample_voise(:,1), delta_func);
delta_sound_2= conv(sample_voise(:,2), delta_func);
delta_signa = [delta_sound_1; delta_sound_2];

% plot(delta_sound_1(1:end/2));
% title('Convolution between delta and sound signals');
% xlabel('Time')
% ylabel('Amplitude');

% passing the sound signal over exp(-2pi*5000t) 
sec_func= exp(-2*pi*5000*time_values) ;
sec_func_sample=linspace(-Fs/2,Fs/2,length(sec_func));
sec_conv_signal_1= conv(sample_voise(:,1), sec_func);
sec_conv_signal_2= conv(sample_voise(:,2), sec_func);
sec_func_signal = [sec_conv_signal_1; sec_conv_signal_2];

% plot(sec_conv_signal_1(1:end/2));
% title('Convolution between exp(-2pi*5000t) and sound signals');
% xlabel('Time')
%   ylabel('Amplitude');



% passing the sound signal over exp(-2pi*1000t)
third_func= exp(-2*pi*1000*time_values);
third_func_sample=linspace(-Fs/2,Fs/2,length(third_func));
third_conv_signal_1= conv(sample_voise(:,1),third_func);
third_conv_signal_2= conv(sample_voise(:,1),third_func);
third_func_signal = [third_conv_signal_1; third_conv_signal_2];

% plot(third_conv_signal_1(1:end/2));
% title('Convolution between exp(-2pi*1000t) and sound signals');
% xlabel('Time')
%   ylabel('Amplitude');


%passing the sound signal over impulse response
imp = [1; zeros(1,1)];
b = 2;
a = [1 -0.5];
impulse_func = filter(b,a,imp);
%stem(0:1,impulse_func)
impulse_conv_signal= conv(sample_voise(:,1),impulse_func);
% 
% plot(impulse_conv_signal);
% title('Convolution between impulse response and sound signals');
% xlabel('Time');
% ylabel('Amplitude');


%Adding Noise
sigma = 0.8;

%Generaing the noise
noise = sigma*randn(size(sample_voise));

%Adding noise to the signal
voice_noisy = sample_voise + noise;

%Play Sound with noise
%sound(voice_noisy, 2*Fs); 

%time domain
time_noise_values = linspace(0,length(voice_noisy)/Fs,length(voice_noisy));
% plot(time_noise_values,voice_noisy);
% xlabel('Time (ms)');
% ylabel('Amplitude');
% title('Audio & Noise Sample in Time Domain');

%frequency domain
Z=fftshift(fft(voice_noisy));
magnitudeZ = abs(Z);        % Magnitude of the FFT
phaseZ = unwrap(angle(Z));  % Phase of the FFT
dB_magZ=mag2db(magnitudeZ);
Frequency_noise_values = linspace(-Fs/2,Fs/2,length(magnitudeZ));


%Magnitude in dB plot 
% plot(Frequency_noise_values, dB_magZ);
% title('Magnitude response of Noise');
% ylabel('Magnitude(dB)');


%Phase plot 
% plot(Frequency_noise_values, phaseZ);
% title('Phase response of Noise');
% xlabel('Frequency in kHz')
% ylabel('radians');

%Receiver:
n = length(magnitudeZ);
sampPerFreq = int64(n/Fs);
limit = (sampPerFreq*(Fs/2 - 3400));
magnitudeZ([1:limit n-limit+1:end])=0; 
real_Z=real(ifft(ifftshift(magnitudeZ)));
%Play the sound file after the filter
%sound(real_Z); 

%Plot the output sound file in time domain
% plot(real_Z);
% xlabel('Time (ms)');
% ylabel('Amplitude');
% title('LPF in Time Domain');

%Plot the output sound file in the frequency domain
%frequency domain
 Y=ifftshift(ifft(magnitudeZ));
 magnitudeY = abs(Y);        % Magnitude of the FFT
 phaseY = unwrap(angle(Y));  % Phase of the FFT
 dB_mag=mag2db(magnitudeY);
 Frequency_noise_values1 = linspace(-Fs/2,Fs/2,length(Frequency_noise_values));

% % %Magnitude plot 
 plot(Frequency_noise_values1, magnitudeY);
 title('Magnitude response of Noise & Audio Signals after the LPF');
 ylabel('Magnitude');
 
%Phase plot 
% plot(Frequency_noise_values1, phaseY);
% title('Phase response of Noise & Audio Signals after the LPF');
% xlabel('Frequency in kHz')
% ylabel('radians');
