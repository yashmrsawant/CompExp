% Noise-Vocoding Simulation in MATLAB
% Load speech and vocode using N channels

clear; close all; clc;

% Load a speech signal
% Use your own .wav file or MATLAB's built-in audio
[signal, fs] = audioread('speech_sample.wav'); % <-- replace with your file
signal = signal(:,1); % mono

% Parameters
N = 6;                    % Number of channels
low_freq = 300;           % Lower frequency bound (Hz)
high_freq = 8000;         % Upper frequency bound (Hz)
order = 4;                % Filter order
t = (0:length(signal)-1)/fs;

% Filter bank cutoff frequencies (log-spaced)
edges = logspace(log10(low_freq), log10(high_freq), N+1);

% Initialize output
vocoded = zeros(size(signal));

% Loop over bands
for n = 1:N
    % Define band edges
    f1 = edges(n);
    f2 = edges(n+1);
    
    % Bandpass filter for current band
    [b, a] = butter(order, [f1, f2]/(fs/2), 'bandpass');
    band = filter(b, a, signal);
    
    % Envelope extraction: rectification + low-pass filter
    envelope = abs(band);
    [blp, alp] = butter(2, 50/(fs/2)); % 50 Hz LPF
    envelope = filter(blp, alp, envelope);
    
    % Generate band-limited noise
    noise = randn(size(signal));
    noise = filter(b, a, noise);  % same bandpass as speech
    
    % Modulate noise with envelope
    modulated = envelope .* noise;
    
    % Sum into final output
    vocoded = vocoded + modulated;
end

% Normalize
vocoded = vocoded / max(abs(vocoded));

% Play original and vocoded
disp('Playing original...');
soundsc(signal, fs); pause(length(signal)/fs + 1);
disp('Playing vocoded...');
soundsc(vocoded, fs);

% Optional: Save vocoded signal
audiowrite(['vocoded_' num2str(N) '_channels.wav'], vocoded, fs);
