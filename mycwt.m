function TF_out = mycwt(Data,freq_range,cycles,ResFreq,srate)

% data
% freq_range [f1 f2]
% cycles
% resFreq: frequnecy resolution eg. 1
% srate: sampling rate

frex = freq_range(1) : ResFreq : freq_range(2);
min_freq = min(frex);
max_freq = max(frex);
num_frex = numel(frex);
s = logspace(log10(cycles(1)),log10(cycles(end)),num_frex) ./ (2*pi*frex);

[nbchan, pnts] = size(Data);
% initialize output time-frequency data

    TF_out = struct();
    tf_cwt = zeros(nbchan,num_frex,pnts);
    tf_amp = zeros(nbchan,num_frex,pnts);
    tf_phase = zeros(nbchan,num_frex,pnts);
    Niter = nbchan*length(frex);
    cont_loop = 1;
    for chi = 1: nbchan
        % loop over frequencies
        for fi=1:length(frex)
            cont_loop = cont_loop + 1;
            waitfillingbar(cont_loop, Niter, 'In Progress');

            wavtime = -2 : 1/srate : 2;
            half_wave = (length(wavtime)-1)/2;
            
            % FFT parameters
            nWave = length(wavtime);            % create wavelet and get its FFT
            nConv = nWave + numel(Data(chi,:)) - 1;
            wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
            waveletX = fft(wavelet,nConv);
            waveletX = waveletX ./ max(waveletX);
            dataX = fft(zscore(Data(chi,:)), nConv);
            
            % run convolution
            as = ifft(waveletX .* dataX);
            as = as(half_wave+1:end-half_wave);
            % put tf data into big matrix
            tf_cwt(chi,fi,:) = 2*abs(as).^2;
            tf_phase(chi,fi,:) = angle(as);
            tf_amp(chi,fi,:) = abs(as);
        end
    end
    TF_out.min_freq = min_freq;
    TF_out.max_freq = max_freq;
    TF_out.num_frex = num_frex;
    TF_out.frex = frex;
    TF_out.range_cycles = cycles;
    TF_out.cwt = tf_cwt;
    TF_out.amp = tf_amp;
    TF_out.phase = tf_phase;
end
