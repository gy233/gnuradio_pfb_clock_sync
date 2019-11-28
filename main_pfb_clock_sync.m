% book: Multirate Signal Processing for Communications Systems
% chapter 13.3, 13.4

% use with /home/guyu/my_gnuradio_projects/cma/verify_block2.grc
clc;clear;close all;

global d_taps_per_filter
global d_dtaps_per_filter

%% RRC
% firdes.root_raised_cosine
samples_per_symbol=2;
nfilts=32;

gain=nfilts;
sampling_freq=nfilts*samples_per_symbol;
symbol_rate=1;
alpha=0.35;
ntaps = nfilts * 11 * samples_per_symbol;    % make nfilts filters of ntaps each
[taps,t,num,den] = root_raised_cosine(gain,sampling_freq,symbol_rate,alpha,ntaps);

% RaisedCosineTransmitFilter
Nsym=11*samples_per_symbol;
sampsPerSym=nfilts;
rctFilt = comm.RaisedCosineTransmitFilter(...
    'Shape', 'Square root', ...
    'RolloffFactor', alpha, ...
    'FilterSpanInSymbols', Nsym, ...
    'OutputSamplesPerSymbol', sampsPerSym);
b1=coeffs(rctFilt);

rcrFilt = comm.RaisedCosineReceiveFilter(...
  'Shape',                  'Square root', ...
  'RolloffFactor',          alpha, ...
  'FilterSpanInSymbols',    Nsym/2, ...
  'InputSamplesPerSymbol',  sampsPerSym*samples_per_symbol, ...
  'DecimationFactor',       1);
b2=coeffs(rcrFilt);

figure
plot(b1.Numerator,'b')
hold on
plot(b2.Numerator,'g','LineWidth',2)
plot(taps,'r')
grid on
legend('RaisedCosineTransmitFilter','RaisedCosineReceiveFilter','firdes.root raised cosine')

%% verify block data
data_length=5000;
filename='pfb_clock_sync_in.bin';
[fid]=fopen(filename,'rb');
input=fread(fid,samples_per_symbol*2*data_length*2,'float32');
input=input(1:2:end)+1i*input(2:2:end);
input=[zeros(d_taps_per_filter+3,1)+1i*zeros(d_taps_per_filter+3,1); input];
fclose(fid);

filename='pfb_clock_sync_out.bin';
[fid]=fopen(filename,'rb');
output_verify=fread(fid,2*data_length,'float32');
output_verify=output_verify(1:2:end)+1i*output_verify(2:2:end);
fclose(fid);

filename='pfb_clock_sync_out_phase.bin'; % d_k
[fid]=fopen(filename,'rb');
pfb_clock_sync_out_phase=fread(fid,data_length,'float32');
fclose(fid);

filename='pfb_clock_sync_out_error.bin'; % d_error
[fid]=fopen(filename,'rb');
pfb_clock_sync_out_error=fread(fid,data_length,'float32');
fclose(fid);

filename='pfb_clock_sync_out_rate.bin'; % d_rate_f
[fid]=fopen(filename,'rb');
pfb_clock_sync_out_rate=fread(fid,data_length,'float32');
fclose(fid);

figure
subplot(121)
scatter(real(input),imag(input));
title('input (constellation)')
subplot(122)
scatter(real(output_verify),imag(output_verify));
title('output verify (constellation)')

%% pfb_clock_sync

%% init pfb_clock_sync
d_max_dev = 1.5;
d_error=0;
init_phase = floor(nfilts/2);
d_out_idx = 0;

d_nfilters = 32;
d_sps = floor(samples_per_symbol);

% Set the damping factor for a critically damped system
d_damping = 2*d_nfilters;

% Set the bandwidth, which will then call update_gains()
d_loop_bw = 2*pi/100;
[d_alpha, d_beta] = update_gains(d_damping, d_loop_bw);

% Store the last filter between calls to work
% The accumulator keeps track of overflow to increment the stride correctly.
% set it here to the fractional difference based on the initial phase
d_k = init_phase;
d_rate = (samples_per_symbol-floor(samples_per_symbol))*d_nfilters;
d_rate_i = floor(d_rate);
d_rate_f = d_rate - d_rate_i;
d_filtnum = floor(d_k);

dtaps = [0,taps(3:end)-taps(1:end-2),0];
dtaps = dtaps * d_nfilters/sum(abs(dtaps));
d_taps=create_taps(taps,0);
d_dtaps=create_taps(dtaps,1);

%% pfb_clock_sync_ccf_impl::general_work
noutput_items = length(output_verify);
d_k_rec=[];
d_error_rec=[];
d_rate_f_rec=[];
i = 1;
count = 0;
output = zeros(size(output_verify));


while(i <= noutput_items)
    
    d_filtnum = floor(d_k);
    
    % Keep the current filter number in [0, d_nfilters]
    % If we've run beyond the last filter, wrap around and go to next sample
    % If we've gone below 0, wrap around and go to previous sample
    while(d_filtnum >= d_nfilters)
        d_k = d_k - d_nfilters;
        d_filtnum = d_filtnum - d_nfilters;
        count = count + 1;
    end
    while(d_filtnum < 0)
        d_k = d_k + d_nfilters;
        d_filtnum = d_filtnum + d_nfilters;
        count = count - 1;
    end
    output(i) = fliplr(d_taps(d_filtnum+1,:))*input(count+1:count+d_taps_per_filter);
    d_k = d_k + d_rate_i + d_rate_f; % update phase
    
    % record
    d_error_rec=[d_error_rec;d_error];
    d_rate_f_rec=[d_rate_f_rec;d_rate_f];
    d_k_rec=[d_k_rec;d_k];
    
    % Update the phase and rate estimates for this symbol
    diff = fliplr(d_dtaps(d_filtnum+1,:)) * input(count+1:count+d_dtaps_per_filter);
    error_r = real(output(i)) * real(diff);
    error_i = imag(output(i)) * imag(diff);
    d_error = (error_i + error_r) / 2.0;       % average error from I&Q channel
    
    
    % Run the control loop to update the current phase (k) and
    % tracking rate estimates based on the error value
    % Interpolating here to update rates for ever sps.
    for s = 1 : d_sps
        d_rate_f = d_rate_f + d_beta*d_error;
        d_k = d_k + d_rate_f + d_alpha*d_error;
    end
    
    % Keep our rate within a good range
    d_rate_f = branchless_clip(d_rate_f, d_max_dev);
    
    i = i + 1;
    count = count + floor(d_sps);
end

d_error_rec=[d_error_rec,pfb_clock_sync_out_error];
d_rate_f_rec=[d_rate_f_rec,pfb_clock_sync_out_rate];
d_k_rec=[d_k_rec,pfb_clock_sync_out_phase];

%% figures
figure
subplot(121)
scatter(real(output_verify),imag(output_verify))
title('GNU radio')
subplot(122)
scatter(real(output),imag(output))
title('mine')
suptitle('scatter')

figure
plot(real(output_verify),'b','LineWidth',2)
hold on
plot(real(output),'r')
legend('GNU radio','mine')

