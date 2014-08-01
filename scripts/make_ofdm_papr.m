clear all; close all; clc;

Nsub=16;

done = 0;
file_number = 76;
count = 1;
count_right = 0;

while(done ~= 1)
    posF = 2*rand(1,Nsub)-1;
    %posF = [ 0.6294    0.8116   -0.7460    0.8268    0.2647   -0.8049   -0.4430...
    %          0.0938    0.9150    0.9298   -0.6848    0.9412    0.9143   -0.0292...
    %          0.6006   -0.7162];
%     posF = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
    %posF = ones(1,16) .* -1;
    
    negF = fliplr(posF);
    xF = [0 posF negF];
    xT = ifft(xF);
    shit = fft(xT);
    
    
    
    % Generate and auto-scale triangle buffer to best operating values
    triangle_buffer = [xT zeros(1,1)];
    %triangle_buffer = xT;
    b = max(triangle_buffer);
    a = min(triangle_buffer);
    swing=4050;  %for full swing use this line
    % swing=3300;
    k=swing/(b-a);
    avg=(k*a+k*b)/2;
    triangle_buffer = triangle_buffer*(k)-avg+(swing/2);
    
    %normalize
    test = triangle_buffer ./ max(abs(triangle_buffer));
    test = test - mean(test);
    
    max_squared = max(abs(test)) ^ 2;
    average_squared = mean(test.^2);
    papr = max_squared / average_squared;
    
    if(papr > 30 && papr < 40)
        papr
        count
        figure;
        subplot(2,1,1);
        plot(posF, '-bs', 'LineWidth', 2);
        grid on;
        subplot(2,1,2);
        plot(triangle_buffer, '-bs', 'LineWidth', 2);grid on;
        pause;
        fid = fopen(['different_subcarriers\tx_signal_' num2str(file_number) ...
            '.txt'], 'w');
        for i = 1:length(xT)
            fprintf(fid, '%i, ', round(triangle_buffer(i)));
        end
        file_number = file_number + 1;
        fprintf(fid, '%i', 0);
        count_right = count_right + 1;
        if(count_right > 2)
            done = 1;
        end
    end
    
    %  fprintf('int hardcode[%i] = {',length(xT)+1);
    %  for i = 1:length(xT)
    %      fprintf('%i, ', round(triangle_buffer(i)));
    %  end
    %  fprintf('0};\n')
    count = count + 1;
end

% fid = fopen(['different_subcarriers\tx_signal_' num2str(file_number) ...
%     '.txt'], 'w');
% for i = 1:length(xT)
%     fprintf(fid, '%i, ', round(triangle_buffer(i)));
% end
% fprintf(fid, '%i', 0);
% 
% figure;
% subplot(2,1,1);
% plot(posF, '-bs', 'LineWidth', 2);
% grid on;
% subplot(2,1,2);
% plot(triangle_buffer, '-bs', 'LineWidth', 2);grid on;
%
%  %testing a bunch of stuff
%
%  triangle_buffer = triangle_buffer ./ max(abs(triangle_buffer));
%
%  %remove DC component before fft
%  triangle_buffer = triangle_buffer - mean(triangle_buffer);
%
%  %fft
%  triangle_buffer =(fft(triangle_buffer(1:end-1)));
%
%  %normalize
%  triangle_buffer(1) = 0;
%  triangle_buffer = triangle_buffer ./ max(abs(triangle_buffer));
%
%
%
%  xF = xF ./ max(abs(xF));
%
%  figure;
%  hold on;
%  plot(abs(xF),'-bs','LineWidth',2);
%  plot(abs(triangle_buffer),'-rs','LineWidth',2);