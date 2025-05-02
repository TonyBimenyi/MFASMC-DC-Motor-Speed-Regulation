% MATLAB script to read and plot motor data from a file
clear all;
close all;

% Define parameters
file_path = 'dataplot.txt'; % Replace with your file's path

% Load the content of the data
file_content = fileread(file_path);

% Extract the CH3, CH1, CH2, CH4 data using regular expressions
ch3_data_match = regexp(file_content, 'CH3_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');
ch1_data_match = regexp(file_content, 'CH1_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');
ch2_data_match = regexp(file_content, 'CH2_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');
ch4_data_match = regexp(file_content, 'CH4_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');

if ~isempty(ch3_data_match) && ~isempty(ch1_data_match) && ~isempty(ch2_data_match) && ~isempty(ch4_data_match)
    % Convert the extracted data to numerical arrays
    ch3_data = str2num(ch3_data_match{1}{1});
    ch1_data = str2num(ch1_data_match{1}{1});
    ch2_data = str2num(ch2_data_match{1}{1});
    ch4_data = str2num(ch4_data_match{1}{1});
    
    % Calculate distributed errors
    error1 = ch2_data - ch1_data;
    error2 = ch2_data - ch3_data;
    error3 = ch2_data - ch4_data;
    
    % Plot the data
    figure('Position', [100, 100, 800, 400]);
    plot(ch2_data, '-', 'LineWidth', 2, 'Color', 'blue', 'MarkerSize', 4.5, 'DisplayName', 'Reference Trajectory');
    hold on;
    plot(ch1_data, '-.m', 'LineWidth', 2, 'MarkerSize', 5.5, 'DisplayName', 'Motor1');
    plot(ch3_data, '-', 'LineWidth', 2, 'Color', 'red', 'MarkerSize', 2.5, 'DisplayName', 'Motor 2');
    plot(ch4_data, '-.', 'LineWidth', 2, 'Color', 'green', 'MarkerSize', 0.5, 'DisplayName', 'Motor 3');
    hold off;
    
    % Set plot limits and labels
    xlim([0, 1660]);
    xlabel('Time step', 'FontSize', 14);
    ylabel('Output', 'FontSize', 14);
    set(gca, 'FontSize', 12);
    
    % Add legend
    legend('show');
    
    % Save the plot (equivalent to plt.show() in Python)
    print('-dpng', 'motor_plot.png');
else
    disp('CH3_Data_OutPut[1707] not found in the file.');
end