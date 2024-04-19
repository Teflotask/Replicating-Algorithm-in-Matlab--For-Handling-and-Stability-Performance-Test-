
classdef Vehiclestabaltesting < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure             matlab.ui.Figure
        ProcessDataButton    matlab.ui.control.Button
        LoadDataButton       matlab.ui.control.Button
        DropDown_TestType    matlab.ui.control.DropDown
        SelectTestTypeLabel  matlab.ui.control.Label
        axesFilterPlot       matlab.ui.control.UIAxes
        axesScatterPlot      matlab.ui.control.UIAxes
    end

    
    properties (Access = public)
        Datatable % Table that stores uploaded data 
    end
   

    
    methods (Access = public)
        
       function results = func(app)


% Function to perform Snake Test
function runSnakeTest(app)
    % Check if data is loaded
    if isempty(app.Datatable)
        uialert(app.UIFigure, 'Load data before running the test.', 'Data Not Loaded');
        return;
    end

    % Data extraction and calculations
    try
        % data extraction
        steeringAngles = app.Datatable.SteeringAngle; % column name
        yawVelocity = app.Datatable.YawVelocity;       % column name
        rollAngle = app.Datatable.RollAngle;           % column name
        acceleration = app.Datatable.Acceleration;     % column name

        % Compute evaluation indices
        avgSteeringAngle = mean(steeringAngles);
        avgYawVelocity = mean(yawVelocity);
        avgRollAngle = mean(rollAngle);
        avgAcceleration = mean(acceleration);

        % Update plots
        updateScatterPlot(app, steeringAngles, rollAngle);
        updateFilterPlot(app, yawVelocity, acceleration);
        
        % Display results in the app (assuming there are UI components to display these results)
        app.LabelAvgSteeringAngle.Text = ['Average Steering Angle: ', num2str(avgSteeringAngle), ' deg'];
        app.LabelAvgYawVelocity.Text = ['Average Yaw Angular Velocity: ', num2str(avgYawVelocity), ' deg/s'];
        app.LabelAvgRollAngle.Text = ['Average Body Roll Angle: ', num2str(avgRollAngle), ' deg'];
        app.LabelAvgAcceleration.Text = ['Average Acceleration: ', num2str(avgAcceleration), ' m/s^2'];

    catch e
        uialert(app.UIFigure, 'Error processing data. Check data format.', 'Data Error');
        disp(e.message);
    end
end

% Function to update scatter plot (Steering Angle vs. Roll Angle)
function updateScatterPlot(app, xData, yData)
    scatter(app.axesScatterPlot, xData, yData, 'filled');
    title(app.axesScatterPlot, 'Scatter Plot: Steering Angle vs. Roll Angle');
    xlabel(app.axesScatterPlot, 'Steering Angle (deg)');
    ylabel(app.axesScatterPlot, 'Roll Angle (deg)');
end

% Function to update filter plot (Yaw Velocity vs. Acceleration)
function updateFilterPlot(app, xData, yData)
    plot(app.axesFilterPlot, xData, yData);
    title(app.axesFilterPlot, 'Filter Plot: Yaw Velocity vs. Acceleration');
    xlabel(app.axesFilterPlot, 'Yaw Velocity (deg/s)');
    ylabel(app.axesFilterPlot, 'Acceleration (m/s^2)');
end

%STEP INPUT TEST

% Function to perform Step Input Test
function runStepInputTest(app)
    % Check if data is loaded
    if isempty(app.Datatable)
        uialert(app.UIFigure, 'Load data before running the test.', 'Data Not Loaded');
        return;
    end

    % Data extraction and preprocessing
    try
        % columns names
        steeringWheelAngle = app.Datatable.SteeringWheelAngle;
        lateralAcceleration = app.Datatable.LateralAcceleration;
        yawVelocity = app.Datatable.YawVelocity;
        sideSlipAngle = app.Datatable.SideSlipAngle;

        % Preprocessing: Remove outliers, zero drift, and apply low-pass filtering
        steeringWheelAngle = preprocessData(steeringWheelAngle);
        lateralAcceleration = preprocessData(lateralAcceleration);
        yawVelocity = preprocessData(yawVelocity);
        sideSlipAngle = preprocessData(sideSlipAngle);

        % Test calculations
        %  Calculate response times, overshoot, and steady-state values
        [responseTime, overshoot, steadyStateValues] = calculateStepInputMetrics(steeringWheelAngle, lateralAcceleration, yawVelocity);

        % Update plots
        updateScatterPlot(app, yawVelocity, lateralAcceleration);
        updateFilterPlot(app, steeringWheelAngle, sideSlipAngle);
        
        % Display results in the app
        displayResults(app, responseTime, overshoot, steadyStateValues);

    catch e
        uialert(app.UIFigure, 'Error processing data. Check data format.', 'Data Error');
        disp(e.message);
    end
end

function dataFiltered = preprocessData(data)
    % preprocessing steps
    dataFiltered = removeOutliers(data);
    dataFiltered = removeZeroDrift(dataFiltered);
    dataFiltered = lowPassFilter(dataFiltered);
end

function [responseTime, overshoot, steadyState] = calculateStepInputMetrics(steeringAngle, lateralAcc, yawVel)
    % calculation logic
    responseTime = mean(steeringAngle); 
    overshoot = max(lateralAcc) - mean(lateralAcc); 
    steadyState = mean([lateralAcc, yawVel], 'all'); 
end


% Function to update scatter plot
%function updateScatterPlot(app, xData, yData)
    scatter(app.axesScatterPlot, xData, yData, 'filled');
    title(app.axesScatterPlot, 'Scatter Plot: Yaw Velocity vs. Lateral Acceleration');
    xlabel(app.axesScatterPlot, 'Yaw Velocity (deg/s)');
    ylabel(app.axesScatterPlot, 'Lateral Acceleration (m/s^2)');
end

% Function to update filter plot
%function updateFilterPlot(app, xData, yData)
    %plot(app.axesFilterPlot, xData, yData);
    %title(app.axesFilterPlot, 'Filter Plot: Steering Wheel Angle vs. Side Slip Angle');
    %xlabel(app.axesFilterPlot, 'Steering Wheel Angle (deg)');
   % ylabel(app.axesFilterPlot, 'Side Slip Angle (deg)');
%end


function displayResults(app, responseTime, overshoot, steadyStateValues)
    % displaying STEPUP INPUT results
    app.LabelResponseTime.Text = ['Response Time: ', num2str(responseTime)];
    app.LabelOvershoot.Text = ['Overshoot: ', num2str(overshoot)];
    app.LabelSteadyStateValues.Text = ['Steady State Values: ', num2str(steadyStateValues)];
end



%PULSE INPUT TEST FUNCTION

% Function to perform Pulse Input Test
function runPulseInputTest(app)
    % Check if data is loaded
    if isempty(app.Datatable)
        uialert(app.UIFigure, 'Load data before running the test.', 'Data Not Loaded');
        return;
    end

    % Data extraction and preprocessing
    try
        % column names, 
        steeringWheelAngle = app.Datatable.SteeringWheelAngle;
        yawVelocity = app.Datatable.YawVelocity;
        lateralAcceleration = app.Datatable.LateralAcceleration;

        % Preprocess data
        steeringWheelAngle = preprocessData(steeringWheelAngle);
        yawVelocity = preprocessData(yawVelocity);
        lateralAcceleration = preprocessData(lateralAcceleration);

        % Fourier Transform to get frequency response
        [amplitudeFreqSteeringYaw, phaseFreqSteeringYaw] = getFrequencyResponse(steeringWheelAngle, yawVelocity);
        [amplitudeFreqSteeringLat, phaseFreqSteeringLat] = getFrequencyResponse(steeringWheelAngle, lateralAcceleration);

        % Calculate evaluation indicators
        [resonantFrequency, phaseLagAngle, resonancePeakLevel] = evaluateFrequencyResponse(amplitudeFreqSteeringYaw, phaseFreqSteeringYaw);

        % Display results
        displayResults(app, resonantFrequency, phaseLagAngle, resonancePeakLevel);

    catch e
        uialert(app.UIFigure, 'Error processing data. Check data format.', 'Data Error');
        disp(e.message);
    end
end

function dataFiltered = preprocessData(data)
    %  preprocessing logic
    dataFiltered = highpass(data, 0.1); % Example to remove low frequency drift
end

function [amplitude, phase] = getFrequencyResponse(inputSignal, responseSignal)
    % Perform Fourier transform and calculate amplitude and phase
    n = length(inputSignal);
    Y = fft(responseSignal);
    P2 = abs(Y/n);
    amplitude = P2(1:n/2+1);
    P1 = angle(Y/n);
    phase = P1(1:n/2+1);
end

function [freqResonance, phaseLag, peakLevel] = evaluateFrequencyResponse(amplitude, phase)
    % Placeholder: Calculate the evaluation metrics from amplitude and phase data
    [peak, idx] = max(amplitude);
    freqResonance = idx; 
    phaseLag = phase(idx);
    peakLevel = 20*log10(peak); % Convert peak amplitude to decibels
end

%function displayResults(app, freqResonance, phaseLag, peakLevel)
    % Displaying results on the app
    %app.LabelFreqResonance.Text = ['Resonance Frequency: ', num2str(freqResonance), ' Hz'];
    %app.LabelPhaseLag.Text = ['Phase Lag Angle: ', num2str(phaseLag), ' degrees'];
   % app.LabelPeakLevel.Text = ['Resonance Peak Level: ', num2str(peakLevel), ' dB'];
%end

%Steering Return Performance Test


% Function to perform Steering Return Performance Test
function runSteeringReturnPerformanceTest(app)
    % Check if data is loaded
    if isempty(app.Datatable)
        uialert(app.UIFigure, 'Load data before running the test.', 'Data Not Loaded');
        return;
    end

    % Data extraction
    try
        %  data column is named 'SteeringAngle'
        steeringAngle = app.Datatable.SteeringAngle;

        % Preprocess data: filtering and smoothing to reduce noise
        steeringAngleFiltered = preprocessSteeringData(steeringAngle);

        % Calculate performance metrics
        [firstPeak, secondPeak, attenuationRatio] = calculateAttenuation(steeringAngleFiltered);

        % Display results
        displaySteeringTestResults(app, firstPeak, secondPeak, attenuationRatio);

    catch e
        uialert(app.UIFigure, 'Error processing data. Check data format.', 'Data Error');
        disp(e.message);
    end
end

function filteredData = preprocessSteeringData(data)
    % Apply a low-pass filter to smooth the data and remove high frequency noise
    filteredData = lowpass(data, 2);  % Low-pass filter at 2 Hz, adjust based on data sampling rate
end

function [firstPeak, secondPeak, ratio] = calculateAttenuation(data)
    % Find peaks in the steering angle data
    [peaks, locs] = findpeaks(data);

    % Ensure there are at least two peaks to compare
    if length(peaks) >= 2
        firstPeak = peaks(1);
        secondPeak = peaks(2);
        ratio = secondPeak / firstPeak;  % Calculate the ratio of the second peak to the first
    else
        firstPeak = NaN;
        secondPeak = NaN;
        ratio = NaN;
        disp('Not enough peaks found to calculate attenuation.');
    end
end

function displaySteeringTestResults(app, firstPeak, secondPeak, ratio)
    % Displaying results on the app
    app.LabelFirstPeak.Text = ['First Peak: ', num2str(firstPeak)];
    app.LabelSecondPeak.Text = ['Second Peak: ', num2str(secondPeak)];
    app.LabelAttenuationRatio.Text = ['Attenuation Ratio: ', num2str(ratio)];
end





%Steering Lightness Test 


% Function to perform Steering Lightness Test
function runSteeringLightnessTest(app)
    % Check if data is loaded
    if isempty(app.Datatable)
        uialert(app.UIFigure, 'Load data before running the test.', 'Data Not Loaded');
        return;
    end

    % Data extraction
    try
        % data column names
        torqueData = app.Datatable.SteeringWheelTorque;  % Torque at the steering wheel
        steeringWheelRadius = app.SteeringWheelRadius;   % Radius of the steering wheel, might be a constant or input

        % Calculate maximum steering wheel torque
        maxTorque = max(abs(torqueData));

        % Calculate maximum force
        maxForce = maxTorque / steeringWheelRadius;

        % Calculate work done by the steering wheel
        workDone = calculateSteeringWheelWork(torqueData);

        % Calculate average friction torque
        averageFrictionTorque = mean(abs(torqueData));

        % Calculate average friction force
        averageFrictionForce = averageFrictionTorque / steeringWheelRadius;

        % Display results
        displaySteeringLightnessResults(app, maxTorque, maxForce, workDone, averageFrictionTorque, averageFrictionForce);

    catch e
        uialert(app.UIFigure, 'Error processing data. Check data format.', 'Data Error');
        disp(e.message);
    end
end

function work = calculateSteeringWheelWork(torqueData)
    % Calculate the work done by the steering wheel along the path
    % angular displacement for each point is small and can be approximated as a small segment of a circle
    theta = 2 * pi / length(torqueData);  % Total rotation divided by number of samples
    work = sum(torqueData * theta);       % Work = torque * angular displacement
end

function displaySteeringLightnessResults(app, maxTorque, maxForce, work, avgTorque, avgForce)
    % Update GUI components with results
    app.LabelMaxTorque.Text = ['Maximum Torque: ', num2str(maxTorque), ' Nm'];
    app.LabelMaxForce.Text = ['Maximum Force: ', num2str(maxForce), ' N'];
    app.LabelWorkDone.Text = ['Work Done: ', num2str(work), ' J'];
    app.LabelAverageTorque.Text = ['Average Friction Torque: ', num2str(avgTorque), ' Nm'];
    app.LabelAverageForce.Text = ['Average Friction Force: ', num2str(avgForce), ' N'];
end





%Steady-State Rotation Test


% Function to perform Steady-State Rotation Test
function runSteadyStateRotationTest(app, testType)
    % Check if data is loaded
    if isempty(app.Datatable)
        uialert(app.UIFigure, 'Load data before running the test.', 'Data Not Loaded');
        return;
    end

    % Data extraction
    try
        % column names
        longitudinalSpeed = app.Datatable.LongitudinalSpeed;
        lateralAcceleration = app.Datatable.LateralAcceleration;
        steeringWheelAngle = app.Datatable.SteeringWheelAngle;

        % Select the test method based on input
        switch testType
            case 'Fixed Radius'
                testResult = performFixedRadiusTest(longitudinalSpeed, lateralAcceleration);
            case 'Constant Speed'
                testResult = performConstantSpeedTest(longitudinalSpeed, steeringWheelAngle);
            case 'Fixed Angle'
                testResult = performFixedAngleTest(steeringWheelAngle, lateralAcceleration);
            otherwise
                uialert(app.UIFigure, 'Invalid test type selected.', 'Test Error');
                return;
        end

        % Display results
        displayTestResults(app, testResult);

    catch e
        uialert(app.UIFigure, 'Error processing data. Check data format.', 'Data Error');
        disp(e.message);
    end
end

function result = performFixedRadiusTest(longitudinalSpeed, lateralAcceleration)
    % function to analyze data for the Fixed Radius test
    result.maxLateralAcc = max(lateralAcceleration);
    result.avgLongitudinalSpeed = mean(longitudinalSpeed);
end

function result = performConstantSpeedTest(longitudinalSpeed, steeringWheelAngle)
    % function to analyze data for the Constant Speed test
    result.avgSpeed = mean(longitudinalSpeed);
    result.constantSteeringAngle = median(steeringWheelAngle);
end

function result = performFixedAngleTest(steeringWheelAngle, lateralAcceleration)
    %  function to analyze data for the Fixed Angle test
    result.constantSteeringAngle = median(steeringAngle);
    result.maxLateralAcc = max(lateralAcceleration);
end

function displayTestResults(app, results)
    % Update GUI components with results from the chosen test method
    app.LabelMaxLateralAcceleration.Text = ['Max Lateral Acceleration: ', num2str(results.maxLateralAcc), ' m/s^2'];
    app.LabelAverageSpeed.Text = ['Average Speed: ', num2str(results.avgSpeed), ' km/h'];
end




%Steering Wheel Center Area Test



% Function to perform Steering Wheel Center Area Test
function runSteeringWheelCenterAreaTest(app)
    % Check if data is loaded
    if isempty(app.Datatable)
        uialert(app.UIFigure, 'Load data before running the test.', 'Data Not Loaded');
        return;
    end

    % Data extraction
    try
        % column names
        steeringStiffness = app.Datatable.SteeringStiffness;  % Total steering stiffness
        steeringFrictionTorque = app.Datatable.SteeringFrictionTorque;
        steeringAngleHysteresis = app.Datatable.SteeringAngleHysteresis;
        yawAngularSpeed = app.Datatable.YawAngularSpeed;
        yawRateResponse = app.Datatable.YawRateResponse;

        % Analyze data specifically around the center area
        [avgSteeringStiffness, centerSteeringStiffness, steeringFriction, angleHysteresis, yawSpeedGain, yawResponseLag] = analyzeCenterArea(steeringStiffness, steeringFrictionTorque, steeringAngleHysteresis, yawAngularSpeed, yawRateResponse); %#ok<*ADMTHDINV>

        % Display results
        displayCenterAreaResults(app, avgSteeringStiffness, centerSteeringStiffness, steeringFriction, angleHysteresis, yawSpeedGain, yawResponseLag);

    catch e
        uialert(app.UIFigure, 'Error processing data. Check data format.', 'Data Error');
        disp(e.message);
    end
end

function [avgStiffness, centerStiffness, frictionTorque, hysteresis, speedGain, responseLag] = analyzeCenterArea(stiffness, frictionTorque, hysteresis, yawSpeed, yawResponse)
    % Analyze the steering characteristics in the center area
    % analysis logic based on the data characteristics and test requirements
    avgStiffness = mean(stiffness);
    centerStiffness = mean(stiffness(abs(stiffness) < 0.5));  % Assume small deviations as 'center area'
    frictionTorque = mean(frictionTorque);
    hysteresis = mean(hysteresis);
    speedGain = max(yawSpeed) / mean(yawSpeed);
    responseLag = mean(diff(yawResponse));  % calculate lag
end

function displayCenterAreaResults(app, avgStiffness, centerStiffness, frictionTorque, hysteresis, speedGain, responseLag)
    % Update GUI components with results from the center area test
    app.LabelAvgStiffness.Text = ['Average Steering Stiffness: ', num2str(avgStiffness), ' Nm/deg'];
    app.LabelCenterStiffness.Text = ['Center Area Steering Stiffness: ', num2str(centerStiffness), ' Nm/deg'];
    app.LabelFrictionTorque.Text = ['Steering Friction Torque: ', num2str(frictionTorque), ' Nm'];
    app.LabelHysteresis.Text = ['Steering Angle Hysteresis: ', num2str(hysteresis), ' degrees'];
    app.LabelYawSpeedGain.Text = ['Yaw Angular Speed Gain: ', num2str(speedGain)];
    app.LabelResponseLag.Text = ['Response Lag Time: ', num2str(responseLag), ' s'];
end






            
            
 %PREPROCESSDATA Function
function Datatable = prepro(app, Datatable)

% Remove rows with missing values
    Datatable = rmmissing(Datatable);

% Normalize or scale specific columns, e.g., scaling longitudinal speed
    if any(strcmp(Datatable.Properties.VariableNames, 'LongitudinalSpeed'))
        Datatable.LongitudinalSpeed = normalize(Datatable.LongitudinalSpeed, 'range');

    end

% Filter high-frequency noise from specific columns, e.g., steering angle
    if any(strcmp(Datatable.Properties.VariableNames, 'SteeringAngle'))
        Datatable.SteeringAngle = lowpass(Datatable.SteeringAngle, 1);  % Using a 1 Hz cutoff
    end
            
% Replace outliers in lateral acceleration using median value
    if any(strcmp(Datatable.Properties.VariableNames, 'LateralAcceleration'))
        accMedian = median(Datatable.LateralAcceleration, 'omitnan');
        isOutlier = isoutlier(dataTable.LateralAcceleration);
        Datatable.LateralAcceleration(isOutlier) = accMedian;
    end



        end
        end
   % end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: LoadDataButton
        function LoadDataButtonPushed(app, event)
          % Open file dialog with no filter
[files, paths] = uigetfile('*.', 'Select files', 'MultiSelect', 'on');

% Check if user canceled selection
if ~iscell(files)
    disp('No files selected.');
    return;
end

% Handle multiple file selections
numFiles = numel(files);
app.Datatable = cell(1, numFiles);  % Pre-allocate cell array for data

for i = 1:numFiles
    filePath = fullfile(paths{i}, files{i});
    try
        % Attempt to read the file based on its extension
        [data, ~] = readtable(filePath);
        app.Datatable{i} = data;
        disp(['Successfully loaded data from: ', filePath]);
    catch ME
        warning(['Error reading file: ', filePath]);
        warning(ME.message);
    end
end

% Check if any files were loaded successfully
if all(cellfun(@isempty, app.Datatable))
    disp('No data loaded from any files.');

% Preprocess data after loading
        app.Datatable = prepro(app, app.Datatable);
        disp('Data preprocessed successfully.');
 end











        end

        % Value changed function: DropDown_TestType
        function DropDown_TestTypeValueChanged(app, event)
              testType = app.DropDown_TestType.Value; 

    % Clear existing data and plots
    %clearDataAndPlots(app);

    % Update UI or internal variables based on selected test
    switch testType
        case 'Snake Test'
            runSnakeTest(app);
        case 'Step Input Test'
            runStepInputTest(app);
        case 'Pulse Input Test'
            runPulseInputTest(app);
        case 'Steering Return Performance Test'
            runSteeringReturnPerformanceTest(app);
        case 'Steering Lightness Test'
            runSteeringLightnessTest(app);
        case 'Steady-State Rotation Test'
            runSteadyStateRotationTest(app);
        case 'Steering Wheel Center Area Test'
            runSteeringWheelCenterAreaTest(app);
        otherwise
            % Handle unexpected case
    %end
       % end

        %function clearDataAndPlots(app)
    % Implement functionality to clear data and plots
    %disp('Clearing data and plots...');
       end
        end

        % Button pushed function: ProcessDataButton
        function ProcessDataButtonPushed(app, event)
   % Callback function for the "Process Data" button
%function ProcessButtonPushed(app, event)
    % Check if the Datatable is loaded
    if isempty(app.Datatable)
        % If not, alert the user and exit the function
        uialert(app.UIFigure, 'No data loaded. Please load data before running tests.', 'Data Not Loaded');
        return;
    end

    % Get the selected test type from the dropdown
    testType = app.DropDown_TestType.Value;

    % Run the appropriate test based on the selected value
    switch testType
        case 'Snake Test'
            runSnakeTest(app);  % Call the function to run Snake Test
        case 'Step Input Test'
            runStepInputTest(app);  % Call the function to run Step Input Test
        case 'Pulse Input Test'
            runPulseInputTest(app);  % Call the function to run Pulse Input Test
        case 'Steering Return Performance Test'
            runSteeringReturnPerformanceTest(app);  % Call the function for Steering Return Performance Test
        case 'Steering Lightness Test'
            runSteeringLightnessTest(app);  % Call the function for Steering Lightness Test
        case 'Steady-State Rotation Test'
            runSteadyStateRotationTest(app);  % Call the function for Steady-State Rotation Test
        case 'Steering Wheel Center Area Test'
            runSteeringWheelCenterAreaTest(app);  % Call the function for Steering Wheel Center Area Test
        otherwise
            uialert(app.UIFigure, 'Test not recognized. Please select a valid test.', 'Test Error');
    end
%end


        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1190 683];
            app.UIFigure.Name = 'MATLAB App';

            % Create axesScatterPlot
            app.axesScatterPlot = uiaxes(app.UIFigure);
            title(app.axesScatterPlot, 'Scatter Plot')
            xlabel(app.axesScatterPlot, 'X')
            ylabel(app.axesScatterPlot, 'Y')
            zlabel(app.axesScatterPlot, 'Z')
            app.axesScatterPlot.Position = [66 159 457 325];

            % Create axesFilterPlot
            app.axesFilterPlot = uiaxes(app.UIFigure);
            title(app.axesFilterPlot, 'Filter Plot')
            xlabel(app.axesFilterPlot, 'X')
            ylabel(app.axesFilterPlot, 'Y')
            zlabel(app.axesFilterPlot, 'Z')
            app.axesFilterPlot.Color = [0.9412 0.9412 0.9412];
            app.axesFilterPlot.Position = [555 147 546 337];

            % Create SelectTestTypeLabel
            app.SelectTestTypeLabel = uilabel(app.UIFigure);
            app.SelectTestTypeLabel.BackgroundColor = [0.0745 0.6235 1];
            app.SelectTestTypeLabel.HorizontalAlignment = 'center';
            app.SelectTestTypeLabel.FontName = 'Arial';
            app.SelectTestTypeLabel.FontWeight = 'bold';
            app.SelectTestTypeLabel.FontColor = [1 1 1];
            app.SelectTestTypeLabel.Position = [84 614 113 34];
            app.SelectTestTypeLabel.Text = 'Select Test Type';

            % Create DropDown_TestType
            app.DropDown_TestType = uidropdown(app.UIFigure);
            app.DropDown_TestType.Items = {'Snake Test', 'Steering Lightness Test', 'Pulse Input Test', 'Step Input Test', 'Steering Return Performance Test ', 'Steady-State Rotation Test ', 'Steering Wheel Center Area Test'};
            app.DropDown_TestType.ItemsData = {'Snake Test', 'Step Input Test', 'Pulse Input Test ', 'Steering Return Performance Test ', 'Steering Lightness Test', 'Steady-State Rotation Test', 'Steering Wheel Center Area Test'};
            app.DropDown_TestType.ValueChangedFcn = createCallbackFcn(app, @DropDown_TestTypeValueChanged, true);
            app.DropDown_TestType.FontName = 'Arial Black';
            app.DropDown_TestType.FontSize = 14;
            app.DropDown_TestType.FontWeight = 'bold';
            app.DropDown_TestType.FontColor = [1 1 1];
            app.DropDown_TestType.BackgroundColor = [0 0 0];
            app.DropDown_TestType.Position = [196 614 309 34];
            app.DropDown_TestType.Value = 'Snake Test';

            % Create LoadDataButton
            app.LoadDataButton = uibutton(app.UIFigure, 'push');
            app.LoadDataButton.ButtonPushedFcn = createCallbackFcn(app, @LoadDataButtonPushed, true);
            app.LoadDataButton.Position = [115 66 145 47];
            app.LoadDataButton.Text = 'Load Data';

            % Create ProcessDataButton
            app.ProcessDataButton = uibutton(app.UIFigure, 'push');
            app.ProcessDataButton.ButtonPushedFcn = createCallbackFcn(app, @ProcessDataButtonPushed, true);
            app.ProcessDataButton.Position = [833 66 145 47];
            app.ProcessDataButton.Text = 'Process Data';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Vehiclestabaltesting

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
