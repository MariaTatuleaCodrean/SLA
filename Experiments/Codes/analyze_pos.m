%% process_velocities_from_tracks.m
% Process tracked_pos.txt and manual_pos.csv, then plot velocity against time.
clear; clc;

%% User settings
% Units: pixels per micrometre for manual_pos.csv x/y coordinates.
manualScale_pixelsPerMicron = 9.95;

% Set these frame numbers for manual track
firstFrame2 = 372; fusionFrame = 803;

% Frame rate
fps = 30; 

%% Locate files
manualFile = fullfile('../Data/manual_pos.csv');
trackFile  = fullfile('../Data/tracked_pos.txt');

%% Load and process tracked_pos.txt
[x, y, z] = parseTrackFile(trackFile);

n = numel(x);
imNumber = (0:n-1)';
time = imNumber ./ fps;

windowSize = max(1, floor(fps / 2));
smooth_x = trailingMeanLikePandas(x, windowSize);
smooth_y = trailingMeanLikePandas(y, windowSize);
smooth_z = trailingMeanLikePandas(z, windowSize);

derivativeStep = 1./fps;  % derivative step, in seconds.

vel_x = smoothDerivative(smooth_x, derivativeStep);
vel_y = smoothDerivative(smooth_y, derivativeStep);
vel_z = smoothDerivative(smooth_z, derivativeStep);
vel = sqrt(vel_x.^2 + vel_y.^2 + vel_z.^2);

trackData = table(imNumber, time, x, y, z, smooth_x, smooth_y, smooth_z, ...
    vel_x, vel_y, vel_z, vel);

%% Load and process manual_pos.csv
manualData = readtable(manualFile);

imageCol = findTableColumn(manualData, {'Image', 'image', 'Frame', 'frame', 'imNumber'});
xCol     = findTableColumn(manualData, {'x', 'X'});
yCol     = findTableColumn(manualData, {'y', 'Y'});

manualFrame = double(manualData.(imageCol));
manualXpix  = double(manualData.(xCol));
manualYpix  = double(manualData.(yCol));
manualTime  = manualFrame ./ fps;

% Manual particle velocity in the image frame, converted from pixels/s to um/s.
manualXvel_pixelsPerSecond = [NaN; diff(manualXpix) ./ diff(manualTime)];
manualYvel_pixelsPerSecond = [NaN; diff(manualYpix) ./ diff(manualTime)];
manualXvel = manualXvel_pixelsPerSecond ./ manualScale_pixelsPerMicron;
manualYvel = manualYvel_pixelsPerSecond ./ manualScale_pixelsPerMicron;

% Add the tracked particle x/y velocity to express the manual particle velocity
% in the lab frame. This aligns by time/frame using interpolation, which is
% safer than relying on row numbers.
trackXvelAtManualFrames = interp1(trackData.time, trackData.vel_x, manualTime, 'linear', NaN);
trackYvelAtManualFrames = interp1(trackData.time, trackData.vel_y, manualTime, 'linear', NaN);

manualXvelReal = manualXvel + trackXvelAtManualFrames;
manualYvelReal = manualYvel + trackYvelAtManualFrames;
manualRealVel  = sqrt(manualXvelReal.^2 + manualYvelReal.^2);

manualSmoothWindow = max(1, floor(fps / 5));
trackSmoothWindow  = max(1, floor(fps));
manualRealVelSmooth = trailingMeanLikePandas(manualRealVel, manualSmoothWindow);
trackVelSmooth      = trailingMeanLikePandas(trackData.vel, trackSmoothWindow);

%% Save time and velocity values for each particle separately

% Particle 1: manually tracked particle
particle1Out = table( ...
    manualTime(:), ...
    manualRealVel(:), ...
    manualRealVelSmooth(:), ...
    'VariableNames', {'time_s', 'velocity_um_per_s', 'velocity_smooth_um_per_s'} ...
);

% Particle 2: particle from tracked_pos.txt
particle2Out = table( ...
    trackData.time(:), ...
    trackData.vel(:), ...
    trackVelSmooth(:), ...
    'VariableNames', {'time_s', 'velocity_um_per_s', 'velocity_smooth_um_per_s'} ...
);

% Optional: remove rows where velocity could not be calculated
particle1Out = rmmissing(particle1Out);
particle2Out = rmmissing(particle2Out);

% Save as separate CSV files
writetable(particle2Out, fullfile('Output/particle1_track_velocity.csv'));
writetable(particle1Out, fullfile('Output/particle2_manual_velocity.csv'));

%% Plot velocities against time
figure('Color', 'w');
hold on;

plot(manualTime, manualRealVelSmooth, 'g', 'LineWidth', 1.5, 'DisplayName', 'Manual particle');
plot(trackData.time, trackVelSmooth, 'b', 'LineWidth', 1.5, 'DisplayName', 'Tracked particle');

hold off;
grid on;
box on;
legend('Location', 'best');
xlabel('Time (s)');
ylabel('Velocity (\mum s^{-1})');
title('Velocities against time');

%% Local helper functions
function [x, y, z] = parseTrackFile(trackFile)

fid = fopen(trackFile, 'r');

x = [];
y = [];
z = [];

while ~feof(fid) % until reached end of file
    line = fgetl(fid); % read next line

    nums = regexp(line, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');

    if numel(nums) >= 4
        x(end+1, 1) = str2double(nums{2});
        y(end+1, 1) = str2double(nums{3});
        z(end+1, 1) = str2double(nums{4});
    end
end

fclose(fid);

end

function derivative = smoothDerivative(v, stepSize)
    v = v(:);
    n = numel(v);
    derivative = zeros(n, 1);

    if n < 5
        derivative(:) = NaN;
        return;
    end

    for k = 3:n-2
        derivative(k) = (2 .* (v(k+1) - v(k-1)) + v(k+2) - v(k-2)) ./ (8 .* stepSize);
    end
end

function out = trailingMeanLikePandas(v, windowSize)
    v = v(:);
    out = NaN(size(v));
    windowSize = max(1, floor(windowSize));

    for k = windowSize:numel(v)
        segment = v(k-windowSize+1:k);
        if all(isfinite(segment))
            out(k) = mean(segment);
        end
    end
end

function name = findTableColumn(tbl, candidates)
    vars = tbl.Properties.VariableNames;
    name = '';

    for i = 1:numel(candidates)
        idx = strcmp(vars, candidates{i});
        if any(idx)
            name = vars{find(idx, 1)};
            return;
        end
    end

    for i = 1:numel(candidates)
        idx = strcmpi(vars, candidates{i});
        if any(idx)
            name = vars{find(idx, 1)};
            return;
        end
    end

    error('Could not find any of these columns in manual_pos.csv: %s', strjoin(candidates, ', '));
end
