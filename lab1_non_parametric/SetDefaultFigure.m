function SetDefaultFigure(DrawMode, MonitorNumber)

% function SetDefaultFigure(DrawMode, MonitorNumber)
%
% Future figures are by default drawn in the desired mode (DrawMode) on the desired monitor (MonitorNumber).
%
% Arguments:
% - DrawMode: {'maximized', 'minimized'}; default: 'maximized'
% - MonitorNumber: number of the monitor;
%                  default: last monitor (in the extended monitor mode the screen used for projection)


set(groot, 'Units', 'pixels');        % graphics root object
MP = get(groot, 'MonitorPositions');  % Number of monitors x 4
N = size(MP, 1);                      % Number of monitors

if nargin < 2
    MonitorNumber = N;
end

if nargin < 1
    DrawMode = 'maximized';
end

if ((nargin > 2) || (length(MonitorNumber) > 1) || (round(MonitorNumber) < 1)  || (round(MonitorNumber) > N))
    error('Wrong argument.');
end

MonitorNumber = round(MonitorNumber);

h = figure('Units', 'pixels');        % new figure
set(h, 'OuterPosition', MP(MonitorNumber,:))
v=ver ('matlab');
if str2num(v.Version) > 9.39
    % This works in R0218a and above. In lower versions the resolution of
    % the screen is selected
    set(h, 'WindowState', DrawMode);      % maximize the figure window
end
drawnow                               % without this the following line produces faulty numbers for figure position
Pos = get(h, 'Position');
close(h);

set(0, 'defaultFigurePosition', Pos);
