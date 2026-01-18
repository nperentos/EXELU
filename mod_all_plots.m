function mod_all_plots(cmd)
% MOD_ALL_PLOTS Apply a command to all subplots (axes) in the current figure
%
% Usage:
%   mod_all_plots('clim([0 150])')
%   mod_all_plots('axis tight')
%   mod_all_plots('set(gca,''YDir'',''normal'')')
%
% Input:
%   cmd - string or char array containing a valid MATLAB command
%         that operates on the current axes (gca)

if nargin < 1
    error('You must provide a command string.');
end

fig = gcf;
ax = findall(fig, 'Type', 'axes');

for k = 1:numel(ax)
    axes(ax(k)); %#ok<LAXES>  % make this axes current
    eval(cmd);
end
end