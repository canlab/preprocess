function hhh = vline(x, in1, in2, varargin)
    % function h = vline(x, linetype, label, [extra line args])
    %
    % Draws a vertical line on the current axes at the location specified by 'x'.  Optional arguments are
    % 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
    % label appears in the same color as the line.
    %
    % The line is held on the current axes, and after plotting the line, the function returns the axes to
    % its prior hold state.
    %
    % The HandleVisibility property of the line object is set to "off", so not only does it not appear on
    % legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
    % return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be
    % overridden by setting the root's ShowHiddenHandles property to on.
    %
    % h = vline(42, 'g', 'The Answer')
    %
    % returns a handle to a green vertical line on the current axes at x = 42, and creates a text object on
    % the current axes, close to the line, which reads "The Answer".
    %
    % vline also supports vector inputs to draw multiple lines at once.  For example, 
    %
    % vline([4 8 12], {'g', 'r', 'b'}, {'l1', 'lab2', 'LABELC'})
    %
    % draws three lines with the appropriate labels and colors.
    %
    % By Brandon Kuczenski for Kensington Labs.
    % brandon_kuczenski@kensingtonlabs.com
    % 8 November 2001
    %
    % Modified by Matthew Davidson to accept extra line arguments
    % 12 February 20008
    % 
    % E.g., to set a black dashed line at x=12 with a width of 2 points:
    % vline(12, 'k--', [], 'LineWidth', 2)

    if length(x) > 1  % vector input
        for i = 1:length(x)
            switch nargin
                case 1
                    linetype = 'r:';
                    label = '';
                case 2
                    if ~iscell(in1)
                        in1 = {in1};
                    end
                    if i>length(in1)
                        linetype = in1{end};
                    else
                        linetype = in1{i};
                    end
                    label = '';
                otherwise
                    if ~iscell(in1)
                        in1 = {in1};
                    end
                    if ~iscell(in2)
                        in2 = {in2};
                    end
                    if i > length(in1)
                        linetype = in1{end};
                    else
                        linetype = in1{i};
                    end
                    if i > length(in2)
                        label = in2{end};
                    else
                        label = in2{i};
                    end
            end
            h(i) = vline(x(i), linetype, label, varargin{:});
        end
    else
        switch nargin
            case 1
                linetype = 'r:';
                label = '';
            case 2
                linetype = in1;
                label = '';
            otherwise
                linetype = in1;
                label = in2;
        end


        g = ishold(gca);
        hold on

        y = get(gca, 'ylim');
        h = plot([x x], y, linetype, varargin{:});
        if ~isempty(label)
            xx = get(gca, 'xlim');
            xrange = xx(2)-xx(1);
            xunit = (x-xx(1))/xrange;
            if xunit < 0.8
                text(x+0.01*xrange, y(1)+0.1*(y(2)-y(1)), label, 'color', get(h, 'color'))
            else
                text(x-.05*xrange, y(1)+0.1*(y(2)-y(1)), label, 'color', get(h, 'color'))
            end
        end

        if g ==  0
            hold off
        end
        set(h, 'tag', 'vline', 'handlevisibility', 'off')
    end

    if nargout
        hhh = h;
    end
end