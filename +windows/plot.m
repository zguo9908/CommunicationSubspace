function matlabGraphicsObjects = plot(windowList, varargin)
% Encapsulates what we repeatedly do with window-plotting. Useful to have
% as a function tool since we do this so often.

ip = inputParser;
ip.addParameter('ylim',[])
ip.parse(varargin{:})
opt = ip.Results;

if isempty(opt.ylim)
    ylim = get(gca,'YLim');
end

% Using patch() for speed of computing
% to see patch documentation, run
% `doc patch`

% Generate vertices
Y = repmat([ylim(1) ylim(2) ylim(2) ylim(1)], 1, size(windowList,1));
Y = Y';
X = repelem(windowList, 1, 2);
X = X';
X = X(:);
vertices = [X, Y]; % list of xy points used to draw windows

% Link vertices into faces (windows)
faces = 1:size(vertices,1);
faces = reshape(faces, 4, [])'; % each row specifies which vertices are linked into a face object (dots that are connected to make a shape)

matlabGraphicsObjects = patch('faces', faces, 'vertices', vertices);
