#!/usr/bin/env octave-cli
% Copyright 2019 Sotiris Papatheodorou
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% Visualize the text file written by calling Octree::writeAllocatedNodes().
% This will show all allocated Nodes (and VoxelBlocks) as colored cubes in 3D.
% Visualizing the cubes may take a while, especially in large volumes.

% NOTES
% Filename format: map_2019-08-22_184424_000009.bin

clear variables



% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
world_names     = {'apartment', 'apartment_64', 'apartment_256', 'maze', 'powerplant'};
world_dims_x    = [10  10  10  20         33];
world_dims_y    = [20  20  20  20         31];
world_dims_z    = [ 3   3   3   2.5       26];
world_offsets_x = [ 0   0   0   0.018794   7.5];
world_offsets_y = [ 0   0   0  -5.61      -8.5];
world_offsets_z = [ 0   0   0   0          0];

plot_path   = false;
interactive = false;
export_plot = true;
export_data = true;



% Patterns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
voxel_volume_pattern    = 'Explored voxel volume: +\d+';
node_volume_pattern     = 'Explored node volume: +\d+';
explored_volume_pattern = 'Explored volume: +\d+';
timestamp_pattern       = '\d{4}-\d{2}-\d{2}_\d{6}';



% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function matched = match_pattern(pattern, line)
  s = regexp(line, pattern);
  matched = !isempty(s);

  global DEBUG;
  if DEBUG && matched
    fprintf('M %s\n', upper(inputname(1)));
  end
end



function p = get_pattern(pattern, line)
  [~, ~, ~, m, ~, ~, ~] = regexp(line, pattern);
  p = m{end};
end



function [filenames, world] = parse_arguments()
	% Get the command line arguments.
	args = argv();
	if isempty(args)
      name = program_invocation_name();
	  fprintf('Usage: %s [-w WORLD] FILE1 [FILE2 ...]\n', name);
	  fprintf('  Use bash globbing with * to select all files of interest, e.g.\n');
	  fprintf('  %s map_1/map_2019-08-22_184424_*\n', name);
      filenames = {};
      world = '';
	  return;
	end

	% Get the world name.
    if strcmp(args{1}, '-w')
      world = args{2};
      args = args(3:end);
    else
      world = 'apartment';
    end

	% Sort the filenames.
	filenames = sort(args);
end



% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the path to the program used to evaluate maps.
[script_dir, ~, ~] = fileparts(program_invocation_name());
voxelcounter_program = [script_dir ...
    '/../../build/Release/se_apps/se-denseslam-ofusion-compute-volume'];
addpath(genpath([script_dir '/octave_functions']));

t = [];
total_volume = [];
voxel_volume = [];
node_volume  = [];
poses = {};

% Parse command line arguments.
[filenames, world] = parse_arguments();
if isempty(filenames)
  return;
end
% Find the world index.
world_ind = find(ismember(world_names, world));
if isempty(world_ind)
	world_ind = 1;
end
% Get the world parameters.
dim_x = world_dims_x(world_ind);
dim_y = world_dims_y(world_ind);
dim_z = world_dims_z(world_ind);
off_x = world_offsets_x(world_ind);
off_y = world_offsets_y(world_ind);
off_z = world_offsets_z(world_ind);
world_volume = dim_x * dim_y * dim_z;

% Iterate over each file.
for i = 1:length(filenames);
  filename = filenames{i};



  % This is a cropped map, skip.
  if strfind(filename, '_cropped.bin')
    continue;
  end



  % This is a map, process.
  if strfind(filename, '.bin')
    % Evaluate the file.
    [status, output] = system([voxelcounter_program ' ' filename ' ' ...
        num2str(dim_x) ' ' num2str(dim_y) ' ' num2str(dim_z) ' ' ...
        num2str(off_x) ' ' num2str(off_y) ' ' num2str(off_z)]);
	if status ~= 0
		fprintf('Error from %s\n', voxelcounter_program);
		fprintf('%s', output);
		continue;
	end

    % Parse the output.
    output_lines = strsplit(output, '\n');
    for l = 1:length(output_lines)
      line = output_lines{l};

      if match_pattern(voxel_volume_pattern, line)
        explored_volume = sscanf(line, 'Explored voxel volume: %f', 1);
        voxel_volume = [voxel_volume explored_volume];

      elseif match_pattern(node_volume_pattern, line)
        explored_volume = sscanf(line, 'Explored node volume: %f', 1);
        node_volume = [node_volume explored_volume];

      elseif match_pattern(explored_volume_pattern, line)
        explored_volume = sscanf(line, 'Explored volume: %f', 1);
        total_volume = [total_volume explored_volume];
      end
    end

    timestamp = str2double(filename(end-9:end-4));
    t = [t timestamp];
  end



  % This is the pose list.
  if plot_path && strfind(filename, '.txt')
    data = importdata(filename);
    num_poses = size(data, 1);
    poses = cell([1 num_poses]);
    for i = 1:num_poses
      poses{i} = data(i, 2:17);
      poses{i} = reshape(poses{i}, 4, 4)';
    end
  end
end



% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw = 2;
figure;
hold on;
grid on;

plot(t, total_volume, 'bo-', 'LineWidth', lw);
plot(t, node_volume,  'go-', 'LineWidth', lw);
plot(t, voxel_volume, 'ro-', 'LineWidth', lw);
xlabel('Time (s)');
ylabel('Explored volume (m^3)');
legend('Total volume', 'Node volume', 'Voxel volume', 'Location', 'southeast');
axis([0 15*60], [0 world_volume]);

if export_plot
  directory = fileparts(filenames{1});
  timestamp = get_pattern(timestamp_pattern, filenames{1});
  image_name = [directory '/' 'data_plot_' timestamp '.png'];
  print(image_name);
end

if export_data
  directory = fileparts(filenames{1});
  timestamp = get_pattern(timestamp_pattern, filenames{1});
  data_file_name = [directory '/' 'data_' timestamp '.csv'];
  % The columns of the .csv file are:
  % timestamp, volume of explored voxels, volume of explored nodes, total
  % explored volume
  csvwrite(data_file_name, [t' voxel_volume' node_volume' total_volume']);
end

if plot_path && ~isempty(poses)
  figure;
  hold on;
  axis equal;
  grid on;

  for i = 1:length(poses)
    plot_axes(poses{i});
  end
  for i = 1:length(poses)-1
    pos_c = poses{i}(1:3, 4);
    pos_n = poses{i+1}(1:3, 4);
    plot3([pos_c(1) pos_n(1)], [pos_c(2) pos_n(2)], [pos_c(3) pos_n(3)], ...
        'm.-', 'LineWidth', lw);
  end

  xlabel('x (m)');
  ylabel('y (m)');
  zlabel('z (m)');
end

if interactive
  ginput();
else
  pause(0.01);
end

