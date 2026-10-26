%Amazon Code with Oxbows

E = AmazonRiverTest.Easting;
N = AmazonRiverTest.Northing;

E = E(:);
N = N(:);

E_original = E;
N_original = N;
num_points = length(E);




% General Parameters

K = 1000;
dt = 1;
num_iterations = 500;

migration_gain = 10;
max_move_fraction = 0.03;

curvature_window = 7;
max_radius = 15000;         

%Point insertion parameters
enable_point_insertion = true;

insertion_threshold = 3;     
diagnostic_threshold = 2;    

insertion_interval = 10;
max_insertions_per_check = 10;
max_points = 1000;

%Oxbow Parameters
enable_cutoffs = true;

cutoff_distance_factor = 2;
min_index_separation = 10;
min_along_channel_factor = 10;

max_cutoffs = 10;
min_active_points_after_cutoff = 20;


connector_samples = 8;       


initial_segment_lengths = hypot(diff(E), diff(N));
%Starting data point spacing
target_spacing = median(initial_segment_lengths);



max_move = max_move_fraction * target_spacing;

curvature_distance = curvature_window * target_spacing;

insertion_distance = insertion_threshold * target_spacing;
cutoff_distance = cutoff_distance_factor * target_spacing;
min_along_channel_distance = min_along_channel_factor * target_spacing;




snapshot_iterations = [10, 50, 100, 200, 300, 500];

snapshots_E = cell(size(snapshot_iterations));
snapshots_N = cell(size(snapshot_iterations));

oxbow_E = {};
oxbow_N = {};

cutoff_log = zeros(0, 4);

total_cutoffs = 0;
%Iterating over time

for iteration = 1:num_iterations

    num_points = length(E);

   
    segment_lengths = hypot(diff(E), diff(N));
    S = [0; cumsum(segment_lengths)];

    

   
    E_new = E;
    N_new = N;
%curvature
    for k = 2:num_points-1

        s_back = S(k) - curvature_distance;
        s_forward = S(k) + curvature_distance;

        
        if s_back < S(1) || s_forward > S(end)
            continue;
        end

        
        x1 = interp1(S, E, s_back, 'linear');
        y1 = interp1(S, N, s_back, 'linear');

        x2 = E(k);
        y2 = N(k);

        x3 = interp1(S, E, s_forward, 'linear');
        y3 = interp1(S, N, s_forward, 'linear');

        
        denominator = 2 * ( ...
            x1 * (y2 - y3) + ...
            x2 * (y3 - y1) + ...
            x3 * (y1 - y2));

        if abs(denominator) < 1e-10
            continue;
        end

        center_x = ( ...
            (x1^2 + y1^2) * (y2 - y3) + ...
            (x2^2 + y2^2) * (y3 - y1) + ...
            (x3^2 + y3^2) * (y1 - y2)) / denominator;

        center_y = ( ...
            (x1^2 + y1^2) * (x3 - x2) + ...
            (x2^2 + y2^2) * (x1 - x3) + ...
            (x3^2 + y3^2) * (x2 - x1)) / denominator;

        radius = hypot(x2 - center_x, y2 - center_y);

        if ~isfinite(radius) || radius <= 0 || radius > max_radius
            continue;
        end

        %Movement direction
        normal_x = (x2 - center_x) / radius;
        normal_y = (y2 - center_y) / radius;

        %Movement amount
        move_amount = migration_gain * (K / radius) * dt;
        move_amount = min(move_amount, max_move);

        E_new(k) = E(k) + move_amount * normal_x;
        N_new(k) = N(k) + move_amount * normal_y;
    end

   
    E = E_new;
    N = N_new;

    
    E(1) = E_original(1);
    N(1) = N_original(1);
    E(end) = E_original(end);
    N(end) = N_original(end);

    %Point insertion
    if enable_point_insertion && ...
            mod(iteration, insertion_interval) == 0 && ...
            length(E) < max_points

        segment_lengths = hypot(diff(E), diff(N));

        long_segments = find(segment_lengths > insertion_distance);

        if ~isempty(long_segments)

            
            [~, order] = sort( ...
                segment_lengths(long_segments), 'descend');

            number_to_insert = min([ ...
                length(order), ...
                max_insertions_per_check, ...
                max_points - length(E)]);

            selected_segments = sort( ...
                long_segments(order(1:number_to_insert)));

            
            E_inserted = zeros(length(E) + number_to_insert, 1);
            N_inserted = zeros(length(N) + number_to_insert, 1);

            output_index = 1;

            for k = 1:length(E)-1

                E_inserted(output_index) = E(k);
                N_inserted(output_index) = N(k);
                output_index = output_index + 1;

                if any(selected_segments == k)
                    E_inserted(output_index) = ...
                        (E(k) + E(k+1)) / 2;

                    N_inserted(output_index) = ...
                        (N(k) + N(k+1)) / 2;

                    output_index = output_index + 1;
                end
            end

            E_inserted(output_index) = E(end);
            N_inserted(output_index) = N(end);

            E = E_inserted;
            N = N_inserted;
        end
    end

    %Oxbow
    if enable_cutoffs && total_cutoffs < max_cutoffs && ...
            length(E) > min_active_points_after_cutoff

        num_points = length(E);

        segment_lengths = hypot(diff(E), diff(N));
        S = [0; cumsum(segment_lengths)];

        
        distance_matrix = hypot( ...
            E - E.', N - N.');

        index_vector = (1:num_points).';

        index_gap = abs(index_vector - index_vector.');
        along_channel_gap = abs(S - S.');

        
        candidate_mask = ...
            distance_matrix <= cutoff_distance & ...
            index_gap >= min_index_separation & ...
            along_channel_gap >= min_along_channel_distance;

        candidate_mask(1:num_points+1:end) = false;

        
        candidate_mask(1, :) = false;
        candidate_mask(:, 1) = false;
        candidate_mask(end, :) = false;
        candidate_mask(:, end) = false;

        
        candidate_mask = triu(candidate_mask, 1);

        [candidate_i, candidate_j] = find(candidate_mask);

        if ~isempty(candidate_i)

            candidate_distances = distance_matrix( ...
                sub2ind(size(distance_matrix), ...
                candidate_i, candidate_j));

            
            [~, order] = sort(candidate_distances, 'ascend');

            for candidate_number = 1:length(order)

                q = order(candidate_number);

                i = candidate_i(q);
                j = candidate_j(q);
                neck_distance = candidate_distances(q);

                
                remaining_points = num_points - (j - i - 1);

                if remaining_points < min_active_points_after_cutoff
                    continue;
                end

                
                oxbow_E{end+1} = E(i:j);
                oxbow_N{end+1} = N(i:j);

                
                A = [E(i), N(i)];
                B = [E(j), N(j)];

                
                tangent_A = [E(i) - E(i-1), ...
                             N(i) - N(i-1)];

                
                tangent_B = [E(j+1) - E(j), ...
                             N(j+1) - N(j)];

                norm_A = norm(tangent_A);
                norm_B = norm(tangent_B);

                if norm_A <= 0 || norm_B <= 0
                    
                    oxbow_E(end) = [];
                    oxbow_N(end) = [];
                    continue;
                end

                tangent_A = tangent_A / norm_A;
                tangent_B = tangent_B / norm_B;

                
                handle_length = min( ...
                    0.25 * neck_distance, target_spacing);

                %Smoothing connector function
                P0 = A;
                P1 = A + handle_length * tangent_A;
                P2 = B - handle_length * tangent_B;
                P3 = B;

                
                t_values = linspace( ...
                    0, 1, connector_samples + 2);

                t_values = t_values(2:end-1);

                connector = zeros(length(t_values), 2);

                for m = 1:length(t_values)
                    t = t_values(m);

                    connector(m, :) = ...
                        (1-t)^3 * P0 + ...
                        3 * (1-t)^2 * t * P1 + ...
                        3 * (1-t) * t^2 * P2 + ...
                        t^3 * P3;
                end

                
                E = [E(1:i); connector(:,1); E(j:end)];
                N = [N(1:i); connector(:,2); N(j:end)];

                total_cutoffs = total_cutoffs + 1;

                cutoff_log(end+1, :) = [ ...
                    iteration, i, j, neck_distance];

                

                
                break;
            end
        end
    end

    

    
    snapshot_index = find( ...
        snapshot_iterations == iteration, 1);

    if ~isempty(snapshot_index)
        snapshots_E{snapshot_index} = E;
        snapshots_N{snapshot_index} = N;
    end
end
%Plots

figure;
hold on;

plot(E_original, N_original, 'k--', ...
    'LineWidth', 1.2, 'DisplayName', 'Original river');

snapshot_colors = lines(length(snapshot_iterations));

for k = 1:length(snapshot_iterations)

    if ~isempty(snapshots_E{k})
        plot(snapshots_E{k}, snapshots_N{k}, ...
            'Color', snapshot_colors(k,:), ...
            'LineWidth', 1.0, ...
            'DisplayName', ...
            sprintf('Iteration %d', snapshot_iterations(k)));
    end
end

% Show abandoned oxbow loops
for k = 1:length(oxbow_E)
    plot(oxbow_E{k}, oxbow_N{k}, ...
        'm-', 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
end


plot(E, N, 'b-', 'LineWidth', 2.0, ...
    'DisplayName', 'Final active river');

axis equal;
grid on;
xlabel('Easting (m)');
ylabel('Northing (m)');
title('River Evolution with Oxbow Cutoffs');
legend('Location', 'bestoutside');
hold off;
