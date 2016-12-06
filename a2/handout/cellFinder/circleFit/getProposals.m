% Joel Cheverie
% 1002924393
function [circles] = getProposals(normals, p, numGuesses)
  % [circles] = getProposals(normals, p, numGuesses)
  % Attempt to produce up to numGuesses circle proposals from
  % the edgel data p and normals.  For typical data sets
  % we will be able to produce numGuesses proposals.  However,
  % on some datasets, say with only a few edgels, we may not be
  % able to generate any proposals.  In this case size(circles,1)
  % can be zero.
  % Input:
  %  normals - N x 2 edgel normals
  %  p         N x 2 edgel positions
  %  numGuesses - attempt to propose this number of circles.
  % Return:
  %   circles a P x 3 array, each row contains [x0(1) x0(2) r]
  %           with 0 <= P <= numGuesses.
  
  [dimension edgelCount] = size(p');
  min_x = min(p(:,1));
  max_x = max(p(:,1));
  min_y = min(p(:,2));
  max_y = max(p(:,2));

  %Need to get the largest radius to constrain the estimation
  x_range = max_x - min_x;
  y_range = max_y - min_y;
  radius = sqrt(x_range^2 + y_range^2);
  radius_ratio = 1.4;
  max_r = ceil(radius);
  
  circles = zeros(0, 3);
  
  % Pair each edgel with another edgel and find their intersection.
  % I am assuming each circle has at least two edgels to choose from.
  for n = 1:numGuesses
    while(1>0)
         % Randomly select two edgels to consider for an intersection.
         random_edgel_proposal = randperm(edgelCount, 2);
         point_index_1 = random_edgel_proposal(1);
         point_index_2 = random_edgel_proposal(2);
      
          % Get the slopes.
          slope_1 = normals(point_index_1, 2) / normals(point_index_1, 1);
          slope_2 = normals(point_index_2, 2) / normals(point_index_2, 1);
          % If the lines are parallel, there will be no intersections, so
          % ignore.
          if(slope_1 == slope_2)
              continue
          end
          
          x_1 = p(point_index_1,1);
          y_1 = p(point_index_1,2);
          x_2 = p(point_index_2,1);
          y_2 = p(point_index_2,2);
          
         b_1 = y_1 - slope_1 * x_1;
         b_2 = y_2 - slope_2 * x_2;
         % Now find where the lines intersect. Take that as the proposed
         % x_c.
         x_cross = (b_2 - b_1) / (slope_1 - slope_2);
         y_cross = slope_1 * x_cross + b_1;
         % check if the intersection point lies in the direction of
         % the two normals. If not, it is outside the circle, so continue.
         if(((normals(point_index_1,:) >= 0) == ([x_cross,y_cross] >= [x_1,y_1])) ~= ...
             ((normals(point_index_2,:) >= 0) == ([x_cross, y_cross] >= [x_2,y_2])))
                continue
         end
        % find the two actual radii
        radius_1_euclid_term = [x_cross, y_cross] - [p(point_index_1, 1), p(point_index_1,2)];
        radius_1 = sqrt(dot(radius_1_euclid_term, radius_1_euclid_term'));
        radius_2_euclid_term = [x_cross, y_cross] - [p(point_index_2, 1), p(point_index_2,2)];
        radius_2 = sqrt(dot(radius_2_euclid_term, radius_2_euclid_term'));
        %Don't consider radii that are too large.
        if(radius_1 >= max_r || radius_2 >= max_r)
            continue
        end
        % if the caclulated radii are very different, then the estimated
        % circle is incorrect, need to try again.
        if (radius_1 / radius_2 > radius_ratio || ...
            radius_2 / radius_1 > radius_ratio)
            continue
        end

            r = mean([radius_1, radius_2]);
            circles(n, :) = [x_cross, y_cross, r];
            break;
      end 
  end
end
