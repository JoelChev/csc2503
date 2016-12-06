% Joel Cheverie
% 1002924393
function [circle] = bestProposal(circles, sigmaGM, normals, p)
% [circle] = bestProposal(circles, sigmaGM, normals, p)
% Chose the best circle from the proposed circles.
% Input
%  circles K x 3 matrix each row containing [x0(1), x0(2), r] for 
%          the center and radius of a proposed circle.
%  sigmaGM - the robust scale parameter 
%  normals - P x 2 edgel normal data
%  p       - P x 2 edgel position data.
% Output
%  circle  1 x 3  best circle params [x0(1), x0(2), r] from the proposals
  
  [dimension, edgelCount] = size(p');
    [circle_dimension, circle_guess_count] = size(circles');
    %initialize an array to accumulate the error. I will be
    %using the minimum value in this array later.
    sum_gm_array = zeros(1,circle_guess_count);
    for i=1:circle_guess_count
        gm_sum = 0;
        %compute e_j for each edgel
        for j = 1:edgelCount
            e_j = (norm(p(j,:)-circles(i,1:2)))^2 - circles(i,3)^2;
            gm = e_j^2/(sigmaGM^2 + e_j^2);
            gm_sum = gm_sum + gm;
        end
        sum_gm_array(1,i) = gm_sum;
    end
    %Use the circle that has minimum GM error over its edgels as the best
    %proposal.
    [minimum_error, index] = min(sum_gm_array);
    circle = circles(index,1:3);
