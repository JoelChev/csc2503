% Joel Cheverie
% 1002924393
function [goodCircle] = isGoodCircle(x0, r, w, ...
                                     circleEstimates, nFound)
  % [goodCircle] = isGoodCircle(x0, r, w, normals, ...
  %                                  circleEstimates, nFound)
  % Decide if the circle with parameters x0 and r, with fitted
  % weights w, is to be added to the current set of circles.
  % Input:
  %  x0 2-vector for circle center
  %  r  radius of circle
  %  w  robust weights of edgels supporting this circle,
  %     weights have been scaled to have max possible value 1.
  %  circleEstimates C x 3 array of circle parameters, the first
  %     nFound rows of which specifies the [x0(1), x0(2) r] parameters
  %     for a previously fitted circle for this data set.
  %  nFound the number of previously fitted circles stored in circleEstimates
  % Output:
  %  goodCircle boolean, true if you wish to add this circle to the
  %             list of previously fitted circles.
  
  x0 = x0(:)';  
  weightThreshold = 0.25;
  angleThreshold = 3*pi;
  goodCircle = true;
  % compare with weight threshold for first circle
  if(nFound == 0)
      if (sum(w) <= 2 * pi * r * weightThreshold)
          goodCircle = false;
      end
      return;
  end
  
 angleSum=0; 
 for i=1:nFound
      radius_1 = circleEstimates(3,i);
      x_1 = circleEstimates(1:2,i);
      distance = norm(x0-x_1);
      
      %This implies the circles are overlapping, so we cannot calculate the
      %intersection angle to make a decision. Just default to it being bad.
      if distance <= (r-radius_1)
          goodCircle = false;
          return;
      end
      %If circles are intersecting, then we calculate the total angle of
      %intersection for every circle that has been found so far. If this
      %value is more than the angleThreshold, reject the circle.
      if distance < (r+radius_1)^2
        d1 =(distance^2 - r^2 + radius_1^2 )/(2*distance);
        angle = 2 * asin(abs(d1/r));
        angleSum = angleSum + angle;
        if angleSum > angleThreshold
            goodCircle = false;
            return;
        end
      end
  end
