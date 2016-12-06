function displayFittedPoints(p, q, cropOffset)
%
% function displayFittedPoints(p, q, cropOffset)
%  
%  plot points p on top of current figure, colored according
%  to integer labels in q (quartile indicators: 0:3).
%  cropOffset maps from form coordinates of points to local window coords.
%

c = 'rmcb';
for n = 1:4
    ptInd = (q == (4-n));
    tmpPts = p(find(ptInd), :);
    plot(tmpPts(:,1)+cropOffset(1), tmpPts(:,2)+cropOffset(2), ...
        [c(n),'o'], 'markerFaceColor', c(n));
end
