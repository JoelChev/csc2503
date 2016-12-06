% Joel Cheverie
% 1002924393
function [x0, r, w, maxW] = fitCircleRobust(pts, initx0, initr, normals, sigmaGM)
%
% function [x0, r, w, maxW] = fitCircleRobust(pts, initx0, initr,
%                                  normals, sigmaGM)
%
%  minimize sum_i  rho( a'x_i + b + x_i'x_i, sigma)
%  w.r.t. a and b, where a = -2x_0 and b = x_0'x_0 - r^2
% Input:
%  pts: an Nx2 matrix containing 2D x_i, i=1...N
%  initx0: 2-vector initial guess for centre
%  initr: initial guess for radius 
%  normals: an Nx2 matrix of corresponding 2D edge normal directions nrml_i  
% Output
%  x0 - 2x1 fitted circle's center position
%  r  - fitted circle's radius
%  w  - N x 1 robust weights for edgels.
%  maxW  - maximum possible robust weight (i.e. for zero error)
%          This may be used by the calling code to scale the weights 
%          to be between 0 and 1.
  


x0 = initx0(:)'; 
[~, K] = size(pts');
a = -2*x0;
b = x0 * x0' - initr^2;
previousQElementSum = 0;
for n=1:50
    X = ones(K, 3);
    Y = zeros(K, 1);

    for i=1:K
        X(i, 1:2) = pts(i,:);
        Y(i) = dot(pts(i,:),pts(i,:));
        error = dot(a, pts(i,:)) + b + Y(i);
        w(i) = (sigmaGM/(sigmaGM^2 + error^2))^2; 
        
        
    end

    % To solve for Q, we need to convert the array into a diagonal matrix.
    P = diag(w);
    
    Q = (X' * P * X)^-1 * (-X'*P*Y);
    a(1) = Q(1);
    a(2) = Q(2);
    a = Q(1:2);
    b = Q(3);
    
    % convergence will be defined to be a difference smaller than 0.05 
    % in the summation of the elements of Q. 
    qElementSum = sum(Q);
    if( (n > 1) && (abs(qElementSum - previousQElementSum) < 0.05)) 
       break;
    end
    previousQElementSum = qElementSum;   
end
x0 = a./(-2); 
r = sqrt(dot(x0,x0')-b);
maxW = max(w);

