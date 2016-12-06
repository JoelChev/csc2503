%Joel Cheverie
%1002924393

function [H, Sa] = linEstH(left, right, NUM_RESCALE)
  % [H, Sa] = linEstH(left, right, NUM_RESCALE)
  % Estimate the homography matrix H from two 3 x n matrices of
  % corresponding points left and right.  
  % Here left(:,k) x (H * right(:,k)) apprx= 0 where x denotes
  % the 3D cross product.
  % NUM_RESCALE (default TRUE) uses Hartley's rescaling. Always use
  % rescaling, unless you wish to show how badly the un-normalized
  % algorithm works.
  % Returns H along with the singular values Sa of the 2nPts x 9 homogeneous
  % linear system for H.
    
  if nargin < 3
    NUM_RESCALE = 1;
  end
  
  nPts = size(left,2);
  if nPts < 4 | nPts ~= size(right,2)
    fprintf(2,'lineEstH: Innappropriate number of left and right points.');
    H = [];
    return;
  end
  
  if size(left,1) == 2
    left = [left; ones(1, nPts)];
  else % Normalize to pixel coords
    left = left./repmat(left(3,:), 3,1);
  end
  if size(right,1) == 2
    right = [right; ones(1, nPts)];
  else % Normalize to pixel coords
    right = right./repmat(right(3,:), 3,1);
  end
  
  imPts = cat(3, left, right);
  
  %% Rescale image data for numerical stability.
  if NUM_RESCALE
    Knum = repmat(eye(3), [1,1,2]);
    %%% Rescale for numerical stability
    mn = sum(imPts(1:2,:,:),2)/nPts;
    mns = reshape(mn, [2 1 2]);
    var = sum(sum((imPts(1:2,:,:)-repmat(mns, [1 nPts 1])).^2,2)/nPts, 1);
    %% Scale image points so that sum of variances of x and y = 2.
    scl = sqrt(2./var(:));
    %% Sanity: varScl =  var .* reshape(scl.^2, [1 1 2]) % Should be 2
    
    %% Scale so x and y variance is roughly 1, translate so image mean (x,y) is zero.
    Knum(1:2,3,:) = -mn;
    Knum(1:2,:,:) = Knum(1:2,:,:).*repmat(reshape(scl, [1 1 2]), [2, 3,1]);
    for kIm = 1:2
      imPts(:,:,kIm) = reshape(Knum(:,:,kIm),3,3) * imPts(:,:,kIm);
    end
    %% Sanity check
    % sum(imPts(1:2,:,:),2)/nPts  % Should be [0 0]'
    % sum(sum(imPts(1:2,:,:).^2,2)/nPts,1) % Should be 2.
  end

  %% Make constraint matrix A. 
  %% You fill in below here.
  
    A = zeros(2*nPts,9);
    
    O = [0 0 0];
    for n = 1:nPts
        X = imPts(:,n,2)';
        x = imPts(1,n,1); y = imPts(2,n,1); w = imPts(3,n,1);
        A(2*n-1,:) = [  O  -w*X  y*X];
        A(2*n,:) = [ w*X   O  -x*X];
        %A(3*n  ,:) = [-y*X  x*X   O ];
    end
    
    [U,Sa,V] = svd(A,0); % 'Economy' decomposition for speed
    Sa = diag(Sa);
    % Extract homography
    H_0 = reshape(V(:,end),3,3)';
    
    
    % Denormalise
    T2 = reshape(Knum(:,:,2),3,3);
    T1 = reshape(Knum(:,:,1),3,3);
    H = inv(T1)*H_0*T2;
    
    %Ensure it has unit length.
    H = reshape(H, 9,1);
    H = H/norm(H);
    H = reshape(H,3,3);
      

  