% Joel Cheverie
% 1002924393
function [depth] = getDepthFromNormals(n, mask)
  % [depth] = getDepthFromNormals(n, mask)
  %
  % Input:
  %    n is an [N, M, 3] matrix of surface normals (or zeros
  %      for no available normal).
  %    mask logical [N,M] matrix which is true for pixels
  %      at which the object is present.
  % Output
  %    depth an [N,M] matrix providing depths which are
  %          orthogonal to the normals n (in the least
  %          squares sense).
  %
  [M, N] = size(mask);
  
  norm_x = squeeze(n(:,:,1));
  norm_y = squeeze(n(:,:,2));
  norm_z = squeeze(n(:,:,3));
  A = sparse(2*M*N, M*N);
  v = zeros(2*M*N, 1);
  for n_index = 1:N
      for m_index = 1:M
          %Only consider pixels where
          %the object exists.
          if mask(m_index,n_index)~=0
              % Figure out which row and column
              % is currently being considered.
              % First for y component, then x component.
              y_offset = M*N;
              row_y = (n_index-1)*M+m_index + y_offset;
              column = row_y - y_offset;
              %Get the individual norm components
              %for the current indices.
              norm_x_mn = norm_x(m_index,n_index);
              norm_y_mn = norm_y(m_index,n_index);
              norm_z_mn = norm_z(m_index,n_index);
              
              % Due to coordinate system we have:
              % -norm_y = norm_z(Z_x,y+1 -Z_x,y)
              A(row_y,column) = -norm_z_mn;
              A(row_y,column+1) = norm_z_mn;
              v(row_y) = -norm_y_mn;
              
              row_x = row_y - y_offset;
              % Likewise:
              % -norm_x = norm_z(Z_x+1,y -Z_x,y)
              A(row_x,column) = -norm_z_mn;
              A(row_x,column+M) = norm_z_mn;
              v(row_x) = -norm_x_mn;
          end
      end
  end
  % Need to turn off the warning otherwise
  % it always pops up.
  warning('off','MATLAB:rankDeficientMatrix');
  z = A\v;
  depth = reshape(z, M, N);
  % Only want to see the depth of the object
  % everything else is irrelevant.
  depth(mask==0) = 0;
