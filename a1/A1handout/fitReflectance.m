% Joel Cheverie
% 1002924393
function [n, albedo] = fitReflectance(im, L)
  % [n, albedo] = fitReflectance(im, L)
  % 
  % Input:
  %   im - nPix x nDirChrome array of brightnesses,
  %   L  - 3 x nDirChrome array of light source directions.
  % Output:
  %   n - nPix x 3 array of surface normals, with n(k,1:3) = (nx, ny, nz)
  %       at the k-th pixel.
  %   albedo - nPix x 1 array of estimated albdedos
    

  g = L'\ im';
  albedo = sqrt(sum(g.^2))';
  % We need to normalize g to recover n.
  n = g' ./repmat(albedo, 1, 3);
  return;


