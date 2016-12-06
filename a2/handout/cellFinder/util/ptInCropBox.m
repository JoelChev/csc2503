function b = ptInCropBox(pt, cropBox)
  b = (pt(1) >= cropBox(1,1) && pt(1) <= cropBox(2,1));
  b = b & (pt(2) >= cropBox(1,2) && pt(2) <= cropBox(2,2));
  return;
  