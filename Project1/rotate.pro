PRO rotate, infn
  cube = readfits(infn)
  rotnth = readfits('rotnth.fits')
  rotated = fltarr((size(cube))[1:(size(cube))[0]])
  FOR i=0, (size(cube))[(size(cube))[0]] - 1 DO rotated[*, *, i] = $
     rot(cube[*,*,i],rotnth[i],1.,1224./2.,1224./2.,/piv,cub=-0.5)
  writefits, 'rotated.fits', rotated
  median_im = median(rotated, dimension = 3)
  writefits, 'median.fits', median_im
END 








