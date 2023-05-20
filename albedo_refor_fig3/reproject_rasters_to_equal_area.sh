for rast in `ls *tif`
do
   gdalwarp -t_srs ESRI:54012 -srcnodata 255 -dstnodata 255 $rast "${rast/%.tif/_ECKERTIV.tif}"
done
