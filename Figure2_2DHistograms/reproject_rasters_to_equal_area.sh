#NOTE: needs to be run in the directory containing the net climate impact, refor albedo radiative forcing, and rasterized biome GeoTIFFs
for rast in `ls *tif`
do
   gdalwarp -t_srs ESRI:54012 -srcnodata 255 -dstnodata 255 $rast "${rast/%.tif/_ECKERTIV.tif}"
done
