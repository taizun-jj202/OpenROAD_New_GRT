cd /home/tsjafri/TAIZUN/OpenROAD/test
# /usr/bin/time -f 'WALL=%E USER=%U SYS=%S' ../build/bin/openroad -exit cugr_jpegsky130hd.tcl | tee runlogs_jpeg_cugr.log
# /usr/bin/time -f 'WALL=%E USER=%U SYS=%S' ../build/bin/openroad -exit sproute_jpegsky130hd.tcl | tee runlogs_jpeg_sproute.log
/usr/bin/time -f 'WALL=%E USER=%U SYS=%S' ../build/bin/openroad -exit fastroute_jpegsky130hd.tcl | tee runlogs_jpeg_fastroute.log
