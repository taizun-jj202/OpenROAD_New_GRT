cd /home/tsjafri/TAIZUN/OpenROAD/test

/usr/bin/time -f 'WALL=%E USER=%U SYS=%S' /home/tsjafri/TAIZUN/OpenROAD/build/bin/openroad_spR_int -no_init -exit fastroute_jpegsky130hd.tcl | tee /home/tsjafri/TAIZUN/TESTS/OpenROAD/jpeg_sky130hd/GR_fastroute.log
/usr/bin/time -f 'WALL=%E USER=%U SYS=%S' /home/tsjafri/TAIZUN/OpenROAD/build/bin/openroad_spR_int -no_init -exit cugr_jpegsky130hd.tcl | tee /home/tsjafri/TAIZUN/TESTS/OpenROAD/jpeg_sky130hd/GR_cugr.log
/usr/bin/time -f 'WALL=%E USER=%U SYS=%S' /home/tsjafri/TAIZUN/OpenROAD/build/bin/openroad_spR_int -no_init -exit sproute_jpegsky130hd.tcl | tee /home/tsjafri/TAIZUN/TESTS/OpenROAD/jpeg_sky130hd/GR_sproute.log
