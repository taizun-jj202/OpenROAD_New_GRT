

# cd /home/tsjafri/TAIZUN/OpenROAD/test
# ################### Commands to run GR on AES design.
# /usr/bin/time -f 'WALL=%E USER=%U SYS=%S' \
#   ../build/bin/openroad_sproute -exit fastroute_aes45.tcl   | tee /home/tsjafri/TAIZUN/TESTS/runlogs_aes_fast.log
# /usr/bin/time -f 'WALL=%E USER=%U SYS=%S' \
#   ../build/bin/openroad_sproute -exit cugr_aes45.tcl        | tee /home/tsjafri/TAIZUN/TESTS/runlogs_aes_cugr.log
# /usr/bin/time -f 'WALL=%E USER=%U SYS=%S' \
#   ../build/bin/openroad_sproute -exit sproute_aes45.tcl     | tee /home/tsjafri/TAIZUN/TESTS/runlogs_aes_sproute.log



# ################### Commands to run GR on JPEG design.
# /usr/bin/time -f 'WALL=%E USER=%U SYS=%S' \
#   ../build/bin/openroad_sproute -exit fastroute_jpegsky130hd.tcl   | tee /home/tsjafri/TAIZUN/TESTS/runlogs_jpeg_fast.log
# /usr/bin/time -f 'WALL=%E USER=%U SYS=%S' \
#   ../build/bin/openroad_sproute -exit cugr_jpegsky130hd.tcl        | tee /home/tsjafri/TAIZUN/TESTS/runlogs_jpeg_cugr.log
# /usr/bin/time -f 'WALL=%E USER=%U SYS=%S' \
#   ../build/bin/openroad_sproute -exit sproute_jpegsky130hd.tcl     | tee /home/tsjafri/TAIZUN/TESTS/runlogs_jpeg_sproute.log



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

cd /home/tsjafri/TAIZUN/OpenROAD/test
/usr/bin/time -f 'WALL=%E USER=%U SYS=%S' \
  ../build/bin/openroad_spR_int -no_init -exit fastroute_aes45.tcl   | tee /home/tsjafri/TAIZUN/TESTS/AES/runlogs_aes_fast.log
# # /usr/bin/time -f 'WALL=%E USER=%U SYS=%S' \
# #   ../build/bin/openroad_spR_int -no_init -exit cugr_aes45.tcl        | tee /home/tsjafri/TAIZUN/TESTS/AES/runlogs_aes_cugr.log
# /usr/bin/time -f 'WALL=%E USER=%U SYS=%S' \
#   ../build/bin/openroad_spR_int -no_init -exit sproute_aes45.tcl     | tee /home/tsjafri/TAIZUN/TESTS/AES/runlogs_aes_sproute.log



# ################### Commands to run GR on JPEG design.
# /usr/bin/time -f 'WALL=%E USER=%U SYS=%S' \
#   ../build/bin/openroad_spR_int -no_init -exit fastroute_jpegsky130hd.tcl   | tee /home/tsjafri/TAIZUN/TESTS/JPEG/runlogs_jpeg_fast.log
# /usr/bin/time -f 'WALL=%E USER=%U SYS=%S' \
#   ../build/bin/openroad_spR_int -no_init -exit cugr_jpegsky130hd.tcl        | tee /home/tsjafri/TAIZUN/TESTS/JPEG/runlogs_jpeg_cugr.log
# /usr/bin/time -f 'WALL=%E USER=%U SYS=%S' \
#   ../build/bin/openroad_spR_int -no_init -exit sproute_jpegsky130hd.tcl     | tee /home/tsjafri/TAIZUN/TESTS/JPEG/runlogs_jpeg_sproute.log
