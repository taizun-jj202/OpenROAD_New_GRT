set router_label "CUGR"
set extra_router_args {-use_cugr}
# set post_gr_db   "/home/tsjafri/TAIZUN/TESTS/aes_cugr_after_gr.db"
# set post_dr_db   "/home/tsjafri/TAIZUN/TESTS/AES/aes_cugr_after_dr.db"
# set dr_drc_rpt   "/home/tsjafri/TAIZUN/TESTS/AES/aes_cugr_drc.rpt"
# set dr_stats_rpt "/home/tsjafri/TAIZUN/TESTS/AES/aes_cugr_dr_stats.rpt"

rename ::global_route ::orig_global_route
proc ::global_route {args} {
  set t0 [clock milliseconds]
  global extra_router_args router_label 
  set final_args [concat $args $extra_router_args]
  # uplevel 1 [list ::orig_global_route] $final_args
  set dt [expr {[clock milliseconds] - $t0}]
  # puts "$router_label global_route took $dt ms"
  # if {$post_gr_db ne ""} {
  #   write_db $post_gr_db
  # }
  # run detailed routing for reporting purposes; GR timing already captured above
  # if {[info commands detailed_route] ne ""} {
  #   detailed_route \
  #     -output_drc $dr_drc_rpt \
  #     -output_guide "" \
  #     -output_maze $dr_stats_rpt
  # }
  # if {$post_dr_db ne ""} {
  #   write_db $post_dr_db
  # }
  exit
}

source "aes_sky130hd.tcl"
