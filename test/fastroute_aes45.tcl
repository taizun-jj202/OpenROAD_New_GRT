set router_label "FastRoute"
set extra_router_args {}
set post_gr_db   "/home/tsjafri/TAIZUN/TESTS/AES/aes_fast_after_gr.db"
set post_dr_db   "/home/tsjafri/TAIZUN/TESTS/AES/aes_fast_after_dr.db"
set dr_drc_rpt   "/home/tsjafri/TAIZUN/TESTS/AES/aes_fast_drc.rpt"
set dr_stats_rpt "/home/tsjafri/TAIZUN/TESTS/AES/aes_fast_dr_stats.rpt"


rename ::global_route ::orig_global_route
proc ::global_route {args} {
  global extra_router_args router_label post_gr_db post_dr_db dr_drc_rpt dr_stats_rpt
  set final_args $args
  if {[llength $extra_router_args] > 0} {
    set final_args [concat $final_args $extra_router_args]
  }
  set t0 [clock milliseconds]
  uplevel 1 [list ::orig_global_route] $final_args
  set dt [expr {[clock milliseconds] - $t0}]
  puts "$router_label global_route took $dt ms"
  if {$post_gr_db ne ""} {
    write_db $post_gr_db
  }
  if {[info commands detailed_route] ne ""} {
    detailed_route \
      -output_drc $dr_drc_rpt \
      -output_guide "" \
      -output_maze $dr_stats_rpt
  }
  if {$post_dr_db ne ""} {
    write_db $post_dr_db
  }
  exit
}

source "aes_nangate45.tcl"
