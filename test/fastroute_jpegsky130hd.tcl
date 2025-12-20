set router_label "FastRoute"
set extra_router_args {}
# set gr_log "/home/tsjafri/TAIZUN/TESTS/ORFS+OpenROAD/OR_jpeg/GR_fastroute.log"
# set dr_log "/home/tsjafri/TAIZUN/TESTS/OpenROAD/jpeg_sky130hd/DR_fastroute.log"

# file mkdir [file dirname $gr_log]
# file mkdir [file dirname $dr_log]

rename ::global_route ::orig_global_route
proc ::global_route {args} {
  # set t0 [clock milliseconds]
  global extra_router_args router_label 
  set final_args $args
  if {[llength $extra_router_args] > 0} {
    set final_args [concat $final_args $extra_router_args]
  }
  # set gr_chan [open $gr_log w]
  # puts $gr_chan "[$router_label] global_route args: $final_args"
  # flush $gr_chan
  # redirect $gr_chan {
  #   uplevel 1 [list ::orig_global_route] $final_args
  # }
  # set dt [expr {[clock milliseconds] - $t0}]
  # puts $gr_chan "$router_label global_route took $dt ms"
  # close $gr_chan
  # if {$post_gr_db ne ""} {
  #   write_db $post_gr_db
  # }
  # if {[info commands detailed_route] ne ""} {
  #   set dr_chan [open $dr_log w]
  #   puts $dr_chan "[$router_label] detailed_route start"
  #   flush $dr_chan
  #   redirect $dr_chan {
  #     detailed_route \
  #       -output_drc $dr_drc_rpt \
  #       -output_guide "" \
  #       -output_maze $dr_stats_rpt
  #   }
  #   puts $dr_chan "[$router_label] detailed_route finished"
  #   close $dr_chan
  # }
  # if {$post_dr_db ne ""} {
  #   write_db $post_dr_db
  # }
  exit
}

# source "jpeg_sky130hd.tcl"
source "jpeg_sky130hd.tcl"
