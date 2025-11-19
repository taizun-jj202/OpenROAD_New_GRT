set router_label "SPRoute"
set extra_router_args {-router sproute}
set post_gr_db "/home/tsjafri/TAIZUN/TESTS/aes_sproute_after_gr.db"

rename ::global_route ::orig_global_route
proc ::global_route {args} {
  global extra_router_args router_label post_gr_db
  set final_args [concat $args $extra_router_args]
  set t0 [clock milliseconds]
  uplevel 1 [list ::orig_global_route] $final_args
  set dt [expr {[clock milliseconds] - $t0}]
  puts "$router_label global_route took $dt ms"
  if {$post_gr_db ne ""} {
    write_db $post_gr_db
  }
  exit
}

source "aes_nangate45.tcl"
