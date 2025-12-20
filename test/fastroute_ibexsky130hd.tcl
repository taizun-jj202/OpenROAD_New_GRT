set router_label "FastRoute"
set extra_router_args {}

rename ::global_route ::orig_global_route
proc ::global_route {args} {
  set t0 [clock milliseconds]
  global extra_router_args router_label 
  set final_args $args
  if {[llength $extra_router_args] > 0} {
    set final_args [concat $final_args $extra_router_args]
  }
  set dt [expr {[clock milliseconds] - $t0}]
  exit
}

source "ibex_sky130hd.tcl"
