set router_label "CUGR"
set extra_router_args {-use_cugr}

rename ::global_route ::orig_global_route
proc ::global_route {args} {
  set t0 [clock milliseconds]
  global extra_router_args router_label 
  set final_args [concat $args $extra_router_args]
  set dt [expr {[clock milliseconds] - $t0}]

  exit
}

source "ibex_sky130hd.tcl"