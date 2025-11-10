#ifndef OUT_OF_GUIDE_H
#define OUT_OF_GUIDE_H

void write_guide_points() {


    ofstream outfile("guidepoints.txt");

    for(int netID = 0; netID < defDB.nets.size(); netID++) {
        auto& net = defDB.nets[netID];
        outfile << net.name << " " << net.guide_points.size() << endl;
        for(auto p : net.guide_points) {
            outfile << p.x << " " << p.y << " " << p.z <<  " " << p.horizontal << endl;
        }
    }
    outfile.close();
}

void read_guide_points() {
    ifstream infile("guidepoints.txt");

    for(int netID = 0; netID < defDB.nets.size(); netID++) {
        auto& net = defDB.nets[netID];
        int npoints;
        string net_name;
        infile >> net_name>> npoints;
        if(net.name != net_name) {
            cout << " net name wrong: " << net.name << " " << net_name << endl;
        }
        for(int i = 0; i < npoints; i++) {
            parser::Routed3DPoint p;
            infile >> p.x >> p.y >> p.z >> p.horizontal;
            net.guide_points.insert(p);
        }
    }
    cout << "read guide points done: " << defDB.nets[0].guide_points.size() << endl;
    infile.close();
}

bool find_nodes_in_routed_points(std::set<parser::Routed3DPoint>& guide_points, int x, int y, int z, bool is_horizontal) {

    if(guide_points.find(parser::Routed3DPoint(x, y, z, is_horizontal)) != guide_points.end()) {
        return true;
    }
    else if(is_horizontal && y - 1 >= 0 && guide_points.find(parser::Routed3DPoint(x, y - 1, z, is_horizontal)) != guide_points.end()) {
        return true;
    }
    else if(!is_horizontal && x - 1 >= 0 && guide_points.find(parser::Routed3DPoint(x - 1, y, z, is_horizontal)) != guide_points.end()) {
        return true;
    }
    else 
        return false;
}

void plot_out_of_guide() {
    bool print = false;

    float* oog = new float [xGrid * yGrid];
    for(int i = 0; i < xGrid * yGrid; i++) 
        oog[i] = 0;

    /*for(auto y_boundary : defDB.yGcellBoundaries) 
        cout << y_boundary << " ";
    cout << endl;*/


    bool* is_horizontal = new bool [numLayers];
	bool m1_horizontal = (lefDB.layers[0].direction == "HORIZONTAL")? true : false; //by default m1 is vertical
	for(int i = 0; i < numLayers; i++) {
		if(i % 2 == 0)
			is_horizontal[i] = m1_horizontal;
		else	
			is_horizontal[i] = !m1_horizontal;
	}
    //find_Gcell
    int hist_size = 100;
    int* oog_hist = new int [hist_size];
    int max_oog = 0;
    int max_oog_netID = 0;
    for(int i = 0; i < hist_size; i++) {
        oog_hist[i] = 0;
    } 

    for(int netID = 0; netID < defDB.nets.size(); netID++) {
        if(netID % 10000 == 0)
            cout << " nets: " << netID << endl;
        auto& net = defDB.nets[netID];
        
        if(netID == 27148)
            print = true;
        /*for(auto p : net.guide_points) {
            if(p.z == 7)
                cout << p.x << " " << p.y << " " << p.z << " " << p.horizontal << endl;
        }*/
        int oog_cnt = 0;
        if(net.paths.size() == 0)
            continue;
        
        for(auto path : net.paths) {
            bool has_end = false;
            if(path.end.x != 0 || path.end.y != 0)
                has_end = true;

            int start_x, start_y, end_x = -1, end_y = -1;
            int layerID = lefDB.layer2idx[path.layerName];
            int z = layerID / 2;
            start_x = parser::find_Gcell(path.begin.x, defDB.xGcellBoundaries);
            start_y = parser::find_Gcell(path.begin.y, defDB.yGcellBoundaries);
            if(has_end) {
                end_x = parser::find_Gcell(path.end.x, defDB.xGcellBoundaries);
                end_y = parser::find_Gcell(path.end.y, defDB.yGcellBoundaries);
            }


            if(!has_end) {
                int x = start_x;
                int y = start_y;
                if(!find_nodes_in_routed_points(net.guide_points, x, y, z, is_horizontal[z])) {
                    oog[y * xGrid + x]++;
                    oog_cnt++;
                    if(print) {
                        cout << "===" << net.name;
                        cout << " oog : " << x << ", " << y << ", " << z << endl; 
                        path.print();
                    }

                }
            }
            else if(start_x == end_x) {
                int x = start_x;
                int big_y = max(start_y, end_y);
                int small_y = min(start_y, end_y);
                for(int y = small_y; y <= big_y; y++) {
                    if(!find_nodes_in_routed_points(net.guide_points, x, y, z, is_horizontal[z])) {
                        oog[y * xGrid + x]++;
                        oog_cnt++;
                        if(print) {
                            cout << "===" << net.name;
                            cout << " oog : " << x << ", " << y << ", " << z << endl; 
                            path.print();
                        }
                    }
                }
            }
            else if(start_y == end_y) {
                int y = start_y;
                int big_x = max(start_x, end_x);
                int small_x = min(start_x, end_x);
                for(int x = small_x; x <= big_x; x++) {
                    if(!find_nodes_in_routed_points(net.guide_points, x, y, z, is_horizontal[z])) {
                        oog[y * xGrid + x]++;
                        oog_cnt++;
                        if(print) {
                            cout << "===" << net.name;
                            cout << " oog : " << x << ", " << y << ", " << z << endl; 
                            path.print();
                        }
                    }
                }
            }
            else {
                cout << "wire is not V or H? " << endl;
                path.print();
                exit(1);
            }
        }
        if(oog_cnt > max_oog) {
            max_oog = oog_cnt;
            max_oog_netID = netID;
        }
        if(oog_cnt >= 100 && netID != 1004) {
            cout << "first 100 net: " << netID << " " << defDB.nets[netID].name << endl;
            break;
        }

        if(oog_cnt < hist_size)
            oog_hist[oog_cnt]++;
        else   
            oog_hist[hist_size - 1]++;


    }
    int total_oog = 0;
    for(int i = 0; i < hist_size; i++) {
        cout << i << ": " << oog_hist[i] << endl;
        total_oog += oog_hist[i];
    }
    cout << "total oog: " << total_oog << endl;
    cout << "max_oog net: " << max_oog_netID << " " << defDB.nets[max_oog_netID].name << " " << max_oog << endl;

    ofstream outfile("heat.txt");
    cout << "oog? " << endl;
    outfile << "oog" << endl;
    outfile << "(" << endl;
    
    for(int x = 0; x < xGrid; x++) {
        for(int y = 0; y < yGrid; y++) {
            int llx  = defDB.xGcellBoundaries.at(x);		
            int lly = defDB.yGcellBoundaries.at(y);
            int urx = defDB.xGcellBoundaries.at(x + 1);
            int ury  = defDB.yGcellBoundaries.at(y + 1);
            outfile << llx << " ";
            outfile << lly << " ";
            outfile << urx << " ";
            outfile << ury << " ";
            outfile << oog[y * xGrid + x] << endl;
        }
    }
    outfile << ")" << endl;
    outfile.close();

    delete[] is_horizontal;
    delete[] oog_hist;
    delete[] oog;
}

#endif