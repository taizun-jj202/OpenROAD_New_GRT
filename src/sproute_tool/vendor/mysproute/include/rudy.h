#ifndef RUDY_H
#define RUDY_h

float plot_rudy(float* rudy, Algo algo) {

    cout << xGrid << " " << yGrid << " grid" << endl;
    float max_rudy = -1;
    for(int i = 0; i < xGrid * yGrid; i++)
        rudy[i] = 0;

    for(int netID = 0; netID < numValidNets; netID++) {
        int min_x = xGrid + 1;
        int min_y = yGrid + 1;
        int max_x = 0;
        int max_y = 0;

        int deg = nets[netID]->deg;
        for(int pinID = 0; pinID < deg; pinID++) {
            min_x = min(min_x, (int) nets[netID]->pinX[pinID]);
            max_x = max(max_x, (int) nets[netID]->pinX[pinID]);
            min_y = min(min_y, (int) nets[netID]->pinY[pinID]);
            max_y = max(max_y, (int) nets[netID]->pinY[pinID]);
        }

        if(max_x - min_x < 0) {
            cout << "Error: " << min_x << " " << max_x << endl;
            exit(1);
        }

        if(max_y - min_y < 0) {
            cout << "Error: " << min_y << " " << max_y << endl;
            exit(1);
        }

        if(max_x - min_x <= SMALL_NET_THRSHD && max_y - min_y <= SMALL_NET_THRSHD) {
            nets[netID]->small = false;
        }
        else    {
            nets[netID]->small = true;
        }

        float wl = max_x - min_x + 1 + max_y - min_y + 1;
        float local_rudy = wl / ((float) (max_x - min_x + 1) * (float) (max_y - min_y + 1));
        if(isnan(local_rudy)) {
            cout << "here: " << min_x << " " << max_x << " " << min_y << " " <<max_y << endl;
            exit(1);
        }
        if(local_rudy <= 0) {
            cout << "Error: negative rudy! " << local_rudy << " " << max_x - min_x << " " << max_y - min_y << endl;
            exit(1);
        }
        for(int x = min_x; x <= max_x; x++) {
            for(int y = min_y; y <= max_y; y++) {
                rudy[y * xGrid + x] += local_rudy;
                if(isnan(rudy[y * xGrid + x])) {
                    cout << min_x << " " << max_x << " " << min_y << " " << max_y << endl;
                    exit(1);
                }
            }
        }
    }

    for(int netID = 0; netID < numInvalidNets; netID++) {
        int x = (int) invalid_nets[netID]->pinX[0];
        int y = (int) invalid_nets[netID]->pinY[0];
        rudy[y * xGrid + x] += 1;
    }

    for(int i = 0; i < xGrid * yGrid; i++)
        max_rudy = max(max_rudy, rudy[i]);

    if(algo == RUDY) {
        float rudy_sum = 0;
        int hist_size = 30;
        int hist[hist_size];
        for(int i = 0; i < hist_size; i++){
            hist[i] = 0;
        }
        
        ofstream outfile("heat.txt");
        cout << "heat? " << endl;
        outfile << "heat" << endl;
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
                outfile << rudy[y * xGrid + x] << endl;
                float rudy1 = rudy[y * xGrid + x];
                rudy_sum += rudy1;
                if((int) rudy1 < hist_size)
                    hist[(int) rudy1]++;
                else hist[hist_size - 1]++;
                
            }
        }
        outfile << ")" << endl;
        outfile.close();

        for(int i = 0; i < hist_size; i++) {
            cout  << hist[i] << endl;
        }

        cout << "avg: " << rudy_sum / (float)(xGrid * yGrid) << endl;
    }

    return max_rudy;

}

bool paint_color(float* color, float* rudy, int x_in, int y_in, int color_cnt, float color_threshold) {

    std::queue<int> WL;
    WL.push(y_in * xGrid + x_in);
    std::set<int> pixel_set;
    while(!WL.empty()) {
        int index = WL.front();
        WL.pop();
        int x = index % xGrid;
        int y = index / xGrid;
        if(color[index] != 0)
            continue; 
        //cout << x << " " << y << " " << color_cnt << " rudy: " << rudy[index] << endl;
        color[index] = color_cnt;
        pixel_set.insert(index);

        if(y < yGrid - 1) {
            if(color[(y + 1) * xGrid + x] == 0 && rudy[(y + 1) * xGrid + x] > color_threshold) {
                WL.push((y + 1) * xGrid + x);
            }
        }

        if(x < xGrid - 1) {
            if(color[y * xGrid + (x + 1)] == 0 && rudy[y * xGrid + (x + 1)] > color_threshold) {
                WL.push(y * xGrid + (x + 1));
            }
        }

        if(y > 0) {
            if(color[(y - 1) * xGrid + x] == 0 && rudy[(y - 1) * xGrid + x] > color_threshold) {
                WL.push((y - 1) * xGrid + x);
            }
        }

        if(x > 0) {
            if(color[y * xGrid + (x - 1)] == 0 && rudy[y * xGrid + (x - 1)] > color_threshold) {
                WL.push(y * xGrid + (x - 1));
            }
        }

        //cout << endl;
    }

    if(pixel_set.size() < 500) {
        for(auto i : pixel_set) {
            color[i] = 0;
        }
        return true;
    }
    else 
        return false;

}

void average_rudy(float* rudy) {
    float* new_rudy = new float [xGrid * yGrid];

    for(int i = 0; i < xGrid * yGrid; i++) {
        new_rudy[i] = 0;
    }

    int size = 2;
    for(int x = size; x < xGrid - size; x += 2 * size + 1) {
        for(int y = size; y < yGrid - size; y += 2 * size + 1) {
            
            float tmp_rudy = 0;
            for(int newx = x - size; newx <= x + size; newx++) {
                for(int newy = y - size; newy <= y + size; newy++) {
                    tmp_rudy += rudy[newy * xGrid + newx];
                }
            }
            float div = (2 * size + 1) * (2 * size + 1);
            tmp_rudy = tmp_rudy / div;
            for(int newx = x - size; newx <= x + size; newx++) {
                for(int newy = y - size; newy <= y + size; newy++) {
                    new_rudy[newy * xGrid + newx] = tmp_rudy;
                }
            }
        }
    }

    for(int i = 0; i < xGrid * yGrid; i++) {
        rudy[i] = new_rudy[i];
    }

    delete[] new_rudy;
}

void find_closest_color(float* color, float* new_color, int x_in, int y_in) {

    std::queue<int> WL;
    WL.push(y_in * xGrid + x_in);
    std::set<int> pixel_visited;
    while(!WL.empty()) {
        int index = WL.front();
        WL.pop();
        int x = index % xGrid;
        int y = index / xGrid;
        if(color[index] != 0) {
            new_color[y_in * xGrid + x_in] = color[index];
            break;
        }

        if(y < yGrid - 1) {
            if(pixel_visited.count((y + 1) * xGrid + x) == 0) {
                WL.push((y + 1) * xGrid + x);
                pixel_visited.insert((y + 1) * xGrid + x);
            }
        }

        if(y > 0) {
            if(pixel_visited.count((y - 1) * xGrid + x) == 0) {
                WL.push((y - 1) * xGrid + x);
                pixel_visited.insert((y - 1) * xGrid + x);
            }
        }

        if(x < xGrid - 1) {
            if(pixel_visited.count(y * xGrid + x + 1) == 0) {
                WL.push(y * xGrid + x + 1);
                pixel_visited.insert(y * xGrid + x + 1);
            }
        }

        if(x > 0) {
            if(pixel_visited.count(y * xGrid + x - 1) == 0) {
                WL.push(y * xGrid + x - 1);
                pixel_visited.insert(y * xGrid + x - 1);
            }
        }

    }

}

void region_stats(Net** nets, int numValidNets, float* region_partition) {

    int bbox = 10;
    int partition_cnt = 30;
    int used_partition[partition_cnt];
    int nets_in_bbox = 0;
    int nets_in_same_color = 0, nets_not_in_same_color = 0;
    int in_big_regular_grid_cnt = 0;
    int big_regular_grid = 5;
    int big_regular_grid_x = xGrid / 5 + 1;
    int big_regular_grid_y = yGrid / 5 + 1;
    
    for(int i = 0; i < numValidNets; i++) {
        int npins = nets[i]->deg;
        int xmax = 0, xmin = xGrid - 1, ymax = 0, ymin = yGrid - 1;
        bool same_big_grid = true;
        int grid_x, grid_y;

        for(int cnt = 0; cnt < partition_cnt; cnt++) {
            used_partition[cnt] = 0;
        }
        for(int pinid = 0; pinid < npins; pinid++) {
            if(pinid == 0) {
                grid_x = nets[i]->pinX[pinid] / big_regular_grid_x;
                grid_y = nets[i]->pinY[pinid] / big_regular_grid_y;
            }
            else {
                if(grid_x != nets[i]->pinX[pinid] / big_regular_grid_x || grid_y != nets[i]->pinY[pinid] / big_regular_grid_y)
                    same_big_grid = false;
            }
            xmax = max(xmax, (int) nets[i]->pinX[pinid]);
            xmin = min(xmin, (int) nets[i]->pinX[pinid]);
            ymax = max(ymax, (int) nets[i]->pinY[pinid]);
            ymin = min(ymin, (int) nets[i]->pinY[pinid]);
            
            int location = (int) nets[i]->pinY[pinid] * xGrid + (int) nets[i]->pinX[pinid];
            if(location > xGrid * yGrid - 1) {
                cout << "location exceeds colored map: " << (int) nets[i]->pinX[pinid] << " " << (int) nets[i]->pinY[pinid] << endl;
            }

            int partition = (int) region_partition[location];
            used_partition[partition]++;
        }

        if(same_big_grid)
            in_big_regular_grid_cnt++;

        if(xmax - xmin <= bbox && ymax - ymin <= bbox) {
            nets_in_bbox++;
        }

        for(int cnt = 1; cnt < partition_cnt; cnt++) { //starting from 1
            if(used_partition[cnt] > 0) {
                if(used_partition[cnt] == npins) {
                    nets_in_same_color++;
                }
                else {
                    nets_not_in_same_color++;
                }
                break;
            }
        }

        
    }

    cout << "valid nets: " << numValidNets << " nets_in_bbox: " << nets_in_bbox << " nets_out_bbox: " << numValidNets - nets_in_bbox << endl;
    cout << " nets_in_big_grid: " << in_big_regular_grid_cnt << " nets_not_in_same_color: " << numValidNets - in_big_regular_grid_cnt << endl;
    cout << " nets_in_same_color: " << nets_in_same_color << " nets_not_in_same_color: " << nets_not_in_same_color << endl;


}

void region_stats_schedule_in_maze(Net** nets, int numValidNets, float* region_partition, galois::LargeArray<bool>& done, vector<vector<int>>& region_part) {

    int bbox = 10;
    int partition_cnt = 30;
    int used_partition[partition_cnt];
    int nets_in_bbox = 0;
    int nets_in_same_color = 0, nets_not_in_same_color = 0;
    int in_big_regular_grid_cnt = 0;
    int big_regular_grid = 5;
    
    int undone_cnt = 0;
    int done_cnt = 0;

    int stdcell_xmin = xGrid;
    int stdcell_xmax = 0;
    int stdcell_ymin = yGrid;
    int stdcell_ymax = 0;
    for(int x = 0; x < xGrid; x++) {
        for(int y = 0; y < yGrid; y++) {
            if(region_partition[y * xGrid + x] != 0) {
                stdcell_xmin = min(stdcell_xmin, x);
                stdcell_xmax = max(stdcell_xmax, x);
                stdcell_ymin = min(stdcell_ymin, y);
                stdcell_ymax = max(stdcell_ymax, y);
            }
        }
    }

    int big_regular_grid_x = (stdcell_xmax - stdcell_xmin) / 5 + 1;
    int big_regular_grid_y = (stdcell_ymax - stdcell_ymin) / 5 + 1;
    
    for(int i = 0; i < numValidNets; i++) {
        if(done[i]) {
            done_cnt++;
            continue;
        }
        else
            undone_cnt++;

        int npins = nets[i]->deg;
        int xmax = 0, xmin = xGrid - 1, ymax = 0, ymin = yGrid - 1;
        bool same_big_grid = true;
        int grid_x, grid_y;

        for(int cnt = 0; cnt < partition_cnt; cnt++) {
            used_partition[cnt] = 0;
        }

        for(int pinid = 0; pinid < npins; pinid++) {
            if(pinid == 0) {
                grid_x = (nets[i]->pinX[pinid] - stdcell_xmin) / big_regular_grid_x;
                grid_y = (nets[i]->pinY[pinid] - stdcell_ymin) / big_regular_grid_y;
            }
            else {
                if(grid_x != (nets[i]->pinX[pinid] - stdcell_xmin) / big_regular_grid_x || grid_y != (nets[i]->pinY[pinid] - stdcell_ymin) / big_regular_grid_y)
                    same_big_grid = false;
            }
            xmax = max(xmax, (int) nets[i]->pinX[pinid]);
            xmin = min(xmin, (int) nets[i]->pinX[pinid]);
            ymax = max(ymax, (int) nets[i]->pinY[pinid]);
            ymin = min(ymin, (int) nets[i]->pinY[pinid]);
            
            int location = (int) nets[i]->pinY[pinid] * xGrid + (int) nets[i]->pinX[pinid];
            if(location > xGrid * yGrid - 1) {
                cout << "location exceeds colored map: " << (int) nets[i]->pinX[pinid] << " " << (int) nets[i]->pinY[pinid] << endl;
            }

            int partition = (int) region_partition[location];
            used_partition[partition]++;
        }

        if(same_big_grid)
            in_big_regular_grid_cnt++;

        if(xmax - xmin <= bbox && ymax - ymin <= bbox) {
            nets_in_bbox++;
        }

        for(int cnt = 1; cnt < partition_cnt; cnt++) { //starting from 1
            if(used_partition[cnt] > 0) {
                if(used_partition[cnt] == npins) {
                    nets_in_same_color++;
                    region_part[cnt].push_back(i); //schedule to region
                }
                else {
                    nets_not_in_same_color++;
                }
                break;
            }
        }
    }

    cout << "valid nets: " << numValidNets << " done: " << done_cnt << " undone: " << undone_cnt << endl;
    cout << "valid nets: " << numValidNets << " nets_in_small_bbox: " << nets_in_bbox << " nets_out_bbox: " << undone_cnt - nets_in_bbox << endl;
    cout << " nets_in_big_grid: " << in_big_regular_grid_cnt << " nets_not_in_big_grid: " << undone_cnt - in_big_regular_grid_cnt << endl;
    cout << " nets_in_same_color: " << nets_in_same_color << " nets_not_in_same_color: " << nets_not_in_same_color << endl;
}

void rudy_region_partition(float* rudy, Algo algo, float* region_partition) {

    plot_rudy(rudy, algo);
    average_rudy(rudy);

    float rudy_sum = 0;
    int color_cnt = 1;
    //float* color = new float [xGrid * yGrid];
    float* color = region_partition;

    for(int i = 0; i < xGrid * yGrid; i++) {
        color[i] = 0;
    }

    float color_threshold = 12.5;
    for(int x = 0; x < xGrid; x++) {
        for(int y = 0; y < yGrid; y++) {

            if(rudy[y * xGrid + x] > color_threshold && color[y * xGrid + x] == 0) {
                bool abort_color = paint_color(color, rudy, x, y, color_cnt, color_threshold);
                if(!abort_color)
                    color_cnt++;
            }
        }
    }

    int stdcell_xmin = xGrid;
    int stdcell_xmax = 0;
    int stdcell_ymin = yGrid;
    int stdcell_ymax = 0;
    float* new_color = new float [xGrid * yGrid]; //including region not colored yet
    for(int x = 0; x < xGrid; x++) {
        for(int y = 0; y < yGrid; y++) {
            new_color[y * xGrid + x] = 0;
            if(color[y * xGrid + x] != 0) {
                stdcell_xmin = min(stdcell_xmin, x);
                stdcell_xmax = max(stdcell_xmax, x);
                stdcell_ymin = min(stdcell_ymin, y);
                stdcell_ymax = max(stdcell_ymax, y);
            }
        }
    }
    
    if(algo == DetPart_Astar_Region) {
        for(int x = 0; x < xGrid; x++) {
            for(int y = 0; y < yGrid; y++) {
                if(color[y * xGrid + x] == 0 && x >= stdcell_xmin && x <= stdcell_xmax && y >= stdcell_ymin && y <= stdcell_ymax)
                    find_closest_color(color, new_color, x, y); //this can be improved by expanding the current color
            }
        }

        for(int x = 0; x < xGrid; x++) {
            for(int y = 0; y < yGrid; y++) {
                if(color[y * xGrid + x] == 0 && x >= stdcell_xmin && x <= stdcell_xmax && y >= stdcell_ymin && y <= stdcell_ymax)
                    color[y * xGrid + x] = new_color[y * xGrid + x];
            }
        }
        cout << "total number of colors: " << color_cnt << endl;
        num_region_partitions = color_cnt;
    }
    else if(algo == DetPart_Astar_Regular_Region) {
        int big_regular_grid = 5;
        int big_regular_grid_x = (stdcell_xmax - stdcell_xmin) / big_regular_grid + 1;
        int big_regular_grid_y = (stdcell_ymax - stdcell_ymin) / big_regular_grid + 1;
        num_region_partitions = big_regular_grid * big_regular_grid + 1;

        for(int x = stdcell_xmin; x <= stdcell_xmax; x++) {
            for(int y = stdcell_ymin; y <= stdcell_ymax; y++) {
                int grid_x = (x - stdcell_xmin) / big_regular_grid_x;
                int grid_y = (y - stdcell_ymin) / big_regular_grid_y;
                color[y * xGrid + x] = grid_y * big_regular_grid + grid_x + 1;
            }
        }

    }

    ofstream outfile("heat.txt");
    outfile << "heat" << endl;
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
            if(color[y * xGrid + x] > 0)
                outfile << color[y * xGrid + x] << endl;
            else
                outfile << 0 << endl;
        }
    }
    outfile << ")" << endl;
    outfile.close();

    //delete [] color;
    delete [] new_color;
}

float plot_pin_density(float* pin_density, Algo algo) {

    cout << xGrid << " " << yGrid << " grid" << endl;
    
    float pin_density_sum = 0;
    for(int i = 0; i < xGrid * yGrid; i++)
        pin_density[i] = 0;

    for(int netID = 0; netID < numValidNets; netID++) {

        int deg = nets[netID]->deg;
        for(int pinID = 0; pinID < deg; pinID++) {
            int x = (int) nets[netID]->pinX[pinID];
            int y = (int) nets[netID]->pinY[pinID];
            if(x != 0 && y != 0 && x != xGrid - 1 && y != yGrid - 1) {
                pin_density[y * xGrid + x] += 1.0;
                pin_density_sum += 1.0;
            }
        }
    }

    for(int netID = 0; netID < numInvalidNets; netID++) {
        int x = (int) invalid_nets[netID]->pinX[0];
        int y = (int) invalid_nets[netID]->pinY[0];
        if(x != 0 && y != 0 && x != xGrid - 1 && y != yGrid - 1) {
            pin_density[y * xGrid + x] += 2.0;
            pin_density_sum += 2.0;
        }
    }

    int* expand_pin_density = new int [xGrid * yGrid]; // Just to plot, not for space reservation
    for(int x = 0; x < xGrid; x++) {
        for(int y = 0; y < yGrid; y++) {
            expand_pin_density[y * xGrid + x] = 0;
        }
    }

    for(int x = 0; x < xGrid; x++) {
        for(int y = 0; y < yGrid; y++) {
            
            for(int local_x = max(0, x - 2); local_x < min(xGrid - 1, x + 2); local_x++) {
                for(int local_y = max(0, y - 2); local_y < min(yGrid - 1, y + 2); local_y++) {
                    expand_pin_density[y * xGrid + x] += pin_density[local_y * xGrid + local_x];
                }
            }
        }
    }

    for(int x = 0; x < xGrid; x++) {
        for(int y = 0; y < yGrid; y++) {

            expand_pin_density[y * xGrid + x] /= 25.0;

        }
    }

    if(algo == PIN_DENSITY) {
        
        int hist_size = 100;
        int hist[hist_size];
        for(int i = 0; i < hist_size; i++){
            hist[i] = 0;
        }
        
        ofstream outfile("heat.txt");
        cout << "heat? " << endl;
        outfile << "heat" << endl;
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
                //expand_pin_density[y * xGrid + x] = (expand_pin_density[y * xGrid + x] < 60)? 0 : expand_pin_density[y * xGrid + x];
                outfile << expand_pin_density[y * xGrid + x]  << endl;
                float pin_density1 = expand_pin_density[y * xGrid + x] ;
                pin_density_sum += pin_density1;
                if((int) pin_density1 < hist_size)
                    hist[(int) pin_density1]++;
                else hist[hist_size - 1]++;

                //if(pin_density1 >= 20)
                //    cout << "high density point > 20 : " << x << " " << y << endl;
                
            }
        }
        outfile << ")" << endl;
        outfile.close();

        for(int i = 0; i < hist_size; i++) {
            cout  << hist[i] << endl;
        }

        cout << "avg: " << pin_density_sum / (float)(xGrid * yGrid) << endl;
    }

    return pin_density_sum / (float)(xGrid * yGrid);

}


void plot_drc_map() {

    ifstream infile("ispd19_test7_metal5.sproute.georpt");

    if(!infile.is_open()) {
        cout << "unable to open rpt file" << endl;
        exit(1);
    }
    
    int* drc_cnt = new int [xGrid * yGrid];

    for(int i = 0; i < xGrid*yGrid; i++)
        drc_cnt[i] = 0;

    int dbuPerMicro = defDB.dbuPerMicro;
    int total_cnt = 0;
    while(!infile.eof()) {
        string tmp;
        infile >> tmp;
        if(tmp == "Metal2") {
            total_cnt++;
            float x1, y1, x2, y2;
            string tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
            infile >> tmp0 >> tmp1;
            infile >> tmp2 >> tmp3 >> x1 >> tmp4 >> y1 >> tmp5 >> tmp6 >> x2 >> tmp7 >> y2 >> tmp8;

            if(tmp0 != ")" || tmp2 != ":" || tmp3 != "(" || tmp4 != "," || tmp5 != ")" || tmp6 != "(" || tmp7 != "," ||
            tmp8 != ")" ) {
                cout << tmp0 << " " << tmp1 << " " << tmp2 <<" " << tmp3 <<" " << tmp4 <<" " << tmp5 <<" " << tmp6;
                cout << " " << tmp7 << " " << tmp8;
                cout << endl;
                cout << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
                exit(1);
            }
            x1 *= dbuPerMicro;
            x2 *= dbuPerMicro;
            y1 *= dbuPerMicro;
            y2 *= dbuPerMicro;

            int x_start = parser::find_Gcell((int) x1, defDB.xGcellBoundaries); //all are inclusive
            int x_end   = parser::find_Gcell((int) x2, defDB.xGcellBoundaries);
            int y_start = parser::find_Gcell((int) y1, defDB.yGcellBoundaries);
            int y_end   = parser::find_Gcell((int) y2, defDB.yGcellBoundaries);

            int x_mid = (x_start + x_end) / 2;
            int y_mid = (y_start + y_end) / 2;
            /*for(int x = x_start; x <= x_end; x++) {
                for(int y = y_start; y <= y_end; y++) {
                    drc_cnt[y * xGrid + x]++;
                }
            }*/
            drc_cnt[y_mid * xGrid + x_mid]++;

        }
    }
    cout << "total number of drcs: " << total_cnt << endl;


    int* expand_drc_cnt = new int [xGrid * yGrid];
    for(int x = 0; x < xGrid; x++) {
        for(int y = 0; y < yGrid; y++) {
            expand_drc_cnt[y * xGrid + x] = 0;
        }
    }
    for(int x = 0; x < xGrid; x++) {
        for(int y = 0; y < yGrid; y++) {
            
            for(int local_x = max(0, x - 2); local_x < min(xGrid - 1, x + 2); local_x++) {
                for(int local_y = max(0, y - 2); local_y < min(yGrid - 1, y + 2); local_y++) {
                    expand_drc_cnt[y * xGrid + x] += drc_cnt[local_y * xGrid + local_x];
                }
            }
        }
    }

    

    int hist_size = 20;
    int hist[hist_size];
    for(int i = 0; i < hist_size; i++){
        hist[i] = 0;
    }
        
    ofstream outfile("heat.txt");
    cout << "heat? " << endl;
    outfile << "heat" << endl;
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
            outfile << expand_drc_cnt[y * xGrid + x] << endl;
            int drc_cnt1 = expand_drc_cnt[y * xGrid + x];
            if(drc_cnt1 < hist_size)
                hist[drc_cnt1]++;
            else hist[hist_size - 1]++;
        }
    }
    outfile << ")" << endl;
    outfile.close();
    infile.close();

    for(int i = 0; i < hist_size; i++) {
        cout  << hist[i] << endl;
    }

}

#endif