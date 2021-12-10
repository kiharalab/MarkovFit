//Created by Eman Alnabati on 12/4/17
//Analyze Results

#include "handle_options.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <sstream>
#include <iterator>
using namespace std;

void rad2degree(double rad[3],double degree[3]){
    for (int i=0; i<3; i++) {
        degree[i]= rad[i] * 180.0 / M_PI; // conversion to degrees
        if( degree[i] < 0 ) degree[i] += 360.0; // convert negative to positive angles
    }
}

void Q2mtx(double q[4],double mtx[3][3]){
    mtx[0][0]=1-2*(q[2]*q[2]+q[3]*q[3]);
    mtx[0][1]=2*(q[1]*q[2]-q[0]*q[3]);
    mtx[0][2]=2*(q[1]*q[3]+q[0]*q[2]);
    
    mtx[1][0]=2*(q[2]*q[1]+q[0]*q[3]);
    mtx[1][1]=1-2*(q[3]*q[3]+q[1]*q[1]);
    mtx[1][2]=2*(q[2]*q[3]-q[0]*q[1]);
    
    mtx[2][0]=2*(q[3]*q[1]-q[0]*q[2]);
    mtx[2][1]=2*(q[3]*q[2]+q[0]*q[1]);
    mtx[2][2]=1-2*(q[1]*q[1]+q[2]*q[2]);
}

void mtx2euler(double q[3],double rot_mat[3][3]){
    q[0] = atan2(rot_mat[2][1], rot_mat[2][2]);
    q[1] = atan2(-rot_mat[2][0],sqrt(pow(rot_mat[2][1],2)+pow(rot_mat[2][2],2)));
    q[2] = atan2(rot_mat[1][0], rot_mat[0][0]);
}

bool sort_CC (const vector <float>& v1, const vector <float>& v2) {
    if (v1[4] < v2[4]) {
        return true;
    } else if (v1[4] > v2[4]) {
        return false;
    }
}

bool sort_OV (const vector <float>& v1, const vector <float>& v2) {
    if (v1[5] < v2[5]) {
        return true;
    } else if (v1[5] > v2[5]) {
        return false;
    }
}

bool sort_CC_OV (const vector <float>& v1, const vector <float>& v2) {
    if (v1[6] > v2[6]) {
        return true;
    } else if (v1[6] < v2[6]) {
        return false;
    }
}

bool integral (float a, float b, float c) {
    return (a==b==c);
}

float find_median (vector <vector<float> > A ,int col) {
    if (A.size() % 2 == 0) {
        return (A[A.size()/2][col] + A[A.size()/2-1][col])/2.0;
    } else {
        return A[A.size()/2][col];
    }
}

float find_MAD (vector <vector<float> > A ,int col, float median) {
    vector <float> B;
    float tmp;
    for (int i=0 ; i < A.size() ; i++) {
        tmp = fabs(A[i][col]-median);
        B.push_back(tmp);
    }
    sort(B.begin(), B.end());
    if (A.size() % 2 == 0) {
        return (B[A.size()/2] + B[A.size()/2-1])/2.0;
    } else {
        return B[A.size()/2];
    }
}

int main(int argc, char** argv)
{
    options options(argc,argv);
    if(!options.parse_successful()) {
        cout << options.usage();
        return 0;
    }
    
    float x = options.get_x_value();
    float y = options.get_y_value();
    float z = options.get_z_value();
    float dist_threshold = options.get_dist_threshold();
    float min_dist = options.get_min_dist();
    string output_filename = options.get_output_filename();
    string input_filename = options.get_input_filename();
    
    //2D array to store all results
    vector <vector<float> > result;
    vector <float> tmp (10);
    ifstream current_file;
    int counter = 0;
    
    //Read  input file and populate the 2D array
    current_file.open(input_filename.c_str());
    if (!current_file.is_open()) {
        cerr << "Could not open file." <<endl;
    } else {
        cout << "Open file " << input_filename.c_str() <<"\n";
    }
    string a, b, c, d, e, f, g, h, l, m, line;
    float distance;
    
    double q[4];
    double mtx[3][3], rot_mtx[3][3];
    double dd[3], rot_d[3], degree[3], rot_degree[3];
    
    tmp.clear();
    while (getline(current_file,line)) {
        istringstream iss(line);
        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};
        
        if(line[0]=='*' && atof(tokens[7].c_str()) >= 0.50) {
            
            tmp.push_back(atof(tokens[2].c_str()));//0
            tmp.push_back(atof(tokens[3].c_str()));//1
            tmp.push_back(atof(tokens[4].c_str()));//2
            if (x != FLT_MAX && y != FLT_MAX && z != FLT_MAX) {
                distance = sqrt(pow(atof(tokens[2].c_str())-x,2)+pow(atof(tokens[3].c_str())-y,2)+pow(atof(tokens[4].c_str())-z,2));}
            else {
                distance = 0.000;
            }
            tmp.push_back(distance);//3
            tmp.push_back(atof(tokens[9].c_str()));//4 CC
            tmp.push_back(atof(tokens[7].c_str()));//5 OV
            
            tmp.push_back(0);//6 CC & OV
            
            q[0] = atof(tokens[25].c_str());
            q[1] = atof(tokens[26].c_str());
            q[2] = atof(tokens[27].c_str());
            q[3] = atof(tokens[28].c_str());
            
            Q2mtx(q,mtx);
            mtx2euler(dd,mtx);
            rad2degree(dd,degree);
            
            if (degree[0] == -0) {
                degree[0] = 0;
            }
            if (degree[1] == -0) {
                degree[1] = 0;
            }
            if (degree[2] == -0) {
                degree[2] = 0;
            }
            
            tmp.push_back(degree[0]);//7
            tmp.push_back(degree[1]);//8
            tmp.push_back(degree[2]);//9
            
            result.push_back(tmp);
            tmp.clear();
            counter++;
        }
    }
    current_file.close();
    current_file.clear();
    
    cout << "Size of search results = " << result.size() << "\n";
    
    float median_CC, median_OV, MAD_CC, MAD_OV, scale_factor, shift_factor;
    
    sort(result.begin(), result.end(), sort_CC);
    median_CC = find_median(result, 4);
    sort(result.begin(), result.end(), sort_OV);
    median_OV = find_median(result, 5);
    shift_factor = median_CC - median_OV;
    MAD_CC = find_MAD(result, 4, median_CC);
    MAD_OV = find_MAD(result, 5, median_OV);
    scale_factor = MAD_CC / MAD_OV;
    
    cout << " median_CC = " << median_CC << " median_OV = " << median_OV << " shift_factor = " << shift_factor << " MAD_CC = " << MAD_CC << " MAD_OV = " << MAD_OV << " scale_factor = " << scale_factor << "\n";
    
    //Combine CC and OV
    for (int i=0 ; i < result.size() ; i++) {
        result[i][6] = (result[i][5] * scale_factor) + shift_factor;
        result[i][6] = (result[i][4] + result[i][6]) / 2.0;
    }
    
    // Sort results
    sort(result.begin(), result.end(), sort_CC_OV);
    
    
    //Stats about search results
    if (x != FLT_MAX && y != FLT_MAX && z != FLT_MAX) {
        int count = 0;
        float smallest_ele = result[0][3];
        int smallest_rank = 1;
        
        //Find rank of smallest distance
        cout << "Rank of top 100 predicted positions with distance < threshold to the native structure after sorting based CC and OVR: \n";
        for (int i=0 ; i< result.size() ; i++) {
            if (result[i][3] <= dist_threshold && count < 100) {
                count++;
                cout << "Rank = " << i+1 << "\t Distance " << result[i][3] << "\t X = " << result[i][0] << " Y = "
                << result[i][1] << " Z = " << result[i][2] << " RX = "
                << result[i][7] << " RY = " << result[i][8] << " RZ = " << result[i][9] << " CC = " << result[i][4] << " OV = "
                << result[i][5] << " CC&OV = " << result[i][6] << "\n";
            }
            
            if (result[i][3] < smallest_ele) {
                smallest_ele = result[i][3];
                smallest_rank = i+1;
            }
        }
        
        cout << "smallest distance = " << smallest_ele << " rank = "<<smallest_rank<< "\n";
    }
    
    //Write result to file
    
    
    string s = output_filename.c_str();
    s.replace(s.end()-4,s.end(),"_b4_clustering.txt");
    cout << s <<endl;
    
    ofstream results_file (s);
    if (!results_file.is_open()) {
        cerr << "Could not create result file!" << endl;
        throw 1;
    }
    
    for (int i=0 ; i< result.size() ; i++) {
        results_file.unsetf(std::ios_base::floatfield);
        results_file << result[i][3] << " " << result[i][0] << " "
        << result[i][1] << " " << result[i][2] << " "
        << result[i][7] << " " << result[i][8] << " " << result[i][9] << " ";
        results_file << std::fixed << result[i][4] << " " << result[i][5]
        << " " << result[i][6];
        
        results_file << endl;
    }
    results_file.close();
    
    //Cluster FFT search results
    bool flag;
    float c1 [3], c2 [3];
    vector<vector<float> > filtered_result;
    filtered_result.clear();
    
    for (int i=0; i<result.size(); i++) {
        c1[0] = (result[i][0]);
        c1[1] = (result[i][1]);
        c1[2] = (result[i][2]);
        flag = true;
        
        for (int j=0; j<filtered_result.size(); j++) {
            c2[0] = (filtered_result[j][0]);
            c2[1] = (filtered_result[j][1]);
            c2[2] = (filtered_result[j][2]);
            
            distance = sqrt(pow(c2[0]-c1[0],2)+pow(c2[1]-c1[1],2)+pow(c2[2]-c1[2],2));
            if(distance < min_dist) { //clash
                flag = false;
                break;
            } //clash
        }//end j loop
        if (flag) {
            filtered_result.push_back(result[i]);
        }
    }
    
    //Print clustered search result
    int count = 0;
    
    cout << "Rank of of clustered search result with distance < threshold to the native structure: \n";
    
    for(int j=0; j< filtered_result.size(); j++) {
        if (filtered_result[j][3]<= dist_threshold) {
            cout<<"rank = "<<j+1<<" : "<<filtered_result[j][3] << " " <<filtered_result[j][0] << " "
            <<filtered_result[j][1] << " "<<filtered_result[j][2] << " "
            <<filtered_result[j][7] << " "<<filtered_result[j][8] << " "
            <<filtered_result[j][9] << " "<<filtered_result[j][4] << " "
            <<filtered_result[j][5] << " "<<filtered_result[j][6] << endl;
            count++;
        }
    }
    
    //Write clustered search result to a file
    ofstream cluster_file (output_filename.c_str());

    if (!cluster_file.is_open()) {
        cerr << "Could not create result file!" << endl;
        throw 1;
    }
    for (int j=0 ; j< filtered_result.size() ; j++) {
        cluster_file.unsetf(std::ios_base::floatfield);
        cluster_file //<< std::fixed
        << filtered_result[j][3] << "\t" << filtered_result[j][0] << "\t" << filtered_result[j][1]
        << "\t" << filtered_result[j][2] << "\t" << filtered_result[j][7] << "\t"
        << filtered_result[j][8] << "\t" <<filtered_result[j][9] << "\t" << filtered_result[j][4]
        << "\t" <<filtered_result[j][5] << "\t" <<filtered_result[j][6] << endl;
    }
    cluster_file.close();
    
    filtered_result.clear();
    //End of chain loop
    
}
