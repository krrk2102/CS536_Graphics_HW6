//
//  main.cpp
//  CG_hw6
//
//  Created by Shangqi Wu on 14/10/31.
//  Copyright (c) 2014 Shangqi Wu. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>
#include <limits>

using namespace std;

class Point;
class Line;
class Curve;

// Define global variant----------------------------------------------------------------------------
int a=0, b=0, c=250, d=250, m=0, n=0, r=0, j=0, k=0, o=250, p=250; // These vars are used to be attributes of the image and window.
float s = 1, L = 0.05;
string input = "./ExtraCredit.ps", output = "./out.xpm";
vector<vector<char> > xpmpix; // This vector represents status of all pixels in the world window.
const double pi = acos((double)-1);

//Define all functions for Cohen-Sutherland algorithm----------------------------------------------
int cstest(Line argl, int a, int b, int c, int d);
Line csclip(Line argl, int a, int b, int c, int d);
// Define all basic functions------------------------------------------------------------------------
void help();
void optana (int argc, char * const argv[]);
string setheader();
string setend();
void outfile(string output);
Curve judgecurves(Curve cur);
int rnd(float arg);

// Define classes needed-----------------------------------------------------------------------------
class Point {
private:
    int x;
    int y;
public:
    void set(int argx, int argy) {
        x = argx;
        y = argy;
    }
    int getx() {
        return x;
    }
    int gety() {
        return y;
    }
    void trans(int m, int n) { // Tanslation
        int ptold[3] = {x, y, 1};
        int matrix[3][3] = {1, 0, m, 0, 1, n, 0, 0, 1};
        int pt[3] = {0, 0, 0};
        for (int i = 0; i < 3; i++) {
            pt[i] += matrix[i][0] * ptold[0];
            pt[i] += matrix[i][1] * ptold[1];
            pt[i] += matrix[i][2] * ptold[2];
        }
        x=pt[0]; y=pt[1];
    }
    void scale(float sx, float sy) { // Scaling
        float ptold[3] = {(float)x, (float)y, 1};
        float matrix[3][3] = {sx, 0, 0, 0, sy, 0, 0, 0, 1};
        float pt[3] = {0, 0, 0};
        for (int i = 0; i < 3  ; i++) {
            pt[i] += matrix[i][0] * ptold[0];
            pt[i] += matrix[i][1] * ptold[1];
            pt[i] += matrix[i][2] * ptold[2];
        }
        x=rnd(pt[0]); y=rnd(pt[1]);
    }
    void rot(int r) { // Rotation
        double ang = (double)r * pi / 180.0;
        float ptold[3] = {(float)x, (float)y, (float)1};
        float pt[3] = {0, 0, 0};
        float matrix[3][3] = {(float)cos(ang), (float)-sin(ang), 0, (float)sin(ang), (float)cos(ang), 0, 0, 0, 1};
        for (int i = 0; i < 3; i++) {
            pt[i] += matrix[i][0] * ptold[0];
            pt[i] += matrix[i][1] * ptold[1];
            pt[i] += matrix[i][2] * ptold[2];
        }
        x=rnd(pt[0]); y=rnd(pt[1]);
    }
    bool operator==(Point a){
        return (x==a.getx() && y==a.gety());
    }
    bool operator!=(Point a){
        return (x!=a.getx() || y!=a.gety());
    }
};

class Line { // present the line by: y = slope*x + d
private:
    Point start;
    Point end;
    int xmax;
    int xmaxy;
    int xmin;
    int xminy;
    int ymax;
    int ymin;
    bool sd_exist;
    float slope;
    float d1;
    bool in;
    void cal() {
        if (start.getx() != end.getx()) {
            if (start.getx()>=end.getx()) {
                xmax = start.getx();
                xmaxy = start.gety();
                xmin = end.getx();
                xminy = end.gety();
            }else {
                xmax = end.getx();
                xmaxy = end.gety();
                xmin = start.getx();
                xminy = start.gety();
            }
            if (start.gety()>=end.gety()) { // set ymax and ymin
                ymax = start.gety();
                ymin = end.gety();
            }else {
                ymax = end.gety();
                ymin = start.gety();
            }
            slope = (float)(start.gety() - end.gety()) / (float)(start.getx() - end.getx());
            d1 = (float)start.gety() - (slope * (float)start.getx());
            sd_exist = true;
        }else {
            if (start.gety()>=end.gety()) {
                ymax = start.gety();
                ymin = end.gety();
                xmaxy = start.gety();
                xminy = end.gety();
            }else {
                ymax = end.gety();
                ymin = start.gety();
                xmaxy = end.gety();
                xminy = start.gety();
            }
            xmax = start.getx();
            xmin = end.getx();
            slope = numeric_limits<float>::max();
            d1 = numeric_limits<float>::max();
            sd_exist = false;
        }
    }
public:
    void setall(int argx1, int argy1, int argx2, int argy2) {
        start.set(argx1, argy1);
        end.set(argx2, argy2);
        cal();
    }
    void setp(Point a, Point b) {
        start = a;
        end = b;
        cal();
    }
    bool sd() {
        return sd_exist;
    }
    Point gets() {
        return start;
    }
    Point gete() {
        return end;
    }
    Point cal_inter_wlr(float x) { // return the intersaction point of x=x
        Point tmp;
        if (sd_exist) {
            float y = slope*x +d1;
            tmp.set(rnd(x), rnd(y));
        }else {
            tmp.set(numeric_limits<int>::max(), numeric_limits<int>::max());
        }
        return tmp;
    }
    Point cal_inter_wtb(float y) { // return the intersaction point of y=y
        Point tmp;
        if (sd_exist) {
            float x = (y - d1)/ slope;
            tmp.set(rnd(x), rnd(y));
        }else {
            tmp.set(start.getx(), (int)y);
        }
        return tmp;
    }
    void setin(bool ifin) {
     in = ifin;
    }
    bool getin() {
     return in;
    }
    float getslope() {
        return slope;
    }
    int getxmin() {
        return xmin;
    }
    int getxmax() {
        return xmax;
    }
    int getymin() {
        return ymin;
    }
    int getymax() {
        return ymax;
    }
    int getxmaxy() {
        return xmaxy;
    }
    int getxminy() {
        return xminy;
    }
    void trans(int m,int n) { // Translation for both ends of the line.
        start.trans(m, n);
        end.trans(m, n);
        cal();
    }
    void scale(float sx, float sy) { // Scaling for both ends of the line.
        start.scale(sx,sy);
        end.scale(sx,sy);
        cal();
    }
    void rot(int r) { //Rotating for both ends of the line.
        start.rot(r);
        end.rot(r);
        cal();
    }
    void toviewport(){
        float sx = (float)(o - j)/(float)(c - a);
        float sy = (float)(p - k)/(float)(d - b);
        start.trans(-a,-b);
        if (((o-j)!=(c-a)) || ((p-k)!=(d-b))) {
            start.scale(sx,sy);
        }
        start.trans(j,k);
        end.trans(-a,-b);
        if (((o-j)!=(c-a)) || ((p-k)!=(d-b))) {
            end.scale(sx,sy);
        }
        end.trans(j,k);
        cal();
    }
    void showall() { // Calculate all the points of the line.
        if (sd_exist == true) {
            // Fllowing codes are of Bresenham Algorithm, presented in L-02_Lines.pdf. Codes are modifiied for this cpp file.
            int dx, dy, D, x, y;
            dx = xmax - xmin;
            dy = ymax - ymin;
            if (0<slope && slope<1) {
                D = 2*dy - dx;
                y = ymin;
                for (x = xmin; x <= xmax; x++) {
                    xpmpix[y][x] = '+';
                    if (D <= 0) {
                        D += 2*dy;
                    }else {
                        D += 2*(dy - dx);
                        y++;
                    }
                }
            }else if (slope > 1) {
                D = 2*dx - dy;
                x = xmin;
                for (y = ymin; y <= ymax; y++) {
                    xpmpix[y][x] = '+';
                    if (D <= 0) {
                        D += 2*dx;
                    }else {
                        D += 2*(dx - dy);
                        x++;
                    }
                }
            }else if (-1<slope && slope<0) {
                D = 2*dy - dx;
                y = ymax;
                for (x = xmin; x <= xmax; x++) {
                    xpmpix[y][x] = '+';
                    if (D <= 0) {
                        D += 2*dy;
                    }else {
                        D += 2*(dy - dx);
                        y--;
                    }
                }
            }else if (slope == 1) {
                y = ymin;
                for (x = xmin; x <= xmax; x++) {
                    xpmpix[y][x]='+';
                    y++;
                }
            }else if (slope == -1) {
                y = ymax;
                for (x = xmin; x <= xmax; x++) {
                    xpmpix[y][x] = '+';
                    y--;
                }
            }else if (slope == 0) {
                y = ymin;
                for (x = xmin; x <= xmax; x++) {
                    xpmpix[y][x] = '+';
                }
            }
            else { // i.e., slope<-1
                D = 2*dx - abs(dy);
                x = xmin;
                for (y = ymax; y >= ymin; y--) {
                    xpmpix[y][x] = '+';
                    if (D <= 0) {
                        D += 2*dx;
                    }else {
                        D += 2*(dx - dy);
                        x++;
                    }
                }
            }
        }else if (sd_exist == false) { // for vertical lines
            int x = xmin;
            for (int y = ymin; y <= ymax; y++) {
                xpmpix[y][x] = '+';
            }
        }
    }
};

class Curve{
private:
    vector<Point> ctrlpts;
    bool in;
    /*Point deCast(vector<Point> vecp, float l) { // Here is code implementing de Casteljau Alorithm.
     Point pt_new;
     if (vecp.size() > 2) {
     vector<Point> new_vecp;
     for (int i = 0; i < vecp.size()-1; i++) {
     int x_new = (int)((1-l)*vecp[i].getx() + l*vecp[i+1].getx());
     int y_new = (int)((1-l)*vecp[i].gety() + l*vecp[i+1].gety());
     pt_new.set(x_new,y_new);
     new_vecp.push_back(pt_new);
     }
     pt_new = deCast(new_vecp, l);
     }else if (vecp.size() == 2) {
     int x_new = (int)((1-l)*vecp[0].getx() + l*vecp[1].getx());
     int y_new = (int)((1-l)*vecp[0].gety() + l*vecp[1].gety());
     pt_new.set(x_new, y_new);
     }
     return pt_new;
     }*/
    vector<float> deCastel(vector<float> vecfx, vector<float> vecfy, float l) { // Here is code implementing de Casteljau Alorithm.
        vector<float> pt_new;
        if (vecfx.size() > 2) {
            vector<float> vecfx_new;
            vector<float> vecfy_new;
            for (int i = 0; i < vecfx.size()-1; i++) {
                float x_new = (1-l)*vecfx[i] + l*vecfx[i+1];
                float y_new = (1-l)*vecfy[i] + l*vecfy[i+1];
                vecfx_new.push_back(x_new);
                vecfy_new.push_back(y_new);
            }
            pt_new = deCastel(vecfx_new, vecfy_new, l);
        }else if (vecfx.size() == 2) {
            float x_new = (1-l)*vecfx[0] + l*vecfx[1];
            float y_new = (1-l)*vecfy[0] + l*vecfy[1];
            pt_new.push_back(x_new);
            pt_new.push_back(y_new);
        }
        return pt_new;
    }
public:
    void set(vector<Point> argp){
        if (argp.size() == 4) {
            ctrlpts.clear();
            ctrlpts = argp;
        }
    }
    Point getp(int i){
        if (i <= 3) {
            return ctrlpts[i];
        }else return ctrlpts[3];
    }
    void trans(int m, int n){
        for (int i = 0; i < 4; i++) {
            ctrlpts[i].trans(m, n);
        }
    }
    void rot(int r){
        for (int i = 0; i < 4; i++) {
            ctrlpts[i].rot(r);
        }
    }
    void scale(float sx, float sy){
        for (int i = 0; i < 4; i++) {
            ctrlpts[i].scale(sx, sy);
        }
    }
    void toviewport(){
        float sx = (float)(o - j)/(float)(c - a);
        float sy = (float)(p - k)/(float)(d - b);
        for (int i = 0; i < 4; i++) {
            ctrlpts[i].trans(-a, -b);
            if (((o-j)!=(c-a)) || ((p-k)!=(d-b))) {
                ctrlpts[i].scale(sx, sy);
            }
            ctrlpts[i].trans(j, k);
        }
    }
    void setin(bool a){
        in = a;
    }
    bool getin(){
        return in;
    }
    void draw(){
        vector<float> buff1(2), buff2(2);
        vector<float> ctrlptsx(4);
        vector<float> ctrlptsy(4);
        for (int i = 0; i < 4; i++) {
            ctrlptsx[i] = ctrlpts[i].getx();
            ctrlptsy[i] = ctrlpts[i].gety();
        }
        Line buffl;
        buff2[0] = ctrlpts[0].getx(); buff2[1] = ctrlpts[0].gety();
        float limit = 1 + 0.5*L;
        for (float u = 0; u < limit; u+=L) {
            buff1 = buff2;
            buff2 = deCastel(ctrlptsx, ctrlptsy, u);
            buffl.setall((int)(buff1[0]), (int)(buff1[1]), (int)(buff2[0]), (int)(buff2[1]));
            buffl = csclip(buffl, j, k, o, p);
            if (buffl.getin() == true) {
                buffl.showall();
            }
        }
        /*Point buff1, buff2;
         buff2.set(ctrlpts[0].getx(), ctrlpts[0].gety());
         Line buffl;
         //float limit = 1 + 0.5*L;
         for (float u = L; u < 1; u+=L) {
         buff1.set(buff2.getx(), buff2.gety());
         buff2 = deCast(ctrlpts, u);
         buffl.setp(buff1, buff2);
         buffl = csclip(buffl, j, k, o, p);
         if (buffl.getin() == true) {
         buffl.showall();
         }
         }
         buff1.set(ctrlpts[3].getx(), ctrlpts[3].gety());
         buffl.setp(buff2, buff1);
         buffl = csclip(buffl, j, k, o, p);
         if (buffl.getin() == true) {
         buffl.showall();
         }*/
    }
};

//Here is main function---------------------------------------------------------------------------
int main(int argc, char * argv[]) {
    // analyze all the input options
    optana(argc, argv);
    int width = 501; // height of the world window
    int height = 501; // width of the world window
    xpmpix.resize(height); // set vector as a 2-d array
    for (int i = 0; i < height; i++) { // Initailize all char to be '-', which stands for white pixel
        xpmpix[i].resize(width);
        for (int j = 0; j < width; j++) { // i.e., all pixels should be white at first
            xpmpix[i][j] = '-'; // stands for point(y,x)
        }
    }
    // read and buffer input *.ps file
    input = "/Users/wushangqi/ExtraCredit.ps"; // dir only for debug*********************************
    vector<Curve> veccurve;
    ifstream infile(input.c_str());
    if (!infile) {
        cout<<"File does not exist, please check your path."<<endl;
        abort();
    }
    string buff;
    bool store = false;
    int buff_pt[2];
    int buff_i = 0;
    vector<Point> vecpoint;
    vector<Line> vecline;
    // Reading input .ps file
    while (infile) {
        infile>>buff;
        if ((store == true) && (buff.compare("stroke") == 0)) {
            Curve buffc;
            buffc.set(vecpoint);
            veccurve.push_back(buffc);
            vecpoint.clear();
        }else if (buff.compare("%%%END") == 0) {
            store = false;
        }else if (buff.compare("%%%BEGIN") == 0) {
            store = true;
        }else if (store == true && (buff.compare("moveto")==0 || buff.compare("curveto")==0)) {
            if (buff_i == 0) {
                buff_i = 0;
            }else {
                cout<<"There must be 2 coordinates for one point. Please check your input file."<<endl;
                abort();
            }
        }else if (store == true && (buff.compare("Line") == 0)) {
            Line buffl;
            buffl.setp(vecpoint[0], vecpoint[1]);
            vecline.push_back(buffl);
            vecpoint.clear();
        }else {
            if (store == true) {
                buff_pt[buff_i] = atoi(buff.c_str());
                buff_i++;
                if (buff_i == 2) {
                    Point pt_buff;
                    pt_buff.set(buff_pt[0], buff_pt[1]);
                    vecpoint.push_back(pt_buff);
                    buff_i = 0;
                }
            }
        }
    }
    infile.close();
    // Do transformations before clip a line.
    int curvesize = (int)veccurve.size();
    int linesize = (int)vecline.size();
    if (s != 1.0) { // Do the scaling first.
        for (int i = 0; i < curvesize; i++) {
            veccurve[i].scale(s,s);
        }
        for (int i = 0; i < linesize; i++) {
            vecline[i].scale(s, s);
        }
    }
    if (r != 0) { // Do rotation second.
        for (int i = 0; i < curvesize; i++) {
            veccurve[i].rot(r);
        }
        for (int i = 0; i < linesize; i++) {
            vecline[i].rot(r);
        }
    }
    if (m != 0 || n != 0) { // Do translation last.
        for (int i = 0; i < curvesize; i++) {
            veccurve[i].trans(m, n);
        }
        for (int i = 0; i < linesize; i++) {
            vecline[i].trans(m, n);
        }
    }
    for (int i = 0; i < linesize; i++) {
        vecline[i] = csclip(vecline[i], a, b, c, d);
    }
    for (int i = 0; i < curvesize; i++) {
        veccurve[i].toviewport();
    }
    for (int i = 0; i < linesize; i++) {
        vecline[i].toviewport();
    }
    for (int i = 0; i < linesize; i++) { // Write pixels for all lines in the world window.
        if (vecline[i].getin() == true) {
            vecline[i].showall();
        }
    }
    for (int i = 0; i < curvesize; i++) { // Draw curves by de Casteljau Algorithm. 
        veccurve[i] = judgecurves(veccurve[i]);
        if (veccurve[i].getin() == true) {
            veccurve[i].draw();
        }
    }
    // Prepare to wirte output file.
    output = "/Users/wushangqi/out.xpm"; // dir only for debug ******************************************
    outfile(output);
    string shell = "display " + output;
    system(shell.c_str());
    return 0;
}

//-------------------------------------------------
int cstest(Line argl, int a, int b, int c, int d) { // Return 1 or 2 or 3 if the line is completely visible. The simple test of C-S algorithm
    int xs = argl.getxmin();
    int ys = argl.getxminy();
    int xe = argl.getxmax();
    int ye = argl.getxmaxy();
    if (((b<=ys&&ys<=d) && (b<=ye&&ye<=d)) && ((xs<a) && (xe>c))) {
        return 1; //Lines go from WL to WR.
    }else if ((a<=xs&&xs<=c) && (a<=xe&&xe<=c) && (ys<b && ye>d)) {
        return 2; //Lines go from WB to WT.
    }else if ((a<=xs&&xs<=c) && (a<=xe&&xe<=c) && (b<=ys&&ys<=d) && (b<=ye&&ye<=d)) {
        return 3; //Lines begin and end within the world window.
    }else return 0; //Lines cannot pass the simple exam.
}

//--------------------------------------------------
Line csclip(Line argl, int a, int b, int c, int d) { //Futher test and clip the line by using Cohen-Sutherland algorithm.
    int xleft = argl.getxmin();
    int xlefty = argl.getxminy();
    int xright = argl.getxmax();
    int xrighty = argl.getxmaxy();
    int ya = argl.cal_inter_wlr((float)a).gety(); // intersaction with WL(x=a)
    int yc = argl.cal_inter_wlr((float)c).gety(); // intersaction with WR(x=c)
    int xb = argl.cal_inter_wtb((float)b).getx(); // intersaction with WB(y=b)
    int xd = argl.cal_inter_wtb((float)d).getx(); // intersaction with WT(y=d)
    int flag = cstest(argl, a, b, c, d);
    if (flag == 0) { //If the line has not passed the simple test, take it to do the following complex test.
        if (xleft < a && xright > c) {
            if ((b<=ya&&ya<=d) && (b<=yc&&yc<=d)) {
                argl.setall(a, ya, c, yc);
                argl.setin(true);
            }else if ((b<=ya&&ya<=d) && yc<b){
                argl.setall(a, ya, xb, b);
                argl.setin(true);
            }else if ((b<=ya&&ya<=d) && yc>d) {
                argl.setall(a, ya, xd, d);
                argl.setin(true);
            }else if (ya<b && (b<=yc&&yc<=d)) {
                argl.setall(xb, b, c, yc);
                argl.setin(true);
            }else if (ya>d && (b<=yc&yc<=d)) {
                argl.setall(xd, d, c, yc);
                argl.setin(true);
            }else if (ya<b && yc>d) {
                argl.setall(xb, b, xd, d);
                argl.setin(true);
            }else if (ya>d && yc<b) {
                argl.setall(xd, d, xb, b);
                argl.setin(true);
            }else argl.setin(false);
        }else if ((a<=xleft && xleft<=c) && xright>c) {
            if ((b<=xlefty&&xlefty<=d) && (b<=yc&&yc<=d)) {
                argl.setall(xleft, xlefty, c, yc);
                argl.setin(true);
            }else if ((b<=xlefty&xlefty<=d) && yc<b) {
                argl.setall(xleft, xlefty, xb, b);
                argl.setin(true);
            }else if ((b<=xlefty&&xlefty<=d) && yc>d) {
                argl.setall(xleft, xlefty, xd, d);
                argl.setin(true);
            }else if (xlefty<b && (b<=yc&&yc<=d)) {
                argl.setall(xb, b, c, yc);
                argl.setin(true);
            }else if (xlefty<b && yc>d) {
                argl.setall(xb, b, xd, d);
                argl.setin(true);
            }else if (xlefty>d && (b<=yc&&yc<=d)) {
                argl.setall(xd, d, c, yc);
                argl.setin(true);
            }else if (xlefty>d && yc<b) {
                argl.setall(xd, d, xb, b);
                argl.setin(true);
            }else argl.setin(false);
        }else if (xleft<a && (a<=xright && xright<=c)) {
            if ((b<=ya&&ya<=d) && (b<=xrighty&&xrighty<=d)) {
                argl.setall(a, ya, xright, xrighty);
                argl.setin(true);
            }else if (ya<b && (b<=xrighty&&xrighty<=d)) {
                argl.setall(xb, b, xright, xrighty);
                argl.setin(true);
            }else if (ya>d && (b<=xrighty&&xrighty<=d)) {
                argl.setall(xd, d, xright, xrighty);
                argl.setin(true);
            }else if ((b<=ya&&ya<=d) && xrighty<b) {
                argl.setall(a, ya, xb, b);
                argl.setin(true);
            }else if ((b<=ya&&ya<=d) && xrighty>d) {
                argl.setall(a, ya, xd, d);
                argl.setin(true);
            }else if (ya<b && xrighty>d) {
                argl.setall(xb, b, xd, d);
                argl.setin(true);
            }else if (ya>d &&xrighty<b) {
                argl.setall(xd, d, xb, b);
                argl.setin(true);
            }else argl.setin(false);
        }else if ((a<=xleft&&xleft<=c) && (a<=xright&&xright<=c)) {
            if (xlefty<b && (b<=xrighty&&xrighty<=d)) {
                argl.setall(xb, b, xright, xrighty);
                argl.setin(true);
            }else if (xlefty<b && xrighty>d) {
                argl.setall(xb, b, xd, d);
                argl.setin(true);
            }else if (xlefty>d && xrighty<b) {
                argl.setall(xd, d, xb, b);
                argl.setin(true);
            }else if (xlefty>d && (b<=xrighty&&xrighty<=d)) {
                argl.setall(xd, d, xright, xrighty);
                argl.setin(true);
            }else if ((b<=xlefty&&xlefty<=d) && xrighty<b) {
                argl.setall(xleft, xlefty, xb, b);
                argl.setin(true);
            }else if ((b<=xlefty&&xlefty<=d) && xrighty>d) {
                argl.setall(xleft, xlefty, xd, d);
                argl.setin(true);
            }else if ((b<=xlefty&&xlefty<=d) && (b<=xrighty&&xrighty<=d)) {
                argl.setin(true);
            }else argl.setin(false);
        }else argl.setin(false);
    }else if (flag == 1) {
        argl.setall(a, ya, c, yc);
        argl.setin(true);
    }else if (flag == 2) {
        argl.setall(xb, b, xd, d);
        argl.setin(true);
    }else if (flag == 3) { // Lines totally within the window do not need to be clipped.
        argl.setin(true);
    }else argl.setin(false); // Mark and ignore all lines out the world window.
    return argl;
}

Curve judgecurves(Curve cur){
    Point buff;
    int count = 0;
    for (int i = 0; i < 4; i++) {
        buff = cur.getp(i);
        int x = buff.getx();
        int y = buff.gety();
        if (((x<a) || (x>c)) && ((y<b) || (d<y))) {
            count++;
        }
    }
    if (count == 4) {
        cur.setin(false);
    }else cur.setin(true);
    return cur;
}

//-------------------------------------------------------------------------------------------
void outfile(string output){
    ofstream out(output.c_str());
    string line = "";
    if (!out) {
        cout<<"Cannot write an output file, please check your output path."<<endl;
    }
    out<<setheader()<<endl;
    int height = 501;
    int width = 501;
    cout<<setheader()<<endl;
    for (int i = height - 1; i >= 0; i--) {
        for (int j = 0; j< width; j++) {
            line += xpmpix[i][j];
        }
        line = "\"" + line + "\"";
        if (i != 0) {
            line = line + ",";
        }
        out<<line<<endl;
        cout<<line<<endl;
        line.clear();
    }
    out<<setend()<<endl;
    cout<<setend()<<endl;
    out.close();
    
}

//-----------------------------------------------------------------------------
int rnd(float arg){ //return a rounded vaule of a float
    if (arg >= 0) {
        return (int)(arg + 0.5);
    }else {
        return (int)(arg - 0.5);
    }
}

//----------------------------------------------------------------------------------------
void help(){
    cout<<"[-f] The next argument is the input \"Postscript\" file."<<endl;
    cout<<"[-s] This next argument is a float specifying the scaling factor in both dimensions about the world origin."<<endl;
    cout<<"[-r] This next argument is an integer specifying the number of degrees for a counter-clockwise rotation about the world origin."<<endl;
    cout<<"[-m] The next argument is an integer specifying a translation in the x dimension."<<endl;
    cout<<"[-n] The next argument is an integer specifying a translation in the y dimension."<<endl;
    cout<<"[-a] The next argument is an integer lower bound in the x dimension of the world window."<<endl;
    cout<<"[-b] The next argument is an integer lower bound in the y dimension of the world window."<<endl;
    cout<<"[-c] The next argument is an integer upper bound in the x dimension of the world window."<<endl;
    cout<<"[-d] The next argument is an integer upper bound in the y dimension of the world window."<<endl;
    cout<<"[-j] The next argument is an integer lower bound in the x dimension of the viewport window."<<endl;
    cout<<"[-k] The next argument is an integer lower bound in the y dimension of the viewport window."<<endl;
    cout<<"[-o] The next argument is an integer upper bound in the x dimension of the viewport window."<<endl;
    cout<<"[-p] The next argument is an integer upper bound in the y dimension of the viewport window."<<endl;
    cout<<"[-L] This option specifies the increment of the parameter u to be used when generating points on the curve. The next argument is a real number between 0 and 1."<<endl;
    cout<<"This program will generate ./out.xpm file automatically."<<endl;
}

//------------------------------------------------------------------------------
void optana(int argc, char * const argv[]){
    // analyze input option and set defualt options
    float temp;
    int opt;
    while ((opt = getopt(argc, argv, "f:a:b:c:d:s:m:n:r:j:k:o:p:L:h"))!= -1) {
        switch (opt) {
            case 'f':{
                input = optarg;
                break;
            }
            case 'a':{
                string astr(optarg);
                temp = atof(astr.c_str());
                a = rnd(temp);
                break;
            }
            case 'b':{
                string bstr(optarg);
                temp = atof(bstr.c_str());
                b = rnd(temp);
                break;
            }
            case 'c':{
                string cstr(optarg);
                temp = atof(cstr.c_str());
                c = rnd(temp);
                break;
            }
            case 'd':{
                string dstr(optarg);
                temp = atof(dstr.c_str());
                d = rnd(temp);
                break;
            }
            case 'm':{
                string mstr(optarg);
                temp = atof(mstr.c_str());
                m = rnd(temp);
                break;
            }
            case 'n':{
                string nstr(optarg);
                temp = atof(nstr.c_str());
                n = rnd(temp);
                break;
            }
            case 'r':{
                string rstr(optarg);
                temp = atof(rstr.c_str());
                r = rnd(temp);
                break;
            }
            case 's':{
                string sstr(optarg);
                s = atof(sstr.c_str());
                break;
            }
            case 'j':{
                string jstr(optarg);
                temp = atof(jstr.c_str());
                j = rnd(temp);
                break;
            }
            case 'k':{
                string kstr(optarg);
                temp = atof(kstr.c_str());
                k = rnd(temp);
                break;
            }
            case 'o':{
                string ostr(optarg);
                temp = atof(ostr.c_str());
                o = rnd(temp);
                break;
            }
            case 'p':{
                string pstr(optarg);
                temp = atof(pstr.c_str());
                p = rnd(temp);
                break;
            }
            case 'L':{
                string Lstr(optarg);
                L = atof(Lstr.c_str());
                break;
            }
            case 'h':help();abort();
                break;
            default:cout<<"Your input is not correct, please enter -h for help."<<endl;abort();
                break;
        }
    }
    if (a > c) { // switch if a is bigger than c
        int tmp = c;
        c = a;
        a = tmp;
    }
    if (b > d) {
        int tmp = d;
        d = b;
        b = tmp;
    }
    if (a == c || b == d) {
        cout<<"Lower bound and upper bound of the world window cannot be the same."<<endl;
        abort();
    }
    if (j > o) {
        int tmp = j;
        j = o;
        o = tmp;
    }
    if (k > p) {
        int tmp = k;
        k = p;
        p = tmp;
    }
    if (j ==o || k ==p) {
        cout<<"Lower bound and upper bound of the viewport cannot be the same."<<endl;
        abort();
    }
    if (j<0 || k<0 || 500<o || 500<p) {
        cout<<"Viewport must within the picture."<<endl;
        abort();
    }
    if (L>1 || L<0) {
        cout<<"The step of parameter cannot exceeds 1."<<endl;
        abort();
    }
}


//-------------------------------------------------------------------------
string setheader() {
    stringstream tmp;
    int w = 501;
    int h = 501;
    tmp<<w;
    string intw;
    tmp>>intw;
    tmp.clear();
    tmp<<h;
    string inth;
    tmp>>inth;
    string str = "/* XPM */\nstatic char *CG_hw6[] = {\n/* width height num_colors chars_per_pixel */\n\"" + intw + " " + inth +" 2 1\",\n/* colors */\n\"- c #ffffff\",\n\"+ c #000000\",\n/* pixels */";
    return str;
}

//-------------------------------------------------------------------------
string setend() {
    string str = "};";
    return str;
}