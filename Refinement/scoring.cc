//#include "options.h"
#include "scoring.h"
//#include "utils.h"
//#include "constants.h"

//#include <iomanip>
#include <ANN/ANN.h>

#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)


const double EFACTOR = 332.0637;
const double probe_radius = 1.5;
const int number_of_dots = 500;
const double Z = 2.0 * M_PI * sqrt(M_PI);

/*********** SASA CALCULATION ROUTINES **********************/
// TAKEN from the BALL C++ library

typedef double * point_double;
typedef int    * point_int;
#define FOURPI (4.*M_PI)
#define TORAD(A)     ((A)*0.017453293)
#define DP_TOL     0.001
#define MAXIMUM(A, B)  (((A) > (B) ? (A) : (B)))
#define MINIMUM(A, B)  (((A) < (B) ? (A) : (B)))
#define UNSP_ICO_DOD      9
#define UNSP_ICO_ARC     10
#define CALLOC(n, size) mycalloc(__FILE__,__LINE__, n, size)
#define REALLOC(ptr, size) myrealloc(__FILE__,__LINE__,  ptr, size)

point_double xpunsp=NULL;
point_int ico_wk=NULL, ico_pt=NULL;
double del_cube;
int n_dot, ico_cube, last_n_dot = 0, last_densit = 0, last_unsp=0;
int last_cubus = 0;
double rg, rh;

void* mycalloc(const char * filename, const int linenr, size_t nelem, size_t elsize)
{
    int * ip;
    ip = (int *) calloc(nelem, elsize);
    if(ip == NULL)
        cerr << "calculateSASAreaCALLOC : failed in file " << filename
        << " at line " << linenr << endl;
    return(ip);
}


void * myrealloc(const char * filename, const int linenr, void * ptr, size_t size)
{
    int * ip;
    ip = (int *) realloc(ptr, size);
    if(ip == NULL)
        cerr << "REALLOC : failed in file " << filename << " at line " << linenr << endl;
    return(ip);
}

double asin_safe(double f)
{
    if((fabs(f) < 1.00))
        return( asin(f) );
    if((fabs(f) - 1.00)  <= DP_TOL)
        cerr << "calculateSASArea: invalid argument" << f << endl;
    return(M_PI * M_PI);
}

void icosaeder_vertices(double *xus)
{
    rh = sqrt(1.-2.*cos(TORAD(72.)))/(1.-cos(TORAD(72.)));
    rg = cos(TORAD(72.))/(1.-cos(TORAD(72.)));
    //icosaeder vertices
    xus[ 0] = 0.;                  xus[ 1] = 0.;                  xus[ 2] = 1.;
    xus[ 3] = rh*cos(TORAD(72.));  xus[ 4] = rh*sin(TORAD(72.));  xus[ 5] = rg;
    xus[ 6] = rh*cos(TORAD(144.)); xus[ 7] = rh*sin(TORAD(144.)); xus[ 8] = rg;
    xus[ 9] = rh*cos(TORAD(216.)); xus[10] = rh*sin(TORAD(216.)); xus[11] = rg;
    xus[12] = rh*cos(TORAD(288.)); xus[13] = rh*sin(TORAD(288.)); xus[14] = rg;
    xus[15] = rh;                  xus[16] = 0;                   xus[17] = rg;
    xus[18] = rh*cos(TORAD(36.));  xus[19] = rh*sin(TORAD(36.));  xus[20] = -rg;
    xus[21] = rh*cos(TORAD(108.)); xus[22] = rh*sin(TORAD(108.)); xus[23] = -rg;
    xus[24] = -rh;                 xus[25] = 0;                   xus[26] = -rg;
    xus[27] = rh*cos(TORAD(252.)); xus[28] = rh*sin(TORAD(252.)); xus[29] = -rg;
    xus[30] = rh*cos(TORAD(324.)); xus[31] = rh*sin(TORAD(324.)); xus[32] = -rg;
    xus[33] = 0.;                  xus[34] = 0.;                  xus[35] = -1.;
}


void divarc(double x1, double y1, double z1, double x2, double y2, double z2,
            int div1, int div2, double *xr, double *yr, double *zr)
{
    double xd, yd, zd, dd, d1, d2, s, x, y, z;
    double phi, sphi, cphi;
    
    xd = y1*z2-y2*z1;
    yd = z1*x2-z2*x1;
    zd = x1*y2-x2*y1;
    dd = sqrt(xd*xd+yd*yd+zd*zd);
    if(dd < DP_TOL)
        cerr << "divarc: rotation axis of length " << dd << ::std::endl;
    
    d1 = x1*x1+y1*y1+z1*z1;
    if(d1 < 0.5)
        cerr << "divarc: vector 1 of sq.length " << d1 << endl;
    
    d2 = x2*x2+y2*y2+z2*z2;
    if(d2 < 0.5)
        cerr << "divarc: vector 2 of sq.length " << d2 << endl;
    
    phi = asin_safe(dd/sqrt(d1*d2));
    phi = phi*((double)div1)/((double)div2);
    sphi = sin(phi); cphi = cos(phi);
    s  = (x1*xd+y1*yd+z1*zd)/dd;
    
    x = xd*s*(1.-cphi)/dd + x1 * cphi + (yd*z1-y1*zd)*sphi/dd;
    y = yd*s*(1.-cphi)/dd + y1 * cphi + (zd*x1-z1*xd)*sphi/dd;
    z = zd*s*(1.-cphi)/dd + z1 * cphi + (xd*y1-x1*yd)*sphi/dd;
    dd = sqrt(x*x+y*y+z*z);
    *xr = x/dd; *yr = y/dd; *zr = z/dd;
}


/**
 * densit...required dots per unit sphere
 * dot distribution on a unit sphere based on an icosaeder
 * great circle average refining of icosahedral face
 */
int ico_dot_arc(int densit)
{ 
    int i, j, k, tl, tl2, tn, tess;
    double a, d, x, y, z, x2, y2, z2, x3, y3, z3;
    double xij, yij, zij, xji, yji, zji, xik, yik, zik, xki, yki, zki,
    xjk, yjk, zjk, xkj, ykj, zkj;
    point_double xus=NULL;
    
    //calculate tessalation level
    a = sqrt((((double) densit)-2.)/10.);
    tess = (int) ceil(a);
    n_dot = 10*tess*tess+2;
    if(n_dot < densit)
        cerr << "ico_dot_arc: error in formula for tessalation level (" << tess
        << "->" << n_dot << ", " << densit << ")" << endl;
    
    xus = (double *) CALLOC(3*n_dot, sizeof(double));
    xpunsp = xus;
    icosaeder_vertices(xus);
    
    if(tess > 1)
    {
        tn = 12;
        a = rh*rh*2.*(1.-cos(TORAD(72.)));
        //calculate tessalation of icosaeder edges
        for(i = 0; i < 11; i++)
        {
            for(j=i+1; j < 12; j++)
            {
                x = xus[3*i]-xus[3*j];
                y = xus[1+3*i]-xus[1+3*j];
                z = xus[2+3*i]-xus[2+3*j];
                d = x*x+y*y+z*z;
                if(fabs(a-d) > DP_TOL)
                    continue;
                for(tl=1; tl<tess; tl++)
                {
                    if(tn >= n_dot)
                        cerr << "ico_dot: tn exceeds dimension of xus" << ::std::endl;
                    divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                           xus[3*j], xus[1+3*j], xus[2+3*j],
                           tl, tess, &xus[3*tn], &xus[1+3*tn], &xus[2+3*tn]);
                    tn++;
                }
            }
        }
        // calculate tessalation of icosaeder faces
        for(i = 0; i < 10; i++)
        {
            for(j=i+1; j < 11; j++)
            {
                x = xus[3*i]-xus[3*j];
                y = xus[1+3*i]-xus[1+3*j];
                z = xus[2+3*i]-xus[2+3*j];
                d = x*x+y*y+z*z;
                if(fabs(a-d) > DP_TOL)
                    continue;
                
                for(k=j+1; k < 12; k++)
                {
                    x = xus[3*i]-xus[3*k];
                    y = xus[1+3*i]-xus[1+3*k];
                    z = xus[2+3*i]-xus[2+3*k];
                    d = x*x+y*y+z*z;
                    if(fabs(a-d) > DP_TOL)
                        continue;
                    x = xus[3*j]-xus[3*k];
                    y = xus[1+3*j]-xus[1+3*k];
                    z = xus[2+3*j]-xus[2+3*k];
                    d = x*x+y*y+z*z;
                    if(fabs(a-d) > DP_TOL)
                        continue;
                    for(tl=1; tl<tess-1; tl++)
                    {
                        divarc(xus[3*j], xus[1+3*j], xus[2+3*j],
                               xus[3*i], xus[1+3*i], xus[2+3*i],
                               tl, tess, &xji, &yji, &zji);
                        divarc(xus[3*k], xus[1+3*k], xus[2+3*k],
                               xus[3*i], xus[1+3*i], xus[2+3*i],
                               tl, tess, &xki, &yki, &zki);
                        
                        for(tl2=1; tl2<tess-tl; tl2++)
                        {
                            divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                                   xus[3*j], xus[1+3*j], xus[2+3*j],
                                   tl2, tess, &xij, &yij, &zij);
                            divarc(xus[3*k], xus[1+3*k], xus[2+3*k],
                                   xus[3*j], xus[1+3*j], xus[2+3*j],
                                   tl2, tess, &xkj, &ykj, &zkj);
                            divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                                   xus[3*k], xus[1+3*k], xus[2+3*k],
                                   tess-tl-tl2, tess, &xik, &yik, &zik);
                            divarc(xus[3*j], xus[1+3*j], xus[2+3*j],
                                   xus[3*k], xus[1+3*k], xus[2+3*k],
                                   tess-tl-tl2, tess, &xjk, &yjk, &zjk);
                            if(tn >= n_dot)
                                cerr << "ico_dot: tn exceeds dimension of xus" << endl;
                            divarc(xki, yki, zki, xji, yji, zji, tl2, tess-tl, &x, &y, &z);
                            divarc(xkj, ykj, zkj, xij, yij, zij, tl, tess-tl2, &x2, &y2, &z2);
                            divarc(xjk, yjk, zjk, xik, yik, zik, tl, tl+tl2, &x3, &y3, &z3);
                            x = x+x2+x3; y = y+y2+y3; z = z+z2+z3;
                            d = sqrt(x*x+y*y+z*z);
                            xus[3*tn] = x/d;
                            xus[1+3*tn] = y/d;
                            xus[2+3*tn] = z/d;
                            tn++;
                        }// cycle tl2
                    }//cycle tl
                }//cycle k
            }// cycle j
        }//cycle i
        if(n_dot != tn)
            cerr << "ico_dot: n_dot(" << n_dot << ") and tn(" << tn << ") differ" << endl;
    }//end of if (tess > 1)
    return n_dot;
}// end of routine ico_dot_arc

/**
 * densit...required dots per unit sphere
 * dot distribution on a unit sphere based on an icosaeder
 * great circle average refining of icosahedral face
 */
int ico_dot_dod(int densit)
{
    int i, j, k, tl, tl2, tn, tess, j1, j2;
    double a, d, x, y, z, x2, y2, z2, x3, y3, z3, ai_d, adod;
    double xij, yij, zij, xji, yji, zji, xik, yik, zik, xki, yki, zki,
    xjk, yjk, zjk, xkj, ykj, zkj;
    point_double xus=NULL;
    // calculate tesselation level
    a = sqrt((((double) densit)-2.)/30.);
    tess = MAXIMUM((int) ceil(a), 1);
    n_dot = 30*tess*tess+2;
    if(n_dot < densit)
        cerr << "ico_dot_dod: error in formula for tessalation level (" << tess << "->"
        << n_dot << ", " << densit << ")" << endl;
    
    xus = (double *) CALLOC(3*n_dot, sizeof(double));
    xpunsp = xus;
    icosaeder_vertices(xus);
    
    tn=12;
    //square of the edge of an icosaeder */
    a = rh*rh*2.*(1.-cos(TORAD(72.)));
    
    //dodecaeder vertices
    for(i = 0; i < 10; i++)
    {
        for(j=i+1; j < 11; j++)
        {
            x = xus[3*i]-xus[3*j];
            y = xus[1+3*i]-xus[1+3*j];
            z = xus[2+3*i]-xus[2+3*j];
            d = x*x+y*y+z*z;
            if(fabs(a-d) > DP_TOL)
                continue;
            for (k=j+1; k < 12; k++)
            {
                x = xus[3*i]-xus[3*k];
                y = xus[1+3*i]-xus[1+3*k];
                z = xus[2+3*i]-xus[2+3*k];
                d = x*x+y*y+z*z;
                if(fabs(a-d) > DP_TOL)
                    continue;
                x = xus[3*j]-xus[3*k];
                y = xus[1+3*j]-xus[1+3*k];
                z = xus[2+3*j]-xus[2+3*k];
                d = x*x+y*y+z*z;
                if(fabs(a-d) > DP_TOL)
                    continue;
                x = xus[  3*i]+xus[  3*j]+xus[  3*k];
                y = xus[1+3*i]+xus[1+3*j]+xus[1+3*k];
                z = xus[2+3*i]+xus[2+3*j]+xus[2+3*k];
                d = sqrt(x*x+y*y+z*z);
                xus[3*tn]=x/d; xus[1+3*tn]=y/d; xus[2+3*tn]=z/d;
                tn++;
            }
        }
    }
    
    if(tess > 1)
    {
        tn = 32;
        //square of the edge of an dodecaeder
        adod = 4.*(cos(TORAD(108.))-cos(TORAD(120.)))/(1.-cos(TORAD(120.)));
        //square of the distance of two adjacent vertices of ico- and dodecaeder
        ai_d = 2.*(1.-sqrt(1.-a/3.));
        
        //calculate tessalation of mixed edges
        for(i = 0; i < 31; i++)
        {
            j1 = 12; j2 = 32; a = ai_d;
            if(i>=12)
            { j1=i+1; a = adod;}
            for(j=j1; j < j2; j++)
            {
                x = xus[3*i]-xus[3*j];
                y = xus[1+3*i]-xus[1+3*j];
                z = xus[2+3*i]-xus[2+3*j];
                d = x*x+y*y+z*z;
                if(fabs(a-d) > DP_TOL)
                    continue;
                for(tl=1; tl<tess; tl++)
                {
                    if(tn >= n_dot)
                        cerr << "ico_dot: tn exceeds dimension of xus" << endl;
                    divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                           xus[3*j], xus[1+3*j], xus[2+3*j],
                           tl, tess, &xus[3*tn], &xus[1+3*tn], &xus[2+3*tn]);
                    tn++;
                }
            }
        }
        //calculate tessalation of pentakisdodecahedron faces
        for(i = 0; i < 12; i++)
        {
            for(j=12; j < 31; j++)
            {
                x = xus[3*i]-xus[3*j];
                y = xus[1+3*i]-xus[1+3*j];
                z = xus[2+3*i]-xus[2+3*j];
                d = x*x+y*y+z*z;
                if(fabs(ai_d-d) > DP_TOL)
                    continue;
                
                for(k=j+1; k < 32; k++)
                {
                    x = xus[3*i]-xus[3*k];
                    y = xus[1+3*i]-xus[1+3*k];
                    z = xus[2+3*i]-xus[2+3*k];
                    d = x*x+y*y+z*z;
                    if(fabs(ai_d-d) > DP_TOL)
                        continue;
                    x = xus[3*j]-xus[3*k];
                    y = xus[1+3*j]-xus[1+3*k];
                    z = xus[2+3*j]-xus[2+3*k];
                    d = x*x+y*y+z*z;
                    if(fabs(adod-d) > DP_TOL)
                        continue;
                    for(tl=1; tl<tess-1; tl++)
                    {
                        divarc(xus[3*j], xus[1+3*j], xus[2+3*j],
                               xus[3*i], xus[1+3*i], xus[2+3*i],
                               tl, tess, &xji, &yji, &zji);
                        divarc(xus[3*k], xus[1+3*k], xus[2+3*k],
                               xus[3*i], xus[1+3*i], xus[2+3*i],
                               tl, tess, &xki, &yki, &zki);
                        
                        for(tl2=1; tl2<tess-tl; tl2++)
                        {
                            divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                                   xus[3*j], xus[1+3*j], xus[2+3*j],
                                   tl2, tess, &xij, &yij, &zij);
                            divarc(xus[3*k], xus[1+3*k], xus[2+3*k],
                                   xus[3*j], xus[1+3*j], xus[2+3*j],
                                   tl2, tess, &xkj, &ykj, &zkj);
                            divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                                   xus[3*k], xus[1+3*k], xus[2+3*k],
                                   tess-tl-tl2, tess, &xik, &yik, &zik);
                            divarc(xus[3*j], xus[1+3*j], xus[2+3*j],
                                   xus[3*k], xus[1+3*k], xus[2+3*k],
                                   tess-tl-tl2, tess, &xjk, &yjk, &zjk);
                            if(tn >= n_dot)
                                cerr << "ico_dot: tn exceeds dimension of xus" << endl;
                            divarc(xki, yki, zki, xji, yji, zji, tl2, tess-tl, &x, &y, &z);
                            divarc(xkj, ykj, zkj, xij, yij, zij, tl, tess-tl2, &x2, &y2, &z2);
                            divarc(xjk, yjk, zjk, xik, yik, zik, tl, tl+tl2, &x3, &y3, &z3);
                            x = x+x2+x3;
                            y = y+y2+y3;
                            z = z+z2+z3;
                            d = sqrt(x*x+y*y+z*z);
                            xus[3*tn] = x/d;
                            xus[1+3*tn] = y/d;
                            xus[2+3*tn] = z/d;
                            tn++;
                        }//cycle tl2
                    }//cycle tl
                }//cycle k
            }//cycle j
        }//cycle i
        //if(n_dot != tn)
        //    cerr << "ico_dot: n_dot(" << n_dot << ") and tn(" << tn << ") differ" << endl;
    }//end of if (tess > 1)
    return n_dot;
} // end of routine ico_dot_dod

int unsp_type(int densit)
{
    int i1, i2;
    i1 = 1;
    while(10*i1*i1+2 < densit)
        i1++;
    i2 = 1;
    while (30*i2*i2+2 < densit)
        i2++;
    if(10*i1*i1-2 < 30*i2*i2-2)
        return UNSP_ICO_ARC;
    else
        return UNSP_ICO_DOD;
}

int make_unsp(int densit, int mode, int * num_dot, int cubus)
{
    int ndot, ico_cube_cb, i, j, k, l, ijk, tn, tl, tl2;
    point_double xus;
    point_int    work;
    double x, y, z;
    
    if(xpunsp)
        free(xpunsp);
    if(ico_wk)
        free(ico_wk);
    
    k=1;
    if(mode < 0)
    {
        k=0;
        mode = -mode;
    }
    if(mode == UNSP_ICO_ARC)
        ndot = ico_dot_arc(densit);
    else if(mode == UNSP_ICO_DOD)
        ndot = ico_dot_dod(densit);
    else
    {
        cerr << "make_unsp: mode " << ((k)?'+':'-') << (int)mode << " not allowed" << endl;
        return 1;
    }
    
    last_n_dot = ndot;
    last_densit = densit;
    last_unsp = mode;
    *num_dot=ndot;
    if (k) return 0;
    
    // in the following the dots of the unit sphere may be resorted
    last_unsp = -last_unsp;
    
    // determine distribution of points in elementary cubes
    if(cubus)
    {
        ico_cube = cubus;
    }
    else
    {
        last_cubus = 0;
        i=1;
        while(i*i*i*2 < ndot)
            i++;
        ico_cube = MAXIMUM(i-1, 0);
    }
    ico_cube_cb = ico_cube*ico_cube*ico_cube;
    del_cube=2./((double)ico_cube);
    work = (int *) CALLOC(ndot, sizeof(int));
    xus = xpunsp;
    for(l=0; l<ndot; l++)
    {
        i = MAXIMUM((int) floor((1.+xus[3*l])/del_cube), 0);
        if(i>=ico_cube)
            i = ico_cube-1;
        j = MAXIMUM((int) floor((1.+xus[1+3*l])/del_cube), 0);
        if(j>=ico_cube)
            j = ico_cube-1;
        k = MAXIMUM((int) floor((1.+xus[2+3*l])/del_cube), 0);
        if(k>=ico_cube)
            k = ico_cube-1;
        ijk = i+j*ico_cube+k*ico_cube*ico_cube;
        work[l] = ijk;
    }
    
    ico_wk = (int *) CALLOC(2*ico_cube_cb+1, sizeof(int));
    ico_pt = ico_wk+ico_cube_cb;
    for(l=0; l<ndot; l++)
        ico_wk[work[l]]++;   // dots per elementary cube
    
    
    // reordering of the coordinate array in accordance with box number
    tn=0;
    for(i = 0; i < ico_cube; i++)
    {
        for(j = 0; j < ico_cube; j++)
        {
            for(k=0; k < ico_cube; k++)
            {
                tl=0;
                tl2 = tn;
                ijk = i+ico_cube*j+ico_cube*ico_cube*k;
                *(ico_pt+ijk) = tn;
                for(l=tl2; l<ndot; l++)
                {
                    if(ijk == work[l])
                    {
                        x = xus[3*l];
                        y = xus[1+3*l];
                        z = xus[2+3*l];
                        xus[3*l] = xus[3*tn];
                        xus[1+3*l] = xus[1+3*tn];
                        xus[2+3*l] = xus[2+3*tn];
                        xus[3*tn] = x;
                        xus[1+3*tn] = y;
                        xus[2+3*tn] = z;
                        ijk = work[l];
                        work[l]=work[tn];
                        work[tn]=ijk;
                        tn++;
                        tl++;
                    }
                }
                *(ico_wk+ijk) = tl;
            }// cycle k
        }// cycle j
    }// cycle i
    free(work);
    return 0;
}

int nsc_(double* co, double* radius, int nat, int densit, double* value_of_area,
         double** at_area, double** lidots, int* nu_dots,
         int** atom_dot_ptr)
{
    int iat, i, ii, iii, ix, iy, iz, ixe, ixs, iye, iys, ize, izs, i_ac;
    int jat, j, jj, jjj, jx, jy, jz;
    int distribution;
    int l;
    int maxnei, nnei, last;
    point_int wkdot=NULL, wkbox=NULL, wkat1=NULL, wkatm=NULL;
    Neighbour *wknb, *ctnb;
    int iii1, iii2, iiat, i_at, j_at;
    double dx, dy, dz, dd, ai, aisq, ajsq, aj, as, a;
    double xi, yi, zi, xs=0., ys=0., zs=0.;
    double dotarea, area;
    point_double xus, atom_area=NULL;
    
    
    int nxbox, nybox, nzbox, nxy, nxyz;
    double xmin, ymin, zmin, xmax, ymax, zmax, ra2max, d, *pco;
    
    distribution = unsp_type(densit);
    if(distribution != -last_unsp || last_cubus != 4 ||
       (densit != last_densit && densit != last_n_dot))
    {
        if(make_unsp(densit, (-distribution), &n_dot, 4))
            return 1;
    }
    xus = xpunsp;
    
    dotarea = FOURPI/(double) n_dot;
    area = 0.;
    
    //start with neighbour list
    //calculate neighbour list with the box algorithm
    if(nat == 0)
    {
        cerr << "nsc_dclm: no surface atoms selected" << endl;
        return 1;
    }
    atom_area = (double *) CALLOC(nat, sizeof(double));
    
    // dimensions of atomic set, cell edge is 2*ra_max
    xmin = co[0]; xmax = xmin; xs=xmin;
    ymin = co[1]; ymax = ymin; ys=ymin;
    zmin = co[2]; zmax = zmin; zs=zmin;
    ra2max = radius[0];
    
    for (iat=1; iat<nat; iat++)
    {
        pco = co+3*iat;
        xmin = MINIMUM(xmin, *pco);
        xmax = MAXIMUM(xmax, *pco);
        ymin = MINIMUM(ymin, *(pco+1));
        ymax = MAXIMUM(ymax, *(pco+1));
        zmin = MINIMUM(zmin, *(pco+2));
        zmax = MAXIMUM(zmax, *(pco+2));
        xs= xs+ *pco; ys = ys+ *(pco+1);
        zs= zs+ *(pco+2);
        ra2max = MAXIMUM(ra2max, radius[iat]);
    }
    xs = xs / (double) nat;
    ys = ys / (double) nat;
    zs = zs / (double) nat;
    ra2max = 2.*ra2max;
    d = xmax-xmin;
    nxbox = (int) MAXIMUM(ceil(d / ra2max), 1.);
    d = (((double)nxbox)*ra2max - d) / 2.;
    xmin = xmin - d;
    xmax = xmax + d;
    d = ymax - ymin;
    nybox = (int) MAXIMUM(ceil(d / ra2max), 1.);
    d = (((double)nybox)*ra2max - d) / 2.;
    ymin = ymin - d;
    ymax = ymax + d;
    d = zmax - zmin;
    nzbox = (int) MAXIMUM(ceil(d / ra2max), 1.);
    d = (((double)nzbox) * ra2max - d) / 2.;
    zmin = zmin - d;
    zmax = zmax + d;
    nxy = nxbox * nybox;
    nxyz = nxy * nzbox;
    //box number of atoms
    wkatm = (int*) CALLOC(3 * nat, sizeof(int));
    wkat1 = wkatm + nat;
    wkdot = (int*) CALLOC(n_dot + nxyz + 1, sizeof(int));
    wkbox = wkdot + n_dot;
    
    for(iat = 0; iat < nat; iat++)
    {
        pco = co + 3 * iat;
        i = (int) MAXIMUM(floor((  *pco  -xmin)/ra2max), 0);
        i = MINIMUM(i, nxbox -1);
        j = (int) MAXIMUM(floor((*(pco + 1) - ymin) / ra2max), 0);
        j = MINIMUM(j, nybox - 1);
        l = (int) MAXIMUM(floor((*(pco + 2) - zmin) / ra2max), 0);
        l = MINIMUM(l, nzbox - 1);
        i = i + j * nxbox + l * nxy;
        wkat1[iat] = i;
        wkbox[i]++;
    }
    
    //sorting of atoms in accordance with box numbers
    j = wkbox[0];
    for(i = 1; i < nxyz; i++)
        j = MAXIMUM(wkbox[i], j);
    for(i = 1; i <= nxyz; i++)
        wkbox[i] += wkbox[i-1];
    
    maxnei = MINIMUM(nat, 27*j);
    wknb = (Neighbour *) CALLOC(maxnei, sizeof(Neighbour));
    for(iat = 0; iat < nat; iat++)
        wkatm[--wkbox[wkat1[iat]]] = iat;
    // calculate surface for all atoms, step cube-wise
    for(iz=0; iz<nzbox; iz++)
    {
        iii = iz * nxy;
        izs = MAXIMUM(iz - 1, 0);
        ize = MINIMUM(iz + 2, nzbox);
        for(iy = 0; iy < nybox; iy++)
        {
            ii = iy * nxbox + iii;
            iys = MAXIMUM(iy - 1,0);
            iye = MINIMUM(iy + 2, nybox);
            for(ix=0; ix<nxbox; ix++)
            {
                i = ii + ix;
                iii1 = wkbox[i];
                iii2 = wkbox[i + 1];
                if(iii1 >= iii2)
                    continue;
                
                ixs = MAXIMUM(ix - 1, 0);
                ixe = MINIMUM(ix + 2, nxbox);
                
                iiat = 0;
                
                //make intermediate atom list */
                for(jz =izs; jz < ize; jz++)
                {
                    jjj = jz * nxy;
                    for(jy = iys; jy <iye; jy++)
                    {
                        jj = jy * nxbox + jjj;
                        for(jx=ixs; jx<ixe; jx++)
                        {
                            j = jj+jx;
                            for(jat = wkbox[j]; jat < wkbox[j + 1]; jat++)
                            {
                                wkat1[iiat] = wkatm[jat];
                                iiat++;
                            }
                        }
                    }
                }
                
                for(iat = iii1; iat < iii2; iat++)
                {
                    i_at = wkatm[iat];
                    ai = radius[i_at];
                    aisq = ai * ai;
                    pco = co + 3 * i_at;
                    xi = *pco;
                    yi = *(pco+1);
                    zi = *(pco+2);
                    for(i = 0; i < n_dot; i++)
                        *(wkdot+i)=0;
                    
                    ctnb = wknb;
                    nnei = 0;
                    for(j = 0; j < iiat; j++)
                    {
                        j_at = *(wkat1 + j);
                        if(j_at == i_at)
                            continue;
                        
                        aj = radius[j_at];
                        ajsq = aj*aj;
                        pco = co+3*j_at;
                        dx = *pco-xi;
                        dy = *(pco+1)-yi;
                        dz = *(pco+2)-zi;
                        dd = dx*dx+dy*dy+dz*dz;
                        
                        as = ai+aj;
                        if(dd > as*as)
                            continue;
                        nnei++;
                        ctnb->x = dx; ctnb->y = dy; ctnb->z = dz;
                        ctnb->dot = (dd+aisq-ajsq)/(2.*ai); /* reference dot product */
                        ctnb++;
                    }
                    
                    // check points on accessibility
                    if(nnei)
                    {
                        last = 0;
                        i_ac = 0;
                        for(l=0; l<n_dot; l++)
                        {
                            if(xus[3*l]*(wknb+last)->x + xus[1+3*l]*(wknb+last)->y
                               + xus[2+3*l]*(wknb+last)->z <= (wknb+last)->dot)
                            {
                                for(j = 0; j < nnei; j++)
                                {
                                    if(xus[3*l]*(wknb+j)->x + xus[1+3*l]*(wknb+j)->y +
                                       xus[2+3*l]*(wknb+j)->z > (wknb+j)->dot)
                                    {
                                        last = j;
                                        break;
                                    }
                                }
                                if(j >= nnei)
                                {
                                    i_ac++;
                                    wkdot[l] = 1;
                                }
                            }// end of cycle j
                        }//end of cycle l
                    }
                    else
                    {
                        i_ac  = n_dot;
                        for(l=0; l < n_dot; l++)
                            wkdot[l] = 1;
                    }
                    
                    a = aisq * dotarea * (double) i_ac;
                    area = area + a;
                    atom_area[i_at] = a;
                }//end of cycle "iat"
            }//end of cycle "ix"
        }// end of cycle "iy"
    }// end of cycle "iz"
    
    free(wkatm); free(wkdot); free(wknb);
    
    *at_area = atom_area;
    *value_of_area = area;
    
    return 0;
}


void get_atomic_sasa(vector<atom>& A)
{
    int N = 0;
    for(size_t i = 0;i<A.size();i++)
    {
        if(A[i].atype.at(0) == 'H')
            continue;
        if(A[i].arad > 0.0)
            N++;
    }
    
    double* coordinates = new double[N * 3];
    double* radii = new double[N];
    int j = 0;
    for(size_t i = 0;i<A.size();i++)
    {
        if(A[i].atype.at(0) == 'H')
            continue;
        if(A[i].arad > 0.0)
        {
            coordinates[j * 3]      = A[i].axyz[0];
            coordinates[j * 3 + 1]  = A[i].axyz[1];
            coordinates[j * 3 + 2]  = A[i].axyz[2];
            radii[j] = A[i].arad + probe_radius;
            j++;
        }
    }
    
    double area;
    int number_of_surface_dots;
    double* atom_areas = 0;
    double* surface_dots = 0;
    int* atom_dots = 0;
    
    // call nsc
    nsc_(coordinates, radii, N, number_of_dots, &area, &atom_areas,
         &surface_dots, &number_of_surface_dots, &atom_dots);
    
    
    j = 0;
    for(size_t i = 0; i<A.size();i++)
    {
        A[i].sas_area = 0.0;
        if(A[i].atype.at(0) == 'H')
            continue;
        if(A[i].arad > 0.0)
        {
            A[i].sas_area = atom_areas[j];
            j++;
        }
    }
    
    if(atom_areas != 0)
        free(atom_areas);
    if(surface_dots != 0)
        free(surface_dots);
    if(xpunsp != 0)
    {
        free(xpunsp);
        xpunsp = 0;
    }
    if(ico_wk != 0)
    {
        free(ico_wk);
        ico_wk = 0;
    }
    
    delete [] coordinates;
    delete [] radii;
}
/********************* SASA ENDS **************************/

void get_pdb_complex(vector<atom>& R, vector<atom>& L, vector<atom>& RL)
{
    RL = R;
    for(size_t i=0;i<L.size();i++)
    {
        RL.push_back(L[i]);
    }
}

void compute_acp(vector<atom>& R, vector<atom>& L, ANNkd_tree* T, double& acp)
{
    int l = 0, m = 0;
    ANNpoint f = annAllocPt(3);
    double eps = 0.;
    
    double r2 = 36.;
    
    for(size_t i=0;i<R.size();i++) // for the fixed protein
    {
        if(R[i].acp_type == 0) continue;
        if(R[i].INTF == 0)
            continue;
        f[0] = R[i].axyz[0];
        f[1] = R[i].axyz[1];
        f[2] = R[i].axyz[2];
        
        ANNidxArray nnIdx = NULL;
        ANNdistArray dists = NULL;
        
        int nnode = T->annkFRSearch(f, r2, 0, nnIdx, dists, eps); // search for n of them.
        if(nnode <= 0)
            continue;
        
        nnIdx = new ANNidx[nnode];
        dists = new ANNdist[nnode];
        T->annkFRSearch(f, r2, nnode, nnIdx, dists, eps); // search for n of them.
        l = R[i].acp_type - 1;
        for(int p1=0; p1<nnode; p1++)
        {
            int j = nnIdx[p1];
            if(L[j].INTF == 0)
                continue;
            if(L[j].acp_type == 0) continue;
            m = L[j].acp_type - 1;
            acp += ZHANG_ACP[l][m];
        }
        delete [] nnIdx;
        delete [] dists;
    }
    acp *= (1./21.); // converting to kcal/mol
}

// hydrogen bonding potential + disulphide
// need to make this orientation dependent
// The disulfide bond contribution is -1 energy units and is summed for all
// interface S-S pairs with distance between 1.9 and 2.1 Å.
// r is the distance between interface donor and acceptor atoms
// stuff taken from FIREDOCK
void compute_hbp(vector<atom>& R, vector<atom>& L, ANNkd_tree *T, double& hbp)
{
    // S-S energy
    double e_ss = 0.;
    hbp = 0.;
    
    double r = 0.;
    double R0 = 2.9;// the optimal distance for hydrogen bonding
    
    ANNpoint f = annAllocPt(3);
    double eps = 0.;
    
    double r2 = 16.0;
    
    for(size_t i=0;i<R.size();i++) // for the fixed protein
    {
        if(R[i].INTF == 0)
            continue;
        
        if(R[i].atype.at(0) == 'H') continue;
        f[0] = R[i].axyz[0];
        f[1] = R[i].axyz[1];
        f[2] = R[i].axyz[2];
        
        ANNidxArray nnIdx = NULL;
        ANNdistArray dists = NULL;
        
        int nnode = T->annkFRSearch(f, r2, 0, nnIdx, dists, eps); // search for n of them.
        if(nnode <= 0)
            continue;
        
        nnIdx = new ANNidx[nnode];
        dists = new ANNdist[nnode];
        T->annkFRSearch(f, r2, nnode, nnIdx, dists, eps); // search for n of them.
        
        for(int p1=0; p1<nnode; p1++)
        {
            int j = nnIdx[p1];
            if(L[j].atype.at(0) == 'H') continue;
            if(L[j].INTF == 0)
                continue;
            
            r = get_distance(R[i].axyz, L[j].axyz);
            if(r > 3.5 || r < 1.9) continue;
            
            if(L[j].atype.at(0) == 'S' && R[i].atype.at(0) == 'S')
            {
                if(r > 1.9 && r < 2.1)
                    e_ss += -1.;
            }
            
            if(R[i].don_acc == 0 || L[j].don_acc == 0) continue;
            
            // check if these are donor acceptor pairs
            if((R[i].don_acc == 2 && L[j].don_acc == 1) || (R[i].don_acc == 1 && L[j].don_acc == 2) ||
               (R[i].don_acc == 1 && L[j].don_acc == 3) || (R[i].don_acc == 3 && L[j].don_acc == 2))
            {
                if(r > 2.74 && r < 3.5)
                    hbp += (5.* pow(R0/r, 12.) - 6.*pow(R0/r, 10.));
            }
        }
        delete [] nnIdx;
        delete [] dists;
    }
    
    // add hbp + e_ss
    hbp += e_ss;
}

// distance cutoff <8A
// energy cutoff 1.0 kcal
// EMIN - kcal/mol
void compute_vdw_energy(vector<atom>& R, vector<atom>& L, ANNkd_tree* T, double& vdw, double& vdw_attr, double& vdw_rep)
{
    double sigma, r, z, eps, tmp;
    double t;
    double A = (1./pow(0.6,12.)) - 2. * (1./pow(0.6, 6.));
    double B = -12.*(1./pow(0.6, 13.)) + 12. * (1./pow(0.6, 7.));
    double f_i = 1.0, f_j = 1.0;
    
    ANNpoint f = annAllocPt(3);
    
    double r2 = 64.0;
    bool reject = false;
    
    for(size_t i=0;i<R.size();i++) // for the fixed protein
    {
        if(R[i].INTF == 0)      continue;
        
        f[0] = R[i].axyz[0];
        f[1] = R[i].axyz[1];
        f[2] = R[i].axyz[2];
        
        if(R[i].atype.at(0) == 'H')
            f_i = 0.4;
        else
            f_i = 1.0;
        
        ANNidxArray nnIdx = NULL;
        ANNdistArray dists = NULL;
        
        int nnode = T->annkFRSearch(f, r2, 0, nnIdx, dists, 0.); // search for n of them.
        if(nnode <= 0)      continue;
        
        nnIdx = new ANNidx[nnode];
        dists = new ANNdist[nnode];
        T->annkFRSearch(f, r2, nnode, nnIdx, dists, 0.); // search for n of them.
//        cout<<"no. of nodes = " << nnode<<endl;
        int count = 0, i_count = 0, l_count = 0, t_count = 0;
        for(int p1=0; p1<nnode; p1++)
        {
            int j = nnIdx[p1];
            
            if(L[j].INTF == 0)      {i_count++;continue;}
            
            r = get_distance(R[i].axyz, L[j].axyz);
            if(r > 8.0)     {l_count++;continue;}// from ZRANK
            if(L[j].atype.at(0) == 'H')
                f_j = 0.4;
            else
                f_j = .8;
            
            sigma = f_i * R[i].chpar.RMIN + f_j * L[j].chpar.RMIN;
            z = sigma/r;
            eps = sqrt(R[i].chpar.EMIN * L[j].chpar.EMIN);
            //printf("%s %s %s %s\n", R[i].atype.c_str(), R[i].residue.c_str(), L[j].atype.c_str(), L[j].residue.c_str());
            tmp = (eps*(pow(z, 12.) - 2.*pow(z, 6.)));
            if(tmp > 1.)
                {vdw += 1.; count++;}
            else
                t_count++;
            
            t = 0.6 * sigma;
            if(((L[j].atype.at(0) == 'H' && R[i].don_acc == 2) || (R[i].atype.at(0) == 'H' && L[j].don_acc == 2)) && r < 2.50)
                reject = true;
            else
                reject = false;
            
            
            if(r > t)
            {
                if(tmp < 0.)
                    vdw_attr += tmp;
                else
                {
                    if(!reject)
                        vdw_rep += tmp;
                }
            }
            else
            {
                if(!reject)
                {
                    tmp = eps*(A + ((r-t)*(B/sigma)));
                    vdw_rep += tmp;
                }
            }
        }
        if (count > 0) {
//            cout<<"no. of interface nodes = " << i_count <<endl;
//            cout<<"no. of nodes within dist. > 8 Angstrom = " << l_count <<endl;
//            cout<<"no. of nodes with tmp <= 1 = " << t_count <<endl;
//            cout<<"no. of nodes within 8 Angstrom = " << count <<endl;
        }

        delete [] nnIdx;
        delete [] dists;
    }
//    cout<<"vdw = " << (vdw) <<endl;

}

void set_polar_group_charges(string res, string type, double& q)
{
    if(res == "ASP")
    {
        if(type == "OD1" || type == "OD2")
            q = -0.5;
        return;
    }
    else if(res == "GLU")
    {
        if(type == "OE1" || type == "OE2")
            q = -0.5;
        return;
    }
    else if(res == "ARG")
    {
        if(type == "NH1" || type == "NH2")
            q = 0.5;
        return;
    }
    else if(res == "LYS")
    {
        if(type == "NZ")
            q = 1.;
        return;
    }
    else if(res == "CTER")
    {
        if(type == "OT1" || type == "OT2")
            q = -1.;
        return;
    }
    else if(res == "NTER")
    {
        if(type == "N")
            q = 1.;
        return;
    }
    else
    {
        q = 0.;
    }
}

// formula taken from pyDOCK Coulombic electrostatics with distance-dependent dielectric constant (eps = 4r),
// explicitly calculated for all intermolecular atom pairs, with q atomic charges from CHARMM19 electrostatics
// weighted sum of short and long range attractive and repulsive forces
// charges of only the most polar groups are considered (ARG, LYS, ASP, GLU, C-terminus, and N-terminus).
// a charge of -0.5 on each of the d-oxygen atoms in ASP and the e-oxygen atoms in Glu,
// a charge of +0.5 on the eta-nitrogen atoms in Arg, and a full positive charge +1.0 on the z-nitrogen atom in Lys.
// use CHARMM19 charges for distances < 5.0 (short range), long range use only fully charged side chains

// 332.17 = Unit conversion factor
void compute_electrostatics_energy(vector<atom>& R, vector<atom>& L, double& elec,
                                   double& elec_sr_attr, double& elec_sr_rep, double& elec_lr_attr, double& elec_lr_rep)
{
    double r, q, e, rmod, q1, q2;
    
    for(size_t i=0;i<R.size();i++)
    {
        if(R[i].INTF == 0)
            continue;
        for(size_t j=0;j<L.size();j++)
        {
            if(L[j].INTF == 0)
                continue;
            
            q = R[i].atomcharge * L[j].atomcharge;
            if(q == 0.)
                continue;
            r = get_squared_distance(R[i].axyz, L[j].axyz);
            
            rmod = MAX(r, 9.);//to avoid singularities as r_ij->0
            e = q/rmod;
            // truncation
            if(e > 1.)
                elec += 1.; // set interatomic potential to 1. kcal
            else if(e < -1.)
                elec += -1.; // set interatomic potential to -1. kcal
            else
                elec += e;
            
            
            if(r < 25.)
            {
                if(e < 0.)
                    elec_sr_attr += e;
                else
                    elec_sr_rep += e;
            }
            else if(r >= 25.) // only charges on the side chain atoms
            {
                set_polar_group_charges(R[i].residue, R[i].atype, q1);
                set_polar_group_charges(L[j].residue, L[j].atype, q2);
                if(q1 == 0. || q2 == 0.)
                    continue;
                q = q1*q2;
                
                e = q/(r);
                if(e < 0.)
                    elec_lr_attr += e;
                else
                    elec_lr_rep += e;
            }
        }
    }
    
    elec *= EFACTOR;
    elec_lr_attr *= EFACTOR;
    elec_lr_rep *= EFACTOR;
    elec_sr_attr *= EFACTOR;
    elec_sr_rep *= EFACTOR;
}


// Solvation computed using Gaussian solvent exclusion
// Themis Lazaridis and Martin Karplus
// Solvation computed using atom wise solvent accesssible area
// solvation parameters in cal/mol/A^2
// divide by 1000
// based on atomic desolvation parameters of Abagyan
void compute_solvation(vector<atom>& AB, double& esasa)
{
    get_atomic_sasa(AB);
    esasa = 0.;
    for(size_t i=0;i<AB.size();i++)
    {
        if(AB[i].atype.at(0) == 'H') continue;
        esasa += (AB[i].sas_area * AB[i].atom_solv);
    }
    esasa /= 1000.0; // conversion to kcal/mol
}


void compute_LZ_solvation(vector<atom>& R, vector<atom>& L, double& solv)
{
    double r, f, n1, n2, t;
    int count1=0,count2=0;
    for(size_t i=0;i<R.size();i++)
    {
        //printf("%s %s %d %d\n", R[i].atype.c_str(), R[i].residue.c_str(), R[i].INTF, R[i].chpar.valset);
        if(R[i].atype.at(0) == 'H') continue;
        if(R[i].INTF == 0)          continue;
        if(!R[i].chpar.valset)      continue;
        
        // This line was added to correct a calculation error, since the DGREF was not being considered
        double tmpsolv = 0;
        count1++;
        count2=0;
        for(size_t j=0;j<L.size();j++)
        {
            //printf("%s %s %s %s\n", R[i].atype.c_str(), R[i].residue.c_str(), L[j].atype.c_str(), L[j].residue.c_str());
            if(L[j].atype.at(0) == 'H') continue;
            if(L[j].INTF == 0)          continue;
            if(!L[j].chpar.valset)      continue;
            
            r = get_distance(R[i].axyz, L[j].axyz);
            if(r > 8.0)     continue;// from ZRANK
            count2++;
            n1 = 0.; n2 = 0;
            
            n1 = (R[i].chpar.DGFREE)/(r*r*R[i].chpar.LAMBDA);
            t = (r - R[i].chpar.RMIN)/R[i].chpar.LAMBDA;
            f = exp(-(t*t));
            f *= L[j].chpar.AVOL;
            n1 *= f;
            
            n2 = (L[j].chpar.DGFREE)/(r*r*L[j].chpar.LAMBDA);
            t = (r - L[j].chpar.RMIN)/L[j].chpar.LAMBDA;
            f = exp(-t*t);
            f *= (R[i].chpar.AVOL);
            n2 *= f;
            // The following line represented the old way of adding, before the calculation correction was done
            // solv += (n1 + n2)/Z;
            // This line was added to correct a calculation error, since the DGREF was not being considered
            tmpsolv += (n1 + n2)/Z;
        }
        solv += (R[i].chpar.DGREF - tmpsolv);
    }
}


/*
 * This function goes through all combinations of receptor/ligand in a multi-chain complex and calculates the
 * soroban score for it.
 */
soroban_score compute_energy(vector< vector<atom> >& proteins, double* weights)
{
    soroban_score score(weights);
    
    // we need to join all vectors into one in order to calculate the solvent accesibility term
    // the outer loop will add the atoms to "all"
   // vector<atom> all;
    // These two vectors hold the current receptor/ligand on each iteration
    vector<atom> R, L;
    //cout<<"compute_energy 1 \n";
    for(size_t receptor_index = 0; receptor_index < (proteins.size() - 1); receptor_index++)
    {
        R = proteins[receptor_index];
        L.clear();
        // create a complex with all the other proteins
        for(size_t ligand_index = receptor_index + 1; ligand_index < proteins.size(); ligand_index++)
        {
            if(ligand_index != receptor_index)
            {
                L.insert(L.end(), proteins[ligand_index].begin(), proteins[ligand_index].end());
            }
        }
        // mark the interface atoms before invoking the scoring functions, since they use this information
        get_interface_residues(R, L);
        int N = (int) L.size(), dim = 3;
        ANNpointArray dataPts = annAllocPts(N, dim);
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<dim;j++)
                dataPts[i][j] = L[i].axyz[j];
        }
        ANNkd_tree *T = new ANNkd_tree(dataPts, N, dim);
        
        // compute overall vdw
        double vdw = 0., vdw_attr = 0., vdw_rep = 0.;
        compute_vdw_energy(R, L, T, vdw, vdw_attr, vdw_rep);
        score.add_vdw(vdw, vdw_attr, vdw_rep);

        //compute overall electrostatics
        double elec = 0., elec_sr_attr = 0., elec_sr_rep = 0., elec_lr_attr = 0., elec_lr_rep = 0.;
        compute_electrostatics_energy(R, L, elec, elec_sr_attr, elec_sr_rep, elec_lr_attr, elec_lr_rep);
        score.add_electrostatics(elec, elec_sr_attr, elec_sr_rep, elec_lr_attr, elec_lr_rep);

        //hydrogen bonding potential + disulphide
        double hbp_ss = 0.;
        //compute_hbp(R, L, T, hbp_ss);
        score.add_hbp(hbp_ss);

        //lz solvation
        double solv = 0.;
        compute_LZ_solvation(R, L, solv);
        score.add_lz_solvation(solv);

        //desolvation based on atom contact potential
        double acp = 0.;
        compute_acp(R, L, T, acp);
        score.add_acp_based_desolvation(acp);

        delete T;
        annDeallocPts(dataPts);
        annClose();
        // add the atoms to all for later solvation calculation
       // all.insert(all.end(), R.begin(), R.end());
    }
    //cout<<"compute_energy 2 \n";

    // add the last chain manually because in the loop we go from 0 to (N-1)
    // and we also need to add the n-th chain for the sasa computation
    R = proteins[proteins.size() - 1];
    //all.insert(all.end(), R.begin(), R.end());
    
    // compute solvation energy based on accessible surface area
    double sasa = 0.;
    //compute_solvation(all, sasa);
    score.add_solvation(sasa);
    //cout<<"compute_energy 3 \n";

    //swap
    //vector<atom>().swap( all );
    vector<atom>().swap( R );
    vector<atom>().swap( L );
    return score;
}
