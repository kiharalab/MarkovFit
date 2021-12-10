typedef struct{
 int t[3];
 double r[3];
 double q[4];
 int code;
 double sco;
} TBL;

bool SearchMAPfft(MRC *,MRC *,double);
bool SearchMAPfftMT(MRC *,MRC *,int);
bool SearchMAPfftMT_OVCC(MRC *,MRC *,int,int); //Overlap or CCC
double GetScore(MRC *, MRC *, int [3],double *, int);
