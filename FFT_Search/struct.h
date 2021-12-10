#define VER 0.43
#define PI 3.141592
#define ATOM 50000
#define RES 5000
#define AMI 20
#define LIN 256
#define radian 57.29577951308232088

//grid
#define GRINO 20000

#define OUTSIDE -999.999

/*def �δĶ�����*/
#define INC 15/*ligand��ž���ٿ�*/
#define MOV 1/*ligand��ư���ٿ�*/
#define START_RT 0 /*ligand��ž��start*/
#define START_MV 0 /*ligand��ư��start*/
#define FIN_RT 360
#define FIN_MV 0
#define ON 1
#define OFF 0
#define TRUE 0
#define FALSE -1
/*#define MATFILE "result030106"*/ /*�������ޥȥꥯ���Υե�����*/
#define MATSIZE 1500
//#define MATLEN 600
#define MATAMI 24 /*�ޥȥꥯ���Υ��ߥλ�����*/

#define CSHNUM 20 /*���ͻĴ��������*/

/*�ޥ�������*/
#define X(a) (a)*(a)
#define L(a,b,c,d,e,f) (sqrt((d-a)*(d-a)+(e-b)*(e-b)+(f-c)*(f-c)))/*�٥��ȥ�Ĺ*/
#define RAS(a) (2.000000*PI*a/360.000000)/*��->�饸����*/
#define RAD(a) (2.000000*PI*a/360.000000)/*��->�饸����*/
#define ANG(a) (360.000000*a/(2.000000*PI))/*�饸����->��*/

#define GAUSS(a,b) 1.00/(sqrt(2.00*PI)*b)*exp(-(a*a)/(2.0*b*b))

/*triangle*/

#define MAXMTX   4

/*��̽���*/
#define TOP 3000

#define NOT_GAUSS 1
#define USE_GAUSS 0

#define VALIABLE 0
#define CONSTANT 1

int FINAL_ATNo;
int FINAL_RENo;
typedef struct{
	        //double x,y,z;
	        float x,y,z;
}COORD;


typedef struct{
        char fname[LIN];
        int NumOfAtom, NumOfRes;
        int TypeAtomOder[RES][17];
        float *Charge; 
        COORD *coord;
        COORD CAcd[RES];
        COORD *CBcd;
        COORD *Cen;
        COORD *Intra;//interaction
        int NumOfIntra;
        COORD *Nonin;//non intra
        int NumOfNonin;
        float *phi,*psi;

        //HETATM
        char **HET_TypeAtom;
        COORD *HET_coord;
        int NumOfHet;
        int NumOfReal;
        int RealResNum[RES];
        //New 2012.10.29-----------------
        double **xyz;
        int *TypeAtomId,*TypeResId;
        char **TypeAtom, **TypeRes, *Chain,**RealNum;
        int *AtomOnRes, *SosuiAtom, *ConservedAtom;
        int *ResOnAtom;
        int *ResNum,*AtomNum;
        double MaxXyz[3],MinXyz[3];
	float *DepthAtom,*DepthRes;
} PDB;


typedef struct{
	char filename[LIN];
	char file1[LIN],file2[LIN];
	double map_t;
	int Nthr;
	double dreso,LocalR;
	double MaxShift,MergeDist;
	double Filter,Dkeep;
	int Nround,Ntabu,Nnb,Nsim,Nbeam;
	int Mode,Amode;
	double Allow;
	double th1,th2;
	double ssize,ang;
	int TopN;
	bool ShowGrid;
    int mapType; //1: experimental, 2: simulated
} CMD;
