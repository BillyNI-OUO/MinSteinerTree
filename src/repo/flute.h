#define POWVFILE "POWV9.dat"   
#define POSTFILE "POST9.dat"   
#define MAXD 10000              
#define D 9                    
#define ROUTING 1              
#define accuracyURACY 3             
#define REMOVE_DUPLICATE_PIN 0
#define LOCAL_REFINEMENT 0     


#define DATATYPE int

typedef struct branch{
    DATATYPE x, y;
    int n;      
} Branch;

typedef struct tree{
    Branch *branch;        
    DATATYPE length;   
    int deg;
    int getDegree(){
        return deg;
    };
} Tree;

typedef struct point{
    DATATYPE x, y;
    int o;
} Point;



int comparePointX(const void *left, const void *right){
    point *pLeft, *pRight;

    pLeft = *(Point**)left;
    pRight = *(Point**)right;

    if(pLeft->x > pRight->x){
        return 1;
    }else if(pLeft->x < pRight->x){
        return -1;
    }else{
        return 0;
    }
}

int comparePointY(const void *left, const void *right){
    point *pLeft, *pRight;

    pLeft = *(Point**)left;
    pRight = *(Point**)right;

    if(pLeft->y > pRight->y){
        return 1;
    }else if(pLeft->y < pRight->y){
        return -1;
    }else{
        return 0;
    }
}

extern void readLUT();
extern DATATYPE flute_wl(int d, DATATYPE x[], DATATYPE y[], int accuracy);
extern Tree flute(int d, DATATYPE x[], DATATYPE y[], int accuracy);
extern DATATYPE wirelength(Tree t);
extern void printTree(Tree t);
extern void plotTree(Tree t);


extern DATATYPE flutes_wl_LD(int d, DATATYPE xCors[], DATATYPE yCors[], int s[]);
extern DATATYPE flutes_wl_MD(int d, DATATYPE xCors[], DATATYPE yCors[], int s[], int accuracy);
extern DATATYPE flutes_wl_RDP(int d, DATATYPE xCors[], DATATYPE yCors[], int s[], int accuracy);
extern Tree flutes_LD(int d, DATATYPE xCors[], DATATYPE yCors[], int s[]);
extern Tree flutes_MD(int d, DATATYPE xCors[], DATATYPE yCors[], int s[], int accuracy);
extern Tree flutes_RDP(int d, DATATYPE xCors[], DATATYPE yCors[], int s[], int accuracy);

#if REMOVE_DUPLICATE_PIN == 1
#define flutes_wl(d, xCors, yCors, s, accuracy) flutes_wl_RDP(d, xCors, yCors, s, accuracy)
#define flutes(d, xCors, yCors, s, accuracy) flutes_RDP(d, xCors, yCors, s, accuracy)
#else
#define flutes_wl(d, xCors, yCors, s, accuracy) flutes_wl_ALLD(d, xCors, yCors, s, accuracy)
#define flutes(d, xCors, yCors, s, accuracy) flutes_ALLD(d, xCors, yCors, s, accuracy)
#endif

#define flutes_wl_ALLD(d, xCors, yCors, s, accuracy) flutes_wl_LMD(d, xCors, yCors, s, accuracy)
#define flutes_ALLD(d, xCors, yCors, s, accuracy) flutes_LMD(d, xCors, yCors, s, accuracy)

#define flutes_wl_LMD(d, xCors, yCors, s, accuracy) \
    (d <= D ? flutes_wl_LD(d, xCors, yCors, s) : flutes_wl_MD(d, xCors, yCors, s, accuracy))
#define flutes_LMD(d, xCors, yCors, s, accuracy) \
    (d <= D ? flutes_LD(d, xCors, yCors, s) : flutes_MD(d, xCors, yCors, s, accuracy))

#define max(x, y) (x > y ? x : y)
#define min(x, y) (x < y ? x : y)
#define ADIFF(x, y) (x > y ? x - y : y - x)
