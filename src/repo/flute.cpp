#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include "flute.h"



#if D <= 7
#define MGROUP 5040 / 4 
#define MPOWV 15        
#elif D == 8
#define MGROUP 40320 / 4 
#define MPOWV 33         
#elif D == 9
#define MGROUP 362880 / 4 
#define MPOWV 79          
#endif

struct Csoln
{
    unsigned char rowcol[D - 2];
    unsigned char parent; 
    unsigned char neighbor[2 * D - 2];
    unsigned char seg[11];
};
struct Csoln *LUT[D + 1][MGROUP]; // storing 4 .. D
int numSoln[D + 1][MGROUP];
unsigned char charNumMap[256];
int numGRPMap[10] = {0, 0, 0, 0, 6, 30, 180, 1260, 10080, 90720};
typedef struct pairArray{
    std::vector<DATATYPE>x;
    std::vector<DATATYPE>y;
}pairArray;

void readLUT();

Tree flute(int d, DATATYPE x[], DATATYPE y[], int accuracy);
Tree flutes_LD(int d, DATATYPE xCors[], DATATYPE yCors[], int s[]);
Tree flutes_MD(int d, DATATYPE xCors[], DATATYPE yCors[], int s[], int accuracy);
Tree flutes_RDP(int d, DATATYPE xCors[], DATATYPE yCors[], int s[], int accuracy);
DATATYPE flute_wl(int d, DATATYPE x[], DATATYPE y[], int accuracy);
DATATYPE flutes_wl_LD(int d, DATATYPE xCors[], DATATYPE yCors[], int s[]);
DATATYPE flutes_wl_MD(int d, DATATYPE xCors[], DATATYPE yCors[], int s[], int accuracy);
DATATYPE flutes_wl_RDP(int d, DATATYPE xCors[], DATATYPE yCors[], int s[], int accuracy);
Tree dmergetree(Tree t1, Tree t2);
Tree hmergetree(Tree t1, Tree t2, int s[]);
Tree vmergetree(Tree t1, Tree t2);
void local_refinement(Tree *tp, int p);
DATATYPE wirelength(Tree t);
void printtree(Tree t);
void plottree(Tree t);
void initCharNumMap();
pairArray treeToPairArray(Tree t);


void initCharNumMap(){
    for (int i = 0; i <= 255; i++){
        if (i >= '0'&& i <= '9'){
            charNumMap[i] = i - '0';
        }
        else if (i >= 'A'){
            charNumMap[i] = i - 'A' + 10;
        }
        else{
            charNumMap[i] = 0;
        }
    }
}

void readLUT(){
    unsigned char line[32], *linep, c;
    FILE *fpwv = NULL, *fprt = NULL;
    struct Csoln *p;
    int d, ns, nn;
    int j;

    fpwv = fopen(POWVFILE, "r");
    if (fpwv == NULL){
        std::cerr << "Error in opening " << POWVFILE << std::endl;
        exit(1);
    }
    
#if ROUTING == 1
    fprt = fopen(POSTFILE, "r");
    if (fprt == NULL){
        std::cerr << "Error in opening " << POSTFILE << std::endl;
        exit(1);
    }
#endif
    initCharNumMap();
    int k, k_index;
    for (d = 4; d <= D; d++){
        fscanf(fpwv, "d=%d\n", &d);
#if ROUTING == 1
        fscanf(fprt, "d=%d\n", &d);
#endif
        
        for (k = 0; k < numGRPMap[d]; k++)
        {
            ns = (int)charNumMap[fgetc(fpwv)];

            if (ns == 0){
                fscanf(fpwv, "%d\n", &k_index);
                LUT[d][k] = LUT[d][k_index];
                numSoln[d][k] = numSoln[d][k_index];
            }
            else{
                fgetc(fpwv);
                numSoln[d][k] = ns;
                p = (struct Csoln *)malloc(ns * sizeof(struct Csoln));
                LUT[d][k] = p;
                for (int i = 1; i <= ns; i++){
                    linep = (unsigned char *)fgets((char *)line, 32, fpwv);
                    p->parent = charNumMap[*(linep++)];
                    j = 0;
                    while ((p->seg[j++] = charNumMap[*(linep++)]) != 0);
                    j = 10;
                    while ((p->seg[j--] = charNumMap[*(linep++)]) != 0);
#if ROUTING == 1
                    nn = 2 * d - 2;
                    fread(line, 1, d - 2, fprt);
                    linep = line;
                    for (j = d; j < nn; j++)
                    {
                        c = charNumMap[*(linep++)];
                        p->rowcol[j - d] = c;
                    }
                    fread(line, 1, nn / 2 + 1, fprt);
                    linep = line;
                    for (int j = 0; j < nn;)
                    {
                        c = *(linep++);
                        
                        p->neighbor[j++] = c / 16;
                        p->neighbor[j++] = c % 16;
                        
                    }
#endif
                    p++;
                }
            }
        }
    }
    if (fpwv)
        fclose(fpwv);
    if (fprt)
        fclose(fprt);
}

DATATYPE flute_wl(int d, DATATYPE x[], DATATYPE y[], int accuracy)
{
    DATATYPE xCors[MAXD], yCors[MAXD];
    DATATYPE l, xu, xl, yu, yl;
    int s[MAXD];
    int i, j, k;
    int minIdex = 0;
    DATATYPE minVal = INT_MAX;
    Point pt[MAXD], *pointDP[MAXD], *tmpp;

    if (d == 2)
        l = ADIFF(x[0], x[1]) + ADIFF(y[0], y[1]);
    else if (d == 3)
    {
        if (x[0] > x[1]){
            xu = max(x[0], x[2]);
            xl = min(x[1], x[2]);
        }
        else{
            xu = max(x[1], x[2]);
            xl = min(x[0], x[2]);
        }
        if (y[0] > y[1]){
            yu = max(y[0], y[2]);
            yl = min(y[1], y[2]);
        }
        else{
            yu = max(y[1], y[2]);
            yl = min(y[0], y[2]);
        }
        l = (xu - xl) + (yu - yl);
    }
    else{
        for (int i = 0; i < d; i++){
            pt[i].x = x[i];
            pt[i].y = y[i];
            pointDP[i] = &pt[i];
        }

        
        for (i = 0; i < d - 1; i++){
            minVal = pointDP[i]->x;
            minIdex = i;
            for (j = i + 1; j < d; j++)
            {
                if (minVal > pointDP[j]->x)
                {
                    minVal = pointDP[j]->x;
                    minIdex = j;
                }
            }
            tmpp = pointDP[i];
            pointDP[i] = pointDP[minIdex];
            pointDP[minIdex] = tmpp;
        }
        //qsort(pointDP, d, sizeof(Point *), comparePointX);

#if REMOVE_DUPLICATE_PIN == 1
        pointDP[d] = &pt[d];
        pointDP[d]->x = pointDP[d]->y = INT_MIN;
        j = 0;
        for (i = 0; i < d; i++){
            for (k = i + 1; pointDP[k]->x == pointDP[i]->x; k++)
                if (pointDP[k]->y == pointDP[i]->y)
                    break;
            if (pointDP[k]->x != pointDP[i]->x)
                pointDP[j++] = pointDP[i];
        }
        d = j;
#endif

        for (i = 0; i < d; i++){
            xCors[i] = pointDP[i]->x;
            pointDP[i]->o = i;
        }

        for (i = 0; i < d - 1; i++){
            minVal = pointDP[i]->y;
            minIdex = i;
            for (j = i + 1; j < d; j++){
                if (minVal > pointDP[j]->y){
                    minVal = pointDP[j]->y;
                    minIdex = j;
                }
            }
            yCors[i] = pointDP[minIdex]->y;
            s[i] = pointDP[minIdex]->o;
            pointDP[minIdex] = pointDP[i];
        }

        //std::qsort(pointDP, d, sizeof(Point *), comparePointY);
        yCors[d - 1] = pointDP[d - 1]->y;
        s[d - 1] = pointDP[d - 1]->o;

        l = flutes_wl(d, xCors, yCors, s, accuracy);
    }
    return l;
}


DATATYPE flutes_wl_RDP(int d, DATATYPE xCors[], DATATYPE yCors[], int s[], int accuracy)
{
    int i, j, ss;

    for (i = 0; i < d - 1; i++)
    {
        if (xCors[s[i]] == xCors[s[i + 1]] && yCors[i] == yCors[i + 1])
        {
            if (s[i] < s[i + 1])
                ss = s[i + 1];
            else
            {
                ss = s[i];
                s[i] = s[i + 1];
            }
            for (j = i + 2; j < d; j++)
            {
                yCors[j - 1] = yCors[j];
                s[j - 1] = s[j];
            }
            for (j = ss + 1; j < d; j++)
                xCors[j - 1] = xCors[j];
            for (j = 0; j <= d - 2; j++)
                if (s[j] > ss)
                    s[j]--;
            i--;
            d--;
        }
    }
    return flutes_wl_ALLD(d, xCors, yCors, s, accuracy);
}

DATATYPE flutes_wl_LD(int d, DATATYPE xCors[], DATATYPE yCors[], int s[])
{
    int k, pi, i, j;
    struct Csoln *rlist;
    DATATYPE dd[2 * D - 2]; // 0..D-2 for v, D-1..2*D-3 for h
    DATATYPE minl, sum, l[MPOWV + 1];

    if (d <= 3)
        minl = xCors[d - 1] - xCors[0] + yCors[d - 1] - yCors[0];
    else
    {
        k = 0;
        if (s[0] < s[2])
            k++;
        if (s[1] < s[2])
            k++;

        for (i = 3; i <= d - 1; i++)
        { // p0=0 alwayCors, skip i=1 for symmetry
            pi = s[i];
            for (j = d - 1; j > i; j--)
                if (s[j] < s[i])
                    pi--;
            k = pi + (i + 1) * k;
        }

        if (k < numGRPMap[d]) // no horizontal flip
            for (i = 1; i <= d - 3; i++)
            {
                dd[i] = yCors[i + 1] - yCors[i];
                dd[d - 1 + i] = xCors[i + 1] - xCors[i];
            }
        else
        {
            k = 2 * numGRPMap[d] - 1 - k;
            for (i = 1; i <= d - 3; i++)
            {
                dd[i] = yCors[i + 1] - yCors[i];
                dd[d - 1 + i] = xCors[d - 1 - i] - xCors[d - 2 - i];
            }
        }

        minl = l[0] = xCors[d - 1] - xCors[0] + yCors[d - 1] - yCors[0];
        rlist = LUT[d][k];
        for (i = 0; rlist->seg[i] > 0; i++)
            minl += dd[rlist->seg[i]];

        l[1] = minl;
        j = 2;
        while (j <= numSoln[d][k])
        {
            rlist++;
            sum = l[rlist->parent];
            for (i = 0; rlist->seg[i] > 0; i++)
                sum += dd[rlist->seg[i]];
            for (i = 10; rlist->seg[i] > 0; i--)
                sum -= dd[rlist->seg[i]];
            minl = min(minl, sum);
            l[j++] = sum;
        }
    }

    return minl;
}

// For medium-degree, i.e., D+1 <= d
DATATYPE flutes_wl_MD(int d, DATATYPE xCors[], DATATYPE yCors[], int s[], int accuracy)
{
    DATATYPE x1[MAXD], x2[MAXD], y1[MAXD], y2[MAXD];
    int si[MAXD], s1[MAXD], s2[MAXD];
    float score[2 * MAXD], penalty[MAXD], pnlty, dx, dy;
    DATATYPE ll, minl, extral;
    int i, r, p, maxbp, nbp, bp, ub, lb, n1, n2, newaccuracy;
    int ms, mins, maxCors, minsi, maxCorsi;
    DATATYPE distx[MAXD], disty[MAXD], xydiff;

    if (s[0] < s[d - 1])
    {
        ms = max(s[0], s[1]);
        for (i = 2; i <= ms; i++)
            ms = max(ms, s[i]);
        if (ms <= d - 3)
        {
            for (i = 0; i <= ms; i++)
            {
                x1[i] = xCors[i];
                y1[i] = yCors[i];
                s1[i] = s[i];
            }
            x1[ms + 1] = xCors[ms];
            y1[ms + 1] = yCors[ms];
            s1[ms + 1] = ms + 1;

            s2[0] = 0;
            for (i = 1; i <= d - 1 - ms; i++)
                s2[i] = s[i + ms] - ms;

            return flutes_wl_LMD(ms + 2, x1, y1, s1, accuracy) + flutes_wl_LMD(d - ms, xCors + ms, yCors + ms, s2, accuracy);
        }
    }
    else{
        ms = min(s[0], s[1]);
        for (i = 2; i <= d - 1 - ms; i++)
            ms = min(ms, s[i]);
        if (ms >= 2)
        {
            x1[0] = xCors[ms];
            y1[0] = yCors[0];
            s1[0] = s[0] - ms + 1;
            for (i = 1; i <= d - 1 - ms; i++)
            {
                x1[i] = xCors[i + ms - 1];
                y1[i] = yCors[i];
                s1[i] = s[i] - ms + 1;
            }
            x1[d - ms] = xCors[d - 1];
            y1[d - ms] = yCors[d - 1 - ms];
            s1[d - ms] = 0;

            s2[0] = ms;
            for (i = 1; i <= ms; i++)
                s2[i] = s[i + d - 1 - ms];

            return flutes_wl_LMD(d + 1 - ms, x1, y1, s1, accuracy) + flutes_wl_LMD(ms + 1, xCors, yCors + d - 1 - ms, s2, accuracy);
        }
    }

    for (r = 0; r < d; r++){
        si[s[r]] = r;
    }

    lb = (d - 2 * accuracy + 2) / 4;
    if (lb < 2){
        lb = 2;
    }
    ub = d - 1 - lb;

#define AAWL 0.6
#define BBWL 0.3
    float CCWL = 7.4 / ((d + 10.) * (d - 3.));
    float DDWL = 4.8 / (d - 1);

    dx = CCWL * (xCors[d - 2] - xCors[1]);
    dy = CCWL * (yCors[d - 2] - yCors[1]);
    for (r = d / 2, pnlty = 0; r >= 0; r--, pnlty += dx)
        penalty[r] = pnlty, penalty[d - 1 - r] = pnlty;
    for (r = d / 2 - 1, pnlty = dy; r >= 0; r--, pnlty += dy)
        penalty[s[r]] += pnlty, penalty[s[d - 1 - r]] += pnlty;
    xydiff = (xCors[d - 1] - xCors[0]) - (yCors[d - 1] - yCors[0]);
    if (s[0] < s[1])
        mins = s[0], maxCors = s[1];
    else
        mins = s[1], maxCors = s[0];
    if (si[0] < si[1])
        minsi = si[0], maxCorsi = si[1];
    else
        minsi = si[1], maxCorsi = si[0];
    for (r = 2; r <= ub; r++)
    {
        if (s[r] < mins)
            mins = s[r];
        else if (s[r] > maxCors)
            maxCors = s[r];
        distx[r] = xCors[maxCors] - xCors[mins];
        if (si[r] < minsi)
            minsi = si[r];
        else if (si[r] > maxCorsi)
            maxCorsi = si[r];
        disty[r] = yCors[maxCorsi] - yCors[minsi] + xydiff;
    }

    if (s[d - 2] < s[d - 1])
        mins = s[d - 2], maxCors = s[d - 1];
    else
        mins = s[d - 1], maxCors = s[d - 2];
    if (si[d - 2] < si[d - 1])
        minsi = si[d - 2], maxCorsi = si[d - 1];
    else
        minsi = si[d - 1], maxCorsi = si[d - 2];
    for (r = d - 3; r >= lb; r--)
    {
        if (s[r] < mins)
            mins = s[r];
        else if (s[r] > maxCors)
            maxCors = s[r];
        distx[r] += xCors[maxCors] - xCors[mins];
        if (si[r] < minsi)
            minsi = si[r];
        else if (si[r] > maxCorsi)
            maxCorsi = si[r];
        disty[r] += yCors[maxCorsi] - yCors[minsi];
    }

    nbp = 0;
    for (r = lb; r <= ub; r++)
    {
        if (si[r] == 0 || si[r] == d - 1)
            score[nbp] = (xCors[r + 1] - xCors[r - 1]) - penalty[r] - AAWL * (yCors[d - 2] - yCors[1]) - DDWL * disty[r];
        else
            score[nbp] = (xCors[r + 1] - xCors[r - 1]) - penalty[r] - BBWL * (yCors[si[r] + 1] - yCors[si[r] - 1]) - DDWL * disty[r];
        nbp++;

        if (s[r] == 0 || s[r] == d - 1)
            score[nbp] = (yCors[r + 1] - yCors[r - 1]) - penalty[s[r]] - AAWL * (xCors[d - 2] - xCors[1]) - DDWL * distx[r];
        else
            score[nbp] = (yCors[r + 1] - yCors[r - 1]) - penalty[s[r]] - BBWL * (xCors[s[r] + 1] - xCors[s[r] - 1]) - DDWL * distx[r];
        nbp++;
    }

    if (accuracy <= 3)
        newaccuracy = 1;
    else
    {
        newaccuracy = accuracy / 2;
        if (accuracy >= nbp)
            accuracy = nbp - 1;
    }

    minl = (DATATYPE)INT_MAX;
    for (i = 0; i < accuracy; i++)
    {
        maxbp = 0;
        for (bp = 1; bp < nbp; bp++)
            if (score[maxbp] < score[bp])
                maxbp = bp;
        score[maxbp] = -9e9;

#define BreakPt(bp) ((bp) / 2 + lb)
#define BreakInX(bp) ((bp) % 2 == 0)
        p = BreakPt(maxbp);
        if (BreakInX(maxbp)){
            n1 = n2 = 0;
            for (r = 0; r < d; r++){
                if (s[r] < p){
                    s1[n1] = s[r];
                    y1[n1] = yCors[r];
                    n1++;
                }
                else if (s[r] > p){
                    s2[n2] = s[r] - p;
                    y2[n2] = yCors[r];
                    n2++;
                }
                else{
                    s1[n1] = p;
                    s2[n2] = 0;
                    if (r == d - 1 || r == d - 2){
                        y1[n1] = y2[n2] = yCors[r - 1];
                        extral = yCors[r] - yCors[r - 1];
                    }
                    if (r == 0 || r == 1){
                        y1[n1] = y2[n2] = yCors[r + 1];
                        extral = yCors[r + 1] - yCors[r];
                    }
                    else{
                        y1[n1] = y2[n2] = yCors[r];
                        extral = 0;
                    }
                    n1++;
                    n2++;
                }
            }
            ll = extral + flutes_wl_LMD(p + 1, xCors, y1, s1, newaccuracy) + flutes_wl_LMD(d - p, xCors + p, y2, s2, newaccuracy);
        }
        else{
            n1 = n2 = 0;
            for (r = 0; r < d; r++){
                if (si[r] < p){
                    s1[si[r]] = n1;
                    x1[n1] = xCors[r];
                    n1++;
                }
                else if (si[r] > p){
                    s2[si[r] - p] = n2;
                    x2[n2] = xCors[r];
                    n2++;
                }
                else{
                    s1[p] = n1;
                    s2[0] = n2;
                    if (r == d - 1 || r == d - 2){
                        x1[n1] = x2[n2] = xCors[r - 1];
                        extral = xCors[r] - xCors[r - 1];
                    }
                    if (r == 0 || r == 1)
                    {
                        x1[n1] = x2[n2] = xCors[r + 1];
                        extral = xCors[r + 1] - xCors[r];
                    }
                    else
                    {
                        x1[n1] = x2[n2] = xCors[r];
                        extral = 0;
                    }
                    n1++;
                    n2++;
                }
            }
            ll = extral + flutes_wl_LMD(p + 1, x1, yCors, s1, newaccuracy) + flutes_wl_LMD(d - p, x2, yCors + p, s2, newaccuracy);
        }
        if (minl > ll)
            minl = ll;
    }
    return minl;
}

Tree flute(int d, DATATYPE x[], DATATYPE y[], int accuracy){
    DATATYPE *xCors, *yCors, minVal;
    int *s;
    int i, j, k, minIdex;
    Point *pt, **pointDP, *tmpp;
    Tree t;
    
    if (d == 2)
    {
        t.deg = 2;
        t.length = ADIFF(x[0], x[1]) + ADIFF(y[0], y[1]);
        t.branch = (Branch *)malloc(2 * sizeof(Branch));
        t.branch[0].x = x[0];
        t.branch[0].y = y[0];
        t.branch[0].n = 1;
        t.branch[1].x = x[1];
        t.branch[1].y = y[1];
        t.branch[1].n = 1;
    }
    else
    {
        xCors = (DATATYPE *)malloc(sizeof(DATATYPE)*(d));
        yCors = (DATATYPE *)malloc(sizeof(DATATYPE)*(d));
        s = (int *)malloc(sizeof(int)*(d));
        pt = (Point *)malloc(sizeof(Point)*(d+1));
        pointDP = (Point **)malloc(sizeof(Point*)*(d+1));
        for (int i = 0; i < d; i++){
            pt[i].x = x[i];
            pt[i].y = y[i];
            pointDP[i] = &pt[i];
        }
        for (i = 0; i < d - 1; i++)
        {
            minVal = pointDP[i]->x;
            minIdex = i;
            for (j = i + 1; j < d; j++)
            {
                if (minVal > pointDP[j]->x)
                {
                    minVal = pointDP[j]->x;
                    minIdex = j;
                }
            }
            tmpp = pointDP[i];
            pointDP[i] = pointDP[minIdex];
            pointDP[minIdex] = tmpp;
        }

        //std::qsort(pointDP, d, sizeof(Point), comparePointX);

#if REMOVE_DUPLICATE_PIN == 1
        pointDP[d] = &pt[d];
        pointDP[d]->x = pointDP[d]->y = INT_MIN;
        j = 0;
        for (int i = 0; i < d; i++)
        {
            for (k = i + 1; pointDP[k]->x == pointDP[i]->x; k++)
                if (pointDP[k]->y == pointDP[i]->y) // pins k and i are the same
                    break;
            if (pointDP[k]->x != pointDP[i]->x)
                pointDP[j++] = pointDP[i];
        }
        d = j;
#endif

        for (i = 0; i < d; i++)
        {
            xCors[i] = pointDP[i]->x;
            pointDP[i]->o = i;
        }

        for (i = 0; i < d - 1; i++)
        {
            minVal = pointDP[i]->y;
            minIdex = i;
            for (j = i + 1; j < d; j++)
            {
                if (minVal > pointDP[j]->y)
                {
                    minVal = pointDP[j]->y;
                    minIdex = j;
                }
            }
            yCors[i] = pointDP[minIdex]->y;
            s[i] = pointDP[minIdex]->o;
            pointDP[minIdex] = pointDP[i];
        }

        //std::qsort(pointDP, d, sizeof(Point), comparePointY);

        yCors[d - 1] = pointDP[d - 1]->y;
        s[d - 1] = pointDP[d - 1]->o;
        
        t = flutes(d, xCors, yCors, s, accuracy);
        
        free(xCors);
        free(yCors);
        free(s);
        free(pt);
        free(pointDP);
        
    }
    return t;
}

// xCors[] and yCors[] are coords in x and y in sorted order
// s[] is a list of nodes in increasing y direction
//   if nodes are indexed in the order of increasing x coord
//   i.e., s[i] = s_i as defined in paper
// The points are (xCors[s[i]], yCors[i]) for i=0..d-1
//             or (xCors[i], yCors[si[i]]) for i=0..d-1

Tree flutes_RDP(int d, DATATYPE xCors[], DATATYPE yCors[], int s[], int accuracy)
{
    int i, j, ss;

    for (i = 0; i < d - 1; i++)
    {
        if (xCors[s[i]] == xCors[s[i + 1]] && yCors[i] == yCors[i + 1])
        {
            if (s[i] < s[i + 1])
                ss = s[i + 1];
            else
            {
                ss = s[i];
                s[i] = s[i + 1];
            }
            for (j = i + 2; j < d; j++)
            {
                yCors[j - 1] = yCors[j];
                s[j - 1] = s[j];
            }
            for (j = ss + 1; j < d; j++)
                xCors[j - 1] = xCors[j];
            for (j = 0; j <= d - 2; j++)
                if (s[j] > ss)
                    s[j]--;
            i--;
            d--;
        }
    }
    return flutes_ALLD(d, xCors, yCors, s, accuracy);
}

// For low-degree, i.e., 2 <= d <= D
Tree flutes_LD(int d, DATATYPE xCors[], DATATYPE yCors[], int s[])
{
    int k, pi, i, j;
    struct Csoln *rlist, *bestrlist;
    DATATYPE dd[2 * D - 2]; // 0..D-2 for v, D-1..2*D-3 for h
    DATATYPE minl, sum, l[MPOWV + 1];
    int hflip;
    Tree t;

    t.deg = d;
    t.branch = (Branch *)malloc((2 * d - 2) * sizeof(Branch));
    if (d == 2)
    {
        minl = xCors[1] - xCors[0] + yCors[1] - yCors[0];
        t.branch[0].x = xCors[s[0]];
        t.branch[0].y = yCors[0];
        t.branch[0].n = 1;
        t.branch[1].x = xCors[s[1]];
        t.branch[1].y = yCors[1];
        t.branch[1].n = 1;
    }
    else if (d == 3)
    {
        minl = xCors[2] - xCors[0] + yCors[2] - yCors[0];
        t.branch[0].x = xCors[s[0]];
        t.branch[0].y = yCors[0];
        t.branch[0].n = 3;
        t.branch[1].x = xCors[s[1]];
        t.branch[1].y = yCors[1];
        t.branch[1].n = 3;
        t.branch[2].x = xCors[s[2]];
        t.branch[2].y = yCors[2];
        t.branch[2].n = 3;
        t.branch[3].x = xCors[1];
        t.branch[3].y = yCors[1];
        t.branch[3].n = 3;
    }
    else
    {
        k = 0;
        if (s[0] < s[2])
            k++;
        if (s[1] < s[2])
            k++;

        for (i = 3; i <= d - 1; i++)
        { // p0=0 alwayCors, skip i=1 for symmetry
            pi = s[i];
            for (j = d - 1; j > i; j--)
                if (s[j] < s[i])
                    pi--;
            k = pi + (i + 1) * k;
        }

        if (k < numGRPMap[d])
        { // no horizontal flip
            hflip = 0;
            for (i = 1; i <= d - 3; i++)
            {
                dd[i] = yCors[i + 1] - yCors[i];
                dd[d - 1 + i] = xCors[i + 1] - xCors[i];
            }
        }
        else
        {
            hflip = 1;
            k = 2 * numGRPMap[d] - 1 - k;
            for (i = 1; i <= d - 3; i++)
            {
                dd[i] = yCors[i + 1] - yCors[i];
                dd[d - 1 + i] = xCors[d - 1 - i] - xCors[d - 2 - i];
            }
        }

        minl = l[0] = xCors[d - 1] - xCors[0] + yCors[d - 1] - yCors[0];
        rlist = LUT[d][k];
        for (i = 0; rlist->seg[i] > 0; i++)
            minl += dd[rlist->seg[i]];
        bestrlist = rlist;
        l[1] = minl;
        j = 2;
        while (j <= numSoln[d][k])
        {
            rlist++;
            sum = l[rlist->parent];
            for (i = 0; rlist->seg[i] > 0; i++)
                sum += dd[rlist->seg[i]];
            for (i = 10; rlist->seg[i] > 0; i--)
                sum -= dd[rlist->seg[i]];
            if (sum < minl)
            {
                minl = sum;
                bestrlist = rlist;
            }
            l[j++] = sum;
        }

        t.branch[0].x = xCors[s[0]];
        t.branch[0].y = yCors[0];
        t.branch[1].x = xCors[s[1]];
        t.branch[1].y = yCors[1];
        for (i = 2; i < d - 2; i++)
        {
            t.branch[i].x = xCors[s[i]];
            t.branch[i].y = yCors[i];
            t.branch[i].n = bestrlist->neighbor[i];
        }
        t.branch[d - 2].x = xCors[s[d - 2]];
        t.branch[d - 2].y = yCors[d - 2];
        t.branch[d - 1].x = xCors[s[d - 1]];
        t.branch[d - 1].y = yCors[d - 1];
        if (hflip)
        {
            if (s[1] < s[0])
            {
                t.branch[0].n = bestrlist->neighbor[1];
                t.branch[1].n = bestrlist->neighbor[0];
            }
            else
            {
                t.branch[0].n = bestrlist->neighbor[0];
                t.branch[1].n = bestrlist->neighbor[1];
            }
            if (s[d - 1] < s[d - 2])
            {
                t.branch[d - 2].n = bestrlist->neighbor[d - 1];
                t.branch[d - 1].n = bestrlist->neighbor[d - 2];
            }
            else
            {
                t.branch[d - 2].n = bestrlist->neighbor[d - 2];
                t.branch[d - 1].n = bestrlist->neighbor[d - 1];
            }
            for (i = d; i < 2 * d - 2; i++)
            {
                t.branch[i].x = xCors[d - 1 - bestrlist->rowcol[i - d] % 16];
                t.branch[i].y = yCors[bestrlist->rowcol[i - d] / 16];
                t.branch[i].n = bestrlist->neighbor[i];
            }
        }
        else
        { // !hflip
            if (s[0] < s[1])
            {
                t.branch[0].n = bestrlist->neighbor[1];
                t.branch[1].n = bestrlist->neighbor[0];
            }
            else
            {
                t.branch[0].n = bestrlist->neighbor[0];
                t.branch[1].n = bestrlist->neighbor[1];
            }
            if (s[d - 2] < s[d - 1])
            {
                t.branch[d - 2].n = bestrlist->neighbor[d - 1];
                t.branch[d - 1].n = bestrlist->neighbor[d - 2];
            }
            else
            {
                t.branch[d - 2].n = bestrlist->neighbor[d - 2];
                t.branch[d - 1].n = bestrlist->neighbor[d - 1];
            }
            for (i = d; i < 2 * d - 2; i++)
            {
                t.branch[i].x = xCors[bestrlist->rowcol[i - d] % 16];
                t.branch[i].y = yCors[bestrlist->rowcol[i - d] / 16];
                t.branch[i].n = bestrlist->neighbor[i];
            }
        }
    }
    t.length = minl;

    return t;
}

// For medium-degree, i.e., D+1 <= d
Tree flutes_MD(int d, DATATYPE xCors[], DATATYPE yCors[], int s[], int accuracy)
{
    DATATYPE x1[MAXD], x2[MAXD], y1[MAXD], y2[MAXD];
    int si[MAXD], s1[MAXD], s2[MAXD];
    float score[2 * MAXD], penalty[MAXD], pnlty, dx, dy;
    DATATYPE ll, minl, coord1, coord2;
    int i, r, p, maxbp, bestbp, bp, nbp, ub, lb, n1, n2, nn1, nn2, newaccuracy;
    Tree t, t1, t2, bestt1, bestt2;
    int ms, mins, maxCors, minsi, maxCorsi;
    DATATYPE distx[MAXD], disty[MAXD], xydiff;

    if (s[0] < s[d - 1])
    {
        ms = max(s[0], s[1]);
        for (i = 2; i <= ms; i++)
            ms = max(ms, s[i]);
        if (ms <= d - 3)
        {
            for (i = 0; i <= ms; i++)
            {
                x1[i] = xCors[i];
                y1[i] = yCors[i];
                s1[i] = s[i];
            }
            x1[ms + 1] = xCors[ms];
            y1[ms + 1] = yCors[ms];
            s1[ms + 1] = ms + 1;

            s2[0] = 0;
            for (i = 1; i <= d - 1 - ms; i++)
                s2[i] = s[i + ms] - ms;

            t1 = flutes_LMD(ms + 2, x1, y1, s1, accuracy);
            t2 = flutes_LMD(d - ms, xCors + ms, yCors + ms, s2, accuracy);
            t = dmergetree(t1, t2);
            free(t1.branch);
            free(t2.branch);

            return t;
        }
    }
    else
    { // (s[0] > s[d-1])
        ms = min(s[0], s[1]);
        for (i = 2; i <= d - 1 - ms; i++)
            ms = min(ms, s[i]);
        if (ms >= 2)
        {
            x1[0] = xCors[ms];
            y1[0] = yCors[0];
            s1[0] = s[0] - ms + 1;
            for (i = 1; i <= d - 1 - ms; i++)
            {
                x1[i] = xCors[i + ms - 1];
                y1[i] = yCors[i];
                s1[i] = s[i] - ms + 1;
            }
            x1[d - ms] = xCors[d - 1];
            y1[d - ms] = yCors[d - 1 - ms];
            s1[d - ms] = 0;

            s2[0] = ms;
            for (i = 1; i <= ms; i++)
                s2[i] = s[i + d - 1 - ms];

            t1 = flutes_LMD(d + 1 - ms, x1, y1, s1, accuracy);
            t2 = flutes_LMD(ms + 1, xCors, yCors + d - 1 - ms, s2, accuracy);
            t = dmergetree(t1, t2);
            free(t1.branch);
            free(t2.branch);

            return t;
        }
    }

    // Find inverse si[] of s[]
    for (r = 0; r < d; r++)
        si[s[r]] = r;

    // Determine breaking directions and positions dp[]
    lb = (d - 2 * accuracy + 2) / 4;
    if (lb < 2)
        lb = 2;
    ub = d - 1 - lb;

// Compute scores
#define AA 0.6 // 2.0*BB
#define BB 0.3
    float CC = 7.4 / ((d + 10.) * (d - 3.));
    float DD = 4.8 / (d - 1);

    // Compute penalty[]
    dx = CC * (xCors[d - 2] - xCors[1]);
    dy = CC * (yCors[d - 2] - yCors[1]);
    for (r = d / 2, pnlty = 0; r >= 2; r--, pnlty += dx)
        penalty[r] = pnlty, penalty[d - 1 - r] = pnlty;
    penalty[1] = pnlty, penalty[d - 2] = pnlty;
    penalty[0] = pnlty, penalty[d - 1] = pnlty;
    for (r = d / 2 - 1, pnlty = dy; r >= 2; r--, pnlty += dy)
        penalty[s[r]] += pnlty, penalty[s[d - 1 - r]] += pnlty;
    penalty[s[1]] += pnlty, penalty[s[d - 2]] += pnlty;
    penalty[s[0]] += pnlty, penalty[s[d - 1]] += pnlty;
    xydiff = (xCors[d - 1] - xCors[0]) - (yCors[d - 1] - yCors[0]);
    if (s[0] < s[1])
        mins = s[0], maxCors = s[1];
    else
        mins = s[1], maxCors = s[0];
    if (si[0] < si[1])
        minsi = si[0], maxCorsi = si[1];
    else
        minsi = si[1], maxCorsi = si[0];
    for (r = 2; r <= ub; r++)
    {
        if (s[r] < mins)
            mins = s[r];
        else if (s[r] > maxCors)
            maxCors = s[r];
        distx[r] = xCors[maxCors] - xCors[mins];
        if (si[r] < minsi)
            minsi = si[r];
        else if (si[r] > maxCorsi)
            maxCorsi = si[r];
        disty[r] = yCors[maxCorsi] - yCors[minsi] + xydiff;
    }

    if (s[d - 2] < s[d - 1])
        mins = s[d - 2], maxCors = s[d - 1];
    else
        mins = s[d - 1], maxCors = s[d - 2];
    if (si[d - 2] < si[d - 1])
        minsi = si[d - 2], maxCorsi = si[d - 1];
    else
        minsi = si[d - 1], maxCorsi = si[d - 2];
    for (r = d - 3; r >= lb; r--)
    {
        if (s[r] < mins)
            mins = s[r];
        else if (s[r] > maxCors)
            maxCors = s[r];
        distx[r] += xCors[maxCors] - xCors[mins];
        if (si[r] < minsi)
            minsi = si[r];
        else if (si[r] > maxCorsi)
            maxCorsi = si[r];
        disty[r] += yCors[maxCorsi] - yCors[minsi];
    }

    nbp = 0;
    for (r = lb; r <= ub; r++)
    {
        if (si[r] <= 1)
            score[nbp] = (xCors[r + 1] - xCors[r - 1]) - penalty[r] - AA * (yCors[2] - yCors[1]) - DD * disty[r];
        else if (si[r] >= d - 2)
            score[nbp] = (xCors[r + 1] - xCors[r - 1]) - penalty[r] - AA * (yCors[d - 2] - yCors[d - 3]) - DD * disty[r];
        else
            score[nbp] = (xCors[r + 1] - xCors[r - 1]) - penalty[r] - BB * (yCors[si[r] + 1] - yCors[si[r] - 1]) - DD * disty[r];
        nbp++;

        if (s[r] <= 1)
            score[nbp] = (yCors[r + 1] - yCors[r - 1]) - penalty[s[r]] - AA * (xCors[2] - xCors[1]) - DD * distx[r];
        else if (s[r] >= d - 2)
            score[nbp] = (yCors[r + 1] - yCors[r - 1]) - penalty[s[r]] - AA * (xCors[d - 2] - xCors[d - 3]) - DD * distx[r];
        else
            score[nbp] = (yCors[r + 1] - yCors[r - 1]) - penalty[s[r]] - BB * (xCors[s[r] + 1] - xCors[s[r] - 1]) - DD * distx[r];
        nbp++;
    }

    if (accuracy <= 3)
        newaccuracy = 1;
    else
    {
        newaccuracy = accuracy / 2;
        if (accuracy >= nbp)
            accuracy = nbp - 1;
    }

    minl = (DATATYPE)INT_MAX;
    bestt1.branch = bestt2.branch = NULL;
    for (i = 0; i < accuracy; i++)
    {
        maxbp = 0;
        for (bp = 1; bp < nbp; bp++)
            if (score[maxbp] < score[bp])
                maxbp = bp;
        score[maxbp] = -9e9;

#define BreakPt(bp) ((bp) / 2 + lb)
#define BreakInX(bp) ((bp) % 2 == 0)
        p = BreakPt(maxbp);
        // Breaking in p
        if (BreakInX(maxbp))
        { // break in x
            n1 = n2 = 0;
            for (r = 0; r < d; r++)
            {
                if (s[r] < p)
                {
                    s1[n1] = s[r];
                    y1[n1] = yCors[r];
                    n1++;
                }
                else if (s[r] > p)
                {
                    s2[n2] = s[r] - p;
                    y2[n2] = yCors[r];
                    n2++;
                }
                else
                { // if (s[r] == p)  i.e.,  r = si[p]
                    s1[n1] = p;
                    s2[n2] = 0;
                    y1[n1] = y2[n2] = yCors[r];
                    nn1 = n1;
                    nn2 = n2;
                    n1++;
                    n2++;
                }
            }

            t1 = flutes_LMD(p + 1, xCors, y1, s1, newaccuracy);
            t2 = flutes_LMD(d - p, xCors + p, y2, s2, newaccuracy);
            ll = t1.length + t2.length;
            coord1 = t1.branch[t1.branch[nn1].n].y;
            coord2 = t2.branch[t2.branch[nn2].n].y;
            if (t2.branch[nn2].y > max(coord1, coord2))
                ll -= t2.branch[nn2].y - max(coord1, coord2);
            else if (t2.branch[nn2].y < min(coord1, coord2))
                ll -= min(coord1, coord2) - t2.branch[nn2].y;
        }
        else
        { // if (!BreakInX(maxbp))
            n1 = n2 = 0;
            for (r = 0; r < d; r++)
            {
                if (si[r] < p)
                {
                    s1[si[r]] = n1;
                    x1[n1] = xCors[r];
                    n1++;
                }
                else if (si[r] > p)
                {
                    s2[si[r] - p] = n2;
                    x2[n2] = xCors[r];
                    n2++;
                }
                else
                { // if (si[r] == p)  i.e.,  r = s[p]
                    s1[p] = n1;
                    s2[0] = n2;
                    x1[n1] = x2[n2] = xCors[r];
                    n1++;
                    n2++;
                }
            }

            t1 = flutes_LMD(p + 1, x1, yCors, s1, newaccuracy);
            t2 = flutes_LMD(d - p, x2, yCors + p, s2, newaccuracy);
            ll = t1.length + t2.length;
            coord1 = t1.branch[t1.branch[p].n].x;
            coord2 = t2.branch[t2.branch[0].n].x;
            if (t2.branch[0].x > max(coord1, coord2))
                ll -= t2.branch[0].x - max(coord1, coord2);
            else if (t2.branch[0].x < min(coord1, coord2))
                ll -= min(coord1, coord2) - t2.branch[0].x;
        }
        if (minl > ll)
        {
            minl = ll;
            free(bestt1.branch);
            free(bestt2.branch);
            bestt1 = t1;
            bestt2 = t2;
            bestbp = maxbp;
        }
        else
        {
            free(t1.branch);
            free(t2.branch);
        }
    }

#if LOCAL_REFINEMENT == 1
    if (BreakInX(bestbp))
    {
        t = hmergetree(bestt1, bestt2, s);
        local_refinement(&t, si[BreakPt(bestbp)]);
    }
    else
    {
        t = vmergetree(bestt1, bestt2);
        local_refinement(&t, BreakPt(bestbp));
    }
#else
    if (BreakInX(bestbp))
    {
        t = hmergetree(bestt1, bestt2, s);
    }
    else
    {
        t = vmergetree(bestt1, bestt2);
    }
#endif

    free(bestt1.branch);
    free(bestt2.branch);

    return t;
}

Tree dmergetree(Tree t1, Tree t2)
{
    int i, d, prev, curr, next, offset1, offset2;
    Tree t;

    t.deg = d = t1.deg + t2.deg - 2;
    t.length = t1.length + t2.length;
    t.branch = (Branch *)malloc((2 * d - 2) * sizeof(Branch));
    offset1 = t2.deg - 2;
    offset2 = 2 * t1.deg - 4;

    for (i = 0; i <= t1.deg - 2; i++)
    {
        t.branch[i].x = t1.branch[i].x;
        t.branch[i].y = t1.branch[i].y;
        t.branch[i].n = t1.branch[i].n + offset1;
    }
    for (i = t1.deg - 1; i <= d - 1; i++)
    {
        t.branch[i].x = t2.branch[i - t1.deg + 2].x;
        t.branch[i].y = t2.branch[i - t1.deg + 2].y;
        t.branch[i].n = t2.branch[i - t1.deg + 2].n + offset2;
    }
    for (i = d; i <= d + t1.deg - 3; i++)
    {
        t.branch[i].x = t1.branch[i - offset1].x;
        t.branch[i].y = t1.branch[i - offset1].y;
        t.branch[i].n = t1.branch[i - offset1].n + offset1;
    }
    for (i = d + t1.deg - 2; i <= 2 * d - 3; i++)
    {
        t.branch[i].x = t2.branch[i - offset2].x;
        t.branch[i].y = t2.branch[i - offset2].y;
        t.branch[i].n = t2.branch[i - offset2].n + offset2;
    }

    prev = t2.branch[0].n + offset2;
    curr = t1.branch[t1.deg - 1].n + offset1;
    next = t.branch[curr].n;
    while (curr != next)
    {
        t.branch[curr].n = prev;
        prev = curr;
        curr = next;
        next = t.branch[curr].n;
    }
    t.branch[curr].n = prev;

    return t;
}

Tree hmergetree(Tree t1, Tree t2, int s[])
{
    int i, prev, curr, next, extra, offset1, offset2;
    int p, ii, n1, n2, nn1, nn2;
    DATATYPE coord1, coord2;
    Tree t;

    t.deg = t1.deg + t2.deg - 1;
    t.length = t1.length + t2.length;
    t.branch = (Branch *)malloc((2 * t.deg - 2) * sizeof(Branch));
    offset1 = t2.deg - 1;
    offset2 = 2 * t1.deg - 3;

    p = t1.deg - 1;
    n1 = n2 = 0;
    for (i = 0; i < t.deg; i++)
    {
        if (s[i] < p)
        {
            t.branch[i].x = t1.branch[n1].x;
            t.branch[i].y = t1.branch[n1].y;
            t.branch[i].n = t1.branch[n1].n + offset1;
            n1++;
        }
        else if (s[i] > p)
        {
            t.branch[i].x = t2.branch[n2].x;
            t.branch[i].y = t2.branch[n2].y;
            t.branch[i].n = t2.branch[n2].n + offset2;
            n2++;
        }
        else
        {
            t.branch[i].x = t2.branch[n2].x;
            t.branch[i].y = t2.branch[n2].y;
            t.branch[i].n = t2.branch[n2].n + offset2;
            nn1 = n1;
            nn2 = n2;
            ii = i;
            n1++;
            n2++;
        }
    }
    for (i = t.deg; i <= t.deg + t1.deg - 3; i++)
    {
        t.branch[i].x = t1.branch[i - offset1].x;
        t.branch[i].y = t1.branch[i - offset1].y;
        t.branch[i].n = t1.branch[i - offset1].n + offset1;
    }
    for (i = t.deg + t1.deg - 2; i <= 2 * t.deg - 4; i++)
    {
        t.branch[i].x = t2.branch[i - offset2].x;
        t.branch[i].y = t2.branch[i - offset2].y;
        t.branch[i].n = t2.branch[i - offset2].n + offset2;
    }
    extra = 2 * t.deg - 3;
    coord1 = t1.branch[t1.branch[nn1].n].y;
    coord2 = t2.branch[t2.branch[nn2].n].y;
    if (t2.branch[nn2].y > max(coord1, coord2))
    {
        t.branch[extra].y = max(coord1, coord2);
        t.length -= t2.branch[nn2].y - t.branch[extra].y;
    }
    else if (t2.branch[nn2].y < min(coord1, coord2))
    {
        t.branch[extra].y = min(coord1, coord2);
        t.length -= t.branch[extra].y - t2.branch[nn2].y;
    }
    else
        t.branch[extra].y = t2.branch[nn2].y;
    t.branch[extra].x = t2.branch[nn2].x;
    t.branch[extra].n = t.branch[ii].n;
    t.branch[ii].n = extra;

    prev = extra;
    curr = t1.branch[nn1].n + offset1;
    next = t.branch[curr].n;
    while (curr != next)
    {
        t.branch[curr].n = prev;
        prev = curr;
        curr = next;
        next = t.branch[curr].n;
    }
    t.branch[curr].n = prev;

    return t;
}

Tree vmergetree(Tree t1, Tree t2)
{
    int i, prev, curr, next, extra, offset1, offset2;
    DATATYPE coord1, coord2;
    Tree t;

    t.deg = t1.deg + t2.deg - 1;
    t.length = t1.length + t2.length;
    t.branch = (Branch *)malloc((2 * t.deg - 2) * sizeof(Branch));
    offset1 = t2.deg - 1;
    offset2 = 2 * t1.deg - 3;

    for (i = 0; i <= t1.deg - 2; i++)
    {
        t.branch[i].x = t1.branch[i].x;
        t.branch[i].y = t1.branch[i].y;
        t.branch[i].n = t1.branch[i].n + offset1;
    }
    for (i = t1.deg - 1; i <= t.deg - 1; i++)
    {
        t.branch[i].x = t2.branch[i - t1.deg + 1].x;
        t.branch[i].y = t2.branch[i - t1.deg + 1].y;
        t.branch[i].n = t2.branch[i - t1.deg + 1].n + offset2;
    }
    for (i = t.deg; i <= t.deg + t1.deg - 3; i++)
    {
        t.branch[i].x = t1.branch[i - offset1].x;
        t.branch[i].y = t1.branch[i - offset1].y;
        t.branch[i].n = t1.branch[i - offset1].n + offset1;
    }
    for (i = t.deg + t1.deg - 2; i <= 2 * t.deg - 4; i++)
    {
        t.branch[i].x = t2.branch[i - offset2].x;
        t.branch[i].y = t2.branch[i - offset2].y;
        t.branch[i].n = t2.branch[i - offset2].n + offset2;
    }
    extra = 2 * t.deg - 3;
    coord1 = t1.branch[t1.branch[t1.deg - 1].n].x;
    coord2 = t2.branch[t2.branch[0].n].x;
    if (t2.branch[0].x > max(coord1, coord2))
    {
        t.branch[extra].x = max(coord1, coord2);
        t.length -= t2.branch[0].x - t.branch[extra].x;
    }
    else if (t2.branch[0].x < min(coord1, coord2))
    {
        t.branch[extra].x = min(coord1, coord2);
        t.length -= t.branch[extra].x - t2.branch[0].x;
    }
    else
        t.branch[extra].x = t2.branch[0].x;
    t.branch[extra].y = t2.branch[0].y;
    t.branch[extra].n = t.branch[t1.deg - 1].n;
    t.branch[t1.deg - 1].n = extra;

    prev = extra;
    curr = t1.branch[t1.deg - 1].n + offset1;
    next = t.branch[curr].n;
    while (curr != next)
    {
        t.branch[curr].n = prev;
        prev = curr;
        curr = next;
        next = t.branch[curr].n;
    }
    t.branch[curr].n = prev;

    return t;
}

void local_refinement(Tree *tp, int p)
{
    int d, dd, i, ii, j, prev, curr, next, root;
    int SteinerPin[2 * MAXD], index[2 * MAXD];
    DATATYPE x[MAXD], xCors[D], yCors[D];
    int ss[D];
    Tree tt;

    d = tp->deg;
    root = tp->branch[p].n;

    prev = root;
    curr = tp->branch[prev].n;
    next = tp->branch[curr].n;
    while (curr != next)
    {
        tp->branch[curr].n = prev;
        prev = curr;
        curr = next;
        next = tp->branch[curr].n;
    }
    tp->branch[curr].n = prev;
    tp->branch[root].n = root;

    for (i = d; i <= 2 * d - 3; i++)
        SteinerPin[i] = -1;
    for (i = 0; i < d; i++)
    {
        next = tp->branch[i].n;
        if (tp->branch[i].x == tp->branch[next].x &&
            tp->branch[i].y == tp->branch[next].y)
            SteinerPin[next] = i;
    }
    SteinerPin[root] = p;

    dd = 0;
    for (i = 0; i < d; i++)
    {
        curr = tp->branch[i].n;
        if (SteinerPin[curr] == i)
            curr = tp->branch[curr].n;
        while (SteinerPin[curr] < 0)
            curr = tp->branch[curr].n;
        if (curr == root)
        {
            x[dd] = tp->branch[i].x;
            if (SteinerPin[tp->branch[i].n] == i && tp->branch[i].n != root)
                index[dd++] = tp->branch[i].n; 
            else
                index[dd++] = i;
        }
    }

    if (4 <= dd && dd <= D)
    {
        ii = dd;
        for (i = 0; i < dd; i++){
            curr = tp->branch[index[i]].n;
            while (SteinerPin[curr] < 0){
                index[ii++] = curr;
                SteinerPin[curr] = INT_MAX;
                curr = tp->branch[curr].n;
            }
        }
        index[ii] = root;

        for (ii = 0; ii < dd; ii++){
            ss[ii] = 0;
            for (j = 0; j < ii; j++){
                if (x[j] < x[ii]){
                    ss[ii]++;
                }
            }
            for (j = ii + 1; j < dd; j++)
                if (x[j] <= x[ii])
                    ss[ii]++;
            xCors[ss[ii]] = x[ii];
            yCors[ii] = tp->branch[index[ii]].y;
        }

        tt = flutes_LD(dd, xCors, yCors, ss);

        
        tp->length += tt.length;
        for (ii = 0; ii < 2 * dd - 3; ii++){
            i = index[ii];
            j = tp->branch[i].n;
            tp->length -= ADIFF(tp->branch[i].x, tp->branch[j].x) + ADIFF(tp->branch[i].y, tp->branch[j].y);
        }

        
        for (ii = 0; ii < dd; ii++){
            tp->branch[index[ii]].n = index[tt.branch[ii].n];
        }
        for (; ii <= 2 * dd - 3; ii++)
        {
            tp->branch[index[ii]].x = tt.branch[ii].x;
            tp->branch[index[ii]].y = tt.branch[ii].y;
            tp->branch[index[ii]].n = index[tt.branch[ii].n];
        }
        free(tt.branch);
    }

    return;
}

DATATYPE wirelength(Tree t){
    int i, j;
    DATATYPE l = 0;

    for (i = 0; i < 2 * t.deg - 2; i++)
    {
        j = t.branch[i].n;
        l += ADIFF(t.branch[i].x, t.branch[j].x) + ADIFF(t.branch[i].y, t.branch[j].y);
    }

    return l;
}

void printTree(Tree t){
    int i;

    for (i = 0; i < t.deg; i++){
        printf(" %-2d:  x=%4g  y=%4g  e=%d\n", i, (float)t.branch[i].x, (float)t.branch[i].y, t.branch[i].n);
    }
    for (i = t.deg; i < 2 * t.deg - 2; i++){
        printf("s%-2d:  x=%4g  y=%4g  e=%d\n", i, (float)t.branch[i].x, (float)t.branch[i].y, t.branch[i].n);
    }
    std::cout << std::endl;
}


void plotTree(Tree t){
    int i;

    for (i = 0; i < 2 * t.deg - 2; i ++){
        printf("%d %d\n", t.branch[i].x, t.branch[i].y);
        printf("%d %d\n\n", t.branch[t.branch[i].n].x,t.branch[t.branch[i].n].y);
    }
}




pairArray treeToPairArray(Tree t){
    pairArray res;
    int from_index, to_index;
    for (int j = 0; j < 2 * t.deg - 2; j++)
    {
        from_index = j;
        to_index = t.branch[j].n;
        res.x.push_back(t.branch[from_index].x);
        res.x.push_back(t.branch[to_index].x);
        res.y.push_back(t.branch[from_index].y);
        res.y.push_back(t.branch[to_index].y);
    }
    return res;
}


PYBIND11_MODULE(_flute, m) {
    pybind11::class_<Tree>(m, "Tree");
    pybind11::class_<Csoln>(m, "Csoln");
    pybind11::class_<pairArray>(m, "PairArray")
    .def_property_readonly("x", [](pairArray &self){return self.x;})
    .def_property_readonly("y", [](pairArray &self){return self.y;});
    m.def("getDegree", &Tree::getDegree);
    m.def("flute", [](int degree, std::vector<int> a, std::vector<int> b, int accuracyuracy) {
            return flute(degree, a.data(), b.data(), accuracyuracy);
        }, "basic flute");
    m.def("flute_wl", [](int degree, std::vector<int> a, std::vector<int> b, int accuracyuracy) {
            return flute_wl(degree, a.data(), b.data(), accuracyuracy);
        }, "basic flute");
    m.def("printtree", &printTree, "printtree");
    m.def("plottree", &plotTree, "plottree");
    m.def("readLUT", &readLUT, "readLUT");
    m.def("treeToPairArray", &treeToPairArray, "treeToPairArray");
    m.def("wirelength", &wirelength, "caculate the wirelength");
}
