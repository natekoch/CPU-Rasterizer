/* proj1B.c
 * Author: Nate Koch
 * Date: 1/19/23
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG 0

double C441(double f)
{
    return ceil(f-0.00001);
}

double F441(double f)
{
    return floor(f+0.00001);
}

typedef struct
{
   double         X[3];
   double         Y[3];
   unsigned char color[3];
} Triangle;

typedef struct
{
   int numTriangles;
   Triangle *triangles;
} TriangleList;

typedef struct
{
    unsigned char R, G, B;
} Pixel;

typedef struct 
{
    int x, y;
    Pixel *data;
} Image;

void WritePNM(Image *img);

void AssignPixels(Image *img, int colMin, int colMax, int rowMin, int rowMax, 
                    unsigned char R, unsigned char G, unsigned char B);

TriangleList* GetTriangles(void);

void RasterizeGoingUpTriangle(Triangle *t, Image *img); 

int main(void)
{
    // Create the Image struct
    Image *img = NULL;

    img = malloc(sizeof(Image));
    img->x = 1000;
    img->y = 1000;
    img->data = malloc(sizeof(Pixel) * img->x * img->y);

    TriangleList* tl = NULL;

    tl = GetTriangles();

    for (int i = 0; i < tl->numTriangles; i++)
        RasterizeGoingUpTriangle(tl->triangles+i, img);
    //RasterizeGoingUpTriangle(tl->triangles+99, img); 

    // Format and create the PNM file
    WritePNM(img);

    free(tl->triangles);
    free(tl);
    
    free(img->data);
    free(img);

    return 0; 
}

void WritePNM(Image *img) 
{
    FILE *pnm = fopen("proj1B_out.pnm", "w");

    if (pnm == NULL) 
    {
        printf("unable to create pnm file.\n");
        return;
    }

    fprintf(pnm, "P6\n");
    fprintf(pnm, "%d %d\n", img->x, img->y);
    fprintf(pnm, "255\n");

    fwrite(img->data, sizeof(Pixel), img->x * img->y, pnm);
    fclose(pnm);
}

void AssignPixels(Image *img, int colMin, int colMax, int rowMin, int rowMax, 
                    unsigned char R, unsigned char G, unsigned char B)
{
    int index;
    for (int row = rowMin; row < rowMax; row++)
    {
        if (rowMin < 0) rowMin = 0;
        if (rowMax > img->y) rowMax = img->y;
        if (colMin < 0) colMin = 0;
        if (colMax > img->x) colMax = img->x;

        int newRow = img->y - row - 1;
        for (int col = colMin; col < colMax; col++) 
        {
            index = newRow*img->x+col;
            img->data[index].R = R;
            img->data[index].G = G;
            img->data[index].B = B;
        }
    }
}

TriangleList* GetTriangles(void)
{
   TriangleList *tl = (TriangleList *) malloc(sizeof(TriangleList));
   tl->numTriangles = 100;
   tl->triangles = (Triangle *) malloc(sizeof(Triangle)*tl->numTriangles);

   unsigned char colors[6][3] = { {255,128,0}, {255, 0, 127}, {0,204,204}, 
                                  {76,153,0}, {255, 204, 204}, {204, 204, 0}};
   for (int i = 0 ; i < 100 ; i++)
   {
       int idxI = i%10;
       int posI = idxI*100;
       int idxJ = i/10;
       int posJ = idxJ*100;
       int firstPt = (i%3);
       tl->triangles[i].X[firstPt] = posI;
       if (i == 50)
           tl->triangles[i].X[firstPt] = -10;
       tl->triangles[i].Y[firstPt] = posJ+10*(idxJ+1);
       tl->triangles[i].X[(firstPt+1)%3] = posI+105;
       tl->triangles[i].Y[(firstPt+1)%3] = posJ;
       tl->triangles[i].X[(firstPt+2)%3] = posI+i;
       tl->triangles[i].Y[(firstPt+2)%3] = posJ;
       if (i == 95)
          tl->triangles[i].Y[firstPt] = 1050;
       tl->triangles[i].color[0] = colors[i%6][0];
       tl->triangles[i].color[1] = colors[i%6][1];
       tl->triangles[i].color[2] = colors[i%6][2];
   }

   return tl;
}

void RasterizeGoingUpTriangle(Triangle *t, Image *img)
{
    if (DEBUG) printf("Triangle: (%f, %f), (%f, %f), (%f, %f)\n", 
            t->X[0], t->Y[0], t->X[1], t->Y[1], t->X[2], t->Y[2]);
    
    int top, left, right;

    // determine top
    top = 0;
    for (int i = 1; i < 3; i++)
    {
        if (t->Y[top] < t->Y[i]) top = i;
    }

    // determine left and right
    int no_val = 1;
    for (int i = 0; i < 3; i++) 
    {
        if (i == top) continue;
        if (no_val) 
        {
            left = i;
            no_val = 0;
            continue;
        }
        if (t->X[left] > t->X[i])
        {
            int temp = left;
            left = i;
            right = temp;
        }
        else
        {
            right = i;
        }
    }

    if (DEBUG) printf("Identified: top = %d, bottom left = %d, bottom right = %d\n",
                top, left, right);

    // Determine rowMin & rowMax
    double rowMin, rowMax;
    rowMin = rowMax = t->Y[0];
    for (int i = 1; i < 3; i++)
    {
        if (rowMin > t->Y[i]) rowMin = t->Y[i];
        if (rowMax < t->Y[i]) rowMax = t->Y[i];
    }
    rowMin = C441(rowMin);
    if (rowMin == -0.0) rowMin += 0.0; // probably not totally necessary
    rowMax = F441(rowMax);

    if (DEBUG) printf("Scanlines go from %.0f to %.0f\n", rowMin, rowMax);

    // TODO add check for out of bounds point
    for (int r = rowMin; r <= rowMax; r++)
    {
        double leftEnd, rightEnd, leftSlope, rightSlope, bLeft, bRight;
        
        // left side
        leftSlope = (t->Y[left] - t->Y[top]) / (t->X[left] - t->X[top]); 
        bLeft = (t->Y[top] - leftSlope * t->X[top]);
        leftEnd = (r - bLeft) / leftSlope;

        // right side
        rightSlope = (t->Y[right] - t->Y[top]) / (t->X[right] - t->X[top]);
        bRight = (t->Y[top] - rightSlope * t->X[top]);
        rightEnd = (r - bRight) / rightSlope;

        int c = C441(leftEnd);
        int floor = F441(rightEnd);
        
        if (DEBUG) printf("Scanline %d: intercepts go from %d to %d\n",
                            r, c, floor);
        
        for (; c <= floor; c++)
        {
            if (c < 0) c = 0;
            if (r < 0) r = 0;
            if (floor > img->x) floor = img->x-1;
            if (c > img->y) c = img->y-1;
            AssignPixels(img, c, floor+1, r, r+1, 
                            t->color[0], t->color[1], t->color[2]);
        }
    }
}