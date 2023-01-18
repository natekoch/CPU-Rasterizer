/* proj1B.c
 * Author: Nate Koch
 * Date: 1/19/23
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

    //for (int i = 0 ; i < tl->numTriangles ; i++)
        //RasterizeGoingUpTriangle(tl->triangles+i, &img);
    for (int i = 0 ; i < 1 ; i++)
        RasterizeGoingUpTriangle(tl->triangles+i, img); 

    // Format and create the PNM file
    WritePNM(img);

    free(img->data);
    free(img);

    return 0; 
}

void WritePNM(Image *img) 
{
    FILE *pnm = fopen("proj1A_out.pnm", "w");

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
    for (int row = rowMin; row < rowMax; row++)
    {
        for (int col = colMin; col < colMax; col++) 
        {
            img->data[row*img->x+col].R = R;
            img->data[row*img->x+col].G = G;
            img->data[row*img->x+col].B = B;
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