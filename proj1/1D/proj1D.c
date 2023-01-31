/* proj1D.c
 * Author: Nate Koch
 * Date: 1/30/23
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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
    double  X[3];
    double  Y[3];
    double  Z[3];
    double  colors[3][3];
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
    double *zbuffer;
} Image;

void WritePNM(Image *img);

void ColorPixel(Image *img, int row, int col, double R, double G, double B, double Z);

TriangleList* Get3DTriangles(void);

void RasterizeTriangle(Triangle *t, Image *img);

int main(void)
{
    // Create the Image struct
    Image *img = NULL;

    img = malloc(sizeof(Image));
    img->x = 1000;
    img->y = 1000;
    img->data = malloc(sizeof(Pixel) * img->x * img->y);
    memset(img->data, 0, sizeof(Pixel) * img->x * img->y); 
    img->zbuffer = malloc(sizeof(double) * img->x * img->y);
    for (int i = 0; i < img->x*img->y; i++) img->zbuffer[i] = -1;

    // Rasterize the triangles
    TriangleList* tl = NULL;

    tl = Get3DTriangles();
    
    
    for (int i = 0; i < tl->numTriangles; i++)
        RasterizeTriangle(tl->triangles+i, img);

    free(tl->triangles);
    free(tl);

    // Format and create the PNM file
    WritePNM(img);    
    
    free(img->data);
    free(img->zbuffer);
    free(img);

    return 0; 
}

void WritePNM(Image *img) 
{
    FILE *pnm = fopen("proj1D_out.pnm", "wb");

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

void ColorPixel(Image *img, int row, int col, double R, double G, double B, double Z)
{
    // check if row or column is off the screen
    if (row < 0 || row >= img->y || 
        col < 0 || col >= img->x) 
        return;

    // set origin bottom left
    int newRow = img->y - row - 1;

    // translate 2D to 1D
    int index = newRow * img->x + col;
    
    if (img->zbuffer[index] > Z) return;
    img->zbuffer[index] = Z; 
      
    // color the pixel
    img->data[index].R = C441(R*255);
    img->data[index].G = C441(G*255);
    img->data[index].B = C441(B*255);
}

char *
ReadTuple3(char *tmp, double *v1, double *v2, double *v3)
{
    tmp++; /* left paren */
    *v1 = atof(tmp);
    while (*tmp != ',')
       tmp++;
    tmp += 2; // comma+space
    *v2 = atof(tmp);
    while (*tmp != ',')
       tmp++;
    tmp += 2; // comma+space
    *v3 = atof(tmp);
    while (*tmp != ')')
       tmp++;
    tmp++; /* right paren */
    return tmp;
}

TriangleList *
Get3DTriangles()
{
   FILE *f = fopen("tris_w_r_rgb.txt", "r");
   if (f == NULL)
   {
       fprintf(stderr, "You must place the tris_w_r_rgb.txt file in the current directory.\n");
       exit(EXIT_FAILURE);
   }
   fseek(f, 0, SEEK_END);
   int numBytes = ftell(f);
   fseek(f, 0, SEEK_SET);
   if (numBytes != 13488634)
   {
       fprintf(stderr, "Your tris_w_r_rgb.txt file is corrupted.  It should be 13488634 bytes, but you have %d.\n", numBytes);
       exit(EXIT_FAILURE);
   }

   char *buffer = (char *) malloc(numBytes);
   if (buffer == NULL)
   {
       fprintf(stderr, "Unable to allocate enough memory to load file.\n");
       exit(EXIT_FAILURE);
   }
   
   fread(buffer, sizeof(char), numBytes, f);

   char *tmp = buffer;
   int numTriangles = atoi(tmp);
   while (*tmp != '\n')
       tmp++;
   tmp++;
 
   if (numTriangles != 42281)
   {
       fprintf(stderr, "Issue with reading file -- can't establish number of triangles.\n");
       exit(EXIT_FAILURE);
   }

   TriangleList *tl = (TriangleList *) malloc(sizeof(TriangleList));
   tl->numTriangles = numTriangles;
   tl->triangles = (Triangle *) malloc(sizeof(Triangle)*tl->numTriangles);

   for (int i = 0 ; i < tl->numTriangles ; i++)
   {
       double x1, y1, z1, x2, y2, z2, x3, y3, z3;
       double r[3], g[3], b[3];
/*
 * Weird: sscanf has a terrible implementation for large strings.
 * When I did the code below, it did not finish after 45 minutes.
 * Reading up on the topic, it sounds like it is a known issue that
 * sscanf fails here.  Stunningly, fscanf would have been faster.
 *     sscanf(tmp, "(%lf, %lf), (%lf, %lf), (%lf, %lf) = (%d, %d, %d)\n%n",
 *              &x1, &y1, &x2, &y2, &x3, &y3, &r, &g, &b, &numRead);
 *
 *  So, instead, do it all with atof/atoi and advancing through the buffer manually...
 */
       tmp = ReadTuple3(tmp, &x1, &y1, &z1);
       tmp += 3; /* space+equal+space */
       tmp = ReadTuple3(tmp, r+0, g+0, b+0);
       tmp += 2; /* comma+space */
       tmp = ReadTuple3(tmp, &x2, &y2, &z2);
       tmp += 3; /* space+equal+space */
       tmp = ReadTuple3(tmp, r+1, g+1, b+1);
       tmp += 2; /* comma+space */
       tmp = ReadTuple3(tmp, &x3, &y3, &z3);
       tmp += 3; /* space+equal+space */
       tmp = ReadTuple3(tmp, r+2, g+2, b+2);
       tmp++;    /* newline */

       tl->triangles[i].X[0] = x1;
       tl->triangles[i].X[1] = x2;
       tl->triangles[i].X[2] = x3;
       tl->triangles[i].Y[0] = y1;
       tl->triangles[i].Y[1] = y2;
       tl->triangles[i].Y[2] = y3;
       tl->triangles[i].Z[0] = z1;
       tl->triangles[i].Z[1] = z2;
       tl->triangles[i].Z[2] = z3;
       tl->triangles[i].colors[0][0] = r[0];
       tl->triangles[i].colors[0][1] = g[0];
       tl->triangles[i].colors[0][2] = b[0];
       tl->triangles[i].colors[1][0] = r[1];
       tl->triangles[i].colors[1][1] = g[1];
       tl->triangles[i].colors[1][2] = b[1];
       tl->triangles[i].colors[2][0] = r[2];
       tl->triangles[i].colors[2][1] = g[2];
       tl->triangles[i].colors[2][2] = b[2];
   }

   free(buffer);
   return tl;
}

void RasterizeTriangle(Triangle *t, Image *img)
{
    for (int loop = 0; loop < 2; loop++)
    {
        int workOnTopPart = (loop == 0 ? 1 : 0);
        int topOrBottom, left, right, middle;
        int swap;
        double rowMin, rowMax;
        double leftEnd, rightEnd, leftSlope, rightSlope, bLeft, bRight;
        double tPropLeft, tPropRight, tPropPixel;
        double zLeftEnd, zRightEnd, zPixel;
        double colorLeftEnd[3], colorRightEnd[3], colorPixel[3];
        int leftIntercept, rightIntercept, topOrBottomScanline, middleScanline;

        // determine top or bottom
        topOrBottom = 0;
        for (int i = 1; i < 3; i++)
        {
            if (workOnTopPart)
            {
                if (t->Y[topOrBottom] < t->Y[i]) topOrBottom = i;
            }
            else
            {
                if (t->Y[topOrBottom] > t->Y[i]) topOrBottom = i;
            }   
        }

        // determine left and right
        left = (topOrBottom == 0 ? 1 : 0);
        for (int i = 0; i < 3; i++) 
        {   
            if (i == topOrBottom) continue; 
            if (t->X[i] < t->X[left])
            {
                // set i to be left
                swap = left;
                left = i;
                right = swap;
            }
            else
            {
                right = i;
            }
        }

        // determine middle
        if (workOnTopPart) middle = (t->Y[left] > t->Y[right] ? left : right);
        else middle = (t->Y[left] < t->Y[right] ? left : right);

        // Determine rowMin & rowMax
        if (workOnTopPart)
        {
            rowMin = t->Y[middle];
            rowMax = t->Y[topOrBottom];
        }
        else 
        {
            rowMin = t->Y[topOrBottom];
            rowMax = t->Y[middle];
        }
        rowMin = C441(rowMin);
        rowMax = F441(rowMax);

        for (int r = rowMin; r <= rowMax; r++)
        {   

            // left side
            if (t->Y[left] == t->Y[topOrBottom]) continue;
            leftSlope   =   (t->Y[left] - t->Y[topOrBottom]) / 
                            (t->X[left] - t->X[topOrBottom]); 
            bLeft       =   (t->Y[topOrBottom] - leftSlope * t->X[topOrBottom]);
            leftEnd     =   (r - bLeft) / leftSlope;
            if (isnan(leftEnd)) leftEnd = t->X[left];

            // right side
            if (t->Y[right] == t->Y[topOrBottom]) continue;
            rightSlope  =   (t->Y[right] - t->Y[topOrBottom]) /
                            (t->X[right] - t->X[topOrBottom]);
            bRight      =   (t->Y[topOrBottom] - rightSlope * t->X[topOrBottom]);
            rightEnd    =   (r - bRight) / rightSlope;
            if (isnan(rightEnd)) rightEnd = t->X[right];

            if (leftEnd > rightEnd) 
            {
                double temp = leftEnd;
                leftEnd = rightEnd;
                rightEnd = temp;

                swap = left;
                left = right;
                right = swap;
            }

            leftIntercept = C441(leftEnd);
            rightIntercept = F441(rightEnd);

            /* validity checks */
            // check if bottom scanline bigger than top
            topOrBottomScanline = F441(t->Y[topOrBottom]);
            middleScanline = C441(t->Y[middle]);

            if ((workOnTopPart) && (r > topOrBottomScanline || r < middleScanline)) continue;
            if ((!workOnTopPart) && (r > middleScanline || r < topOrBottomScanline)) continue;

            // check intercepts
            //if (leftIntercept > rightIntercept) printf("r: %d", r);
           
            // interpolate z at left and right ends
            // F(X) = F(A) + t*(F(B) - F(A))
            // t = (X-A)/(B-A)
            tPropLeft = (r - t->Y[left]) / (t->Y[topOrBottom] - t->Y[left]);
            zLeftEnd = t->Z[left] + tPropLeft * (t->Z[topOrBottom] - t->Z[left]);

            tPropRight = (r - t->Y[right]) / (t->Y[topOrBottom] - t->Y[right]);
            zRightEnd = t->Z[right] + tPropRight * (t->Z[topOrBottom] - t->Z[right]);

            // interpolate color at left and right ends 
            colorLeftEnd[0] = t->colors[left][0] + tPropLeft * (t->colors[topOrBottom][0] - t->colors[left][0]);
            colorLeftEnd[1] = t->colors[left][1] + tPropLeft * (t->colors[topOrBottom][1] - t->colors[left][1]);
            colorLeftEnd[2] = t->colors[left][2] + tPropLeft * (t->colors[topOrBottom][2] - t->colors[left][2]);

            colorRightEnd[0] = t->colors[right][0] + tPropRight * (t->colors[topOrBottom][0] - t->colors[right][0]);  
            colorRightEnd[1] = t->colors[right][1] + tPropRight * (t->colors[topOrBottom][1] - t->colors[right][1]);
            colorRightEnd[2] = t->colors[right][2] + tPropRight * (t->colors[topOrBottom][2] - t->colors[right][2]);
           
            for (int c = leftIntercept; c <= rightIntercept; c++)
            { 
                tPropPixel = (c - leftEnd) / (rightEnd - leftEnd);
                if (rightEnd != leftEnd)
                {
                    zPixel = zLeftEnd + tPropPixel * (zRightEnd - zLeftEnd);
                    colorPixel[0] = colorLeftEnd[0] + tPropPixel * (colorRightEnd[0] - colorLeftEnd[0]);
                    colorPixel[1] = colorLeftEnd[1] + tPropPixel * (colorRightEnd[1] - colorLeftEnd[1]);
                    colorPixel[2] = colorLeftEnd[2] + tPropPixel * (colorRightEnd[2] - colorLeftEnd[2]);
                }
                else
                {
                    zPixel = zRightEnd;
                    colorPixel[0] = colorRightEnd[0];
                    colorPixel[1] = colorRightEnd[1];
                    colorPixel[2] = colorRightEnd[2];
                }

                ColorPixel(img, r, c, colorPixel[0], colorPixel[1], colorPixel[2], zPixel);
            }
        }
    }
}
