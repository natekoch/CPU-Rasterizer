/* proj1C.c
 * Author: Nate Koch
 * Date: 1/24/23
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

void swap(int v1, int v2);

void WritePNM(Image *img);

void ColorPixel(Image *img, int row, int col, 
                    unsigned char R, unsigned char G, unsigned char B);

TriangleList* GetTriangles(int small_read);

void RasterizeTriangle(Triangle *t, Image *img);

int main(void)
{
    // Create the Image struct
    Image *img = NULL;

    img = malloc(sizeof(Image));
    img->x = 1786;
    img->y = 1344;
    img->data = malloc(sizeof(Pixel) * img->x * img->y);

    // Create the duck
    TriangleList* tl = NULL;

    tl = GetTriangles(0);

    for (int i = 0; i < tl->numTriangles; i++)
        RasterizeTriangle(tl->triangles+i, img);

    free(tl->triangles);
    free(tl);

    // Format and create the PNM file
    WritePNM(img);    
    
    free(img->data);
    free(img);

    return 0; 
}

void swap(int v1, int v2)
{
    int temp = v1;
    v1 = v2;
    v2 = temp;
}

void WritePNM(Image *img) 
{
    FILE *pnm = fopen("proj1C_out.pnm", "wb");

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

void ColorPixel(Image *img, int row, int col, 
                    unsigned char R, unsigned char G, unsigned char B)
{
    // check if row or column is off the screen
    if (row < 0 || row >= img->y || 
        col < 0 || col >= img->x) 
        return;

    // set origin bottom left
    int newRow = img->y - row - 1;

    // translate 2D to 1D
    int index = newRow * img->x + col;
    
    // color the pixel
    img->data[index].R = R;
    img->data[index].G = G;
    img->data[index].B = B;
}

TriangleList * GetTriangles(int small_read)
{
   FILE *f = fopen("tris.txt", "r");
   if (f == NULL)
   {
       fprintf(stderr, "You must place the tris.txt file in the current directory.\n");
       exit(EXIT_FAILURE);
   }
   fseek(f, 0, SEEK_END);
   int numBytes = ftell(f);
   fseek(f, 0, SEEK_SET);
   if (numBytes != 241511792)
   {
       fprintf(stderr, "Your tris.txt file is corrupted.  It should be 241511792 bytes, but you only have %d.\n", numBytes);
       exit(EXIT_FAILURE);
   }

   if (small_read == 1)
   {
       numBytes = 10000;
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
 
   if (numTriangles != 2566541)
   {
       fprintf(stderr, "Issue with reading file -- can't establish number of triangles.\n");
       exit(EXIT_FAILURE);
   }

   if (small_read == 1)
       numTriangles = 100;

   TriangleList *tl = (TriangleList *) malloc(sizeof(TriangleList));
   tl->numTriangles = numTriangles;
   tl->triangles = (Triangle *) malloc(sizeof(Triangle)*tl->numTriangles);

   for (int i = 0 ; i < tl->numTriangles ; i++)
   {
       double x1, y1, x2, y2, x3, y3;
       int    r, g, b;
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
       tmp++,
       x1 = atof(tmp);
       while (*tmp != ',')
          tmp++;
       tmp += 2; // comma+space
       y1 = atof(tmp);
       while (*tmp != ')')
          tmp++;
       tmp += 4; // right-paren+comma+space+left-paren
       x2 = atof(tmp);
       while (*tmp != ',')
          tmp++;
       tmp += 2; // comma+space
       y2 = atof(tmp);
       while (*tmp != ')')
          tmp++;
       tmp += 4; // right-paren+comma+space+left-paren
       x3 = atof(tmp);
       while (*tmp != ',')
          tmp++;
       tmp += 2; // comma+space
       y3 = atof(tmp);
       while (*tmp != ')')
          tmp++;
       tmp += 5; // right-paren+space+equal+space+left-paren
       r = atoi(tmp);
       while (*tmp != ',')
          tmp++;
       tmp += 2; // comma+space
       g = atoi(tmp);
       while (*tmp != ',')
          tmp++;
       tmp += 2; // comma+space
       b = atoi(tmp);
       while (*tmp != '\n')
          tmp++;
       tmp++; // onto next line
       
       tl->triangles[i].X[0] = x1;
       tl->triangles[i].X[1] = x2;
       tl->triangles[i].X[2] = x3;
       tl->triangles[i].Y[0] = y1;
       tl->triangles[i].Y[1] = y2;
       tl->triangles[i].Y[2] = y3;
       tl->triangles[i].color[0] = r;
       tl->triangles[i].color[1] = g;
       tl->triangles[i].color[2] = b;
       //printf("Read triangle %f, %f, %f, %f, %f, %f, %d, %d, %d\n", x1, y1, x2, y2, x3, y3, r, g, b);
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
            if (isnan(leftEnd)) leftEnd = t->Y[left];

            // right side
            if (t->Y[right] == t->Y[topOrBottom]) continue;
            rightSlope  =   (t->Y[right] - t->Y[topOrBottom]) /
                            (t->X[right] - t->X[topOrBottom]);
            bRight      =   (t->Y[topOrBottom] - rightSlope * t->X[topOrBottom]);
            rightEnd    =   (r - bRight) / rightSlope;
            if (isnan(rightEnd)) rightEnd = t->Y[right];

            if (leftEnd > rightEnd) 
            {
                double temp = leftEnd;
                leftEnd = rightEnd;
                rightEnd = temp;
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
            if (leftIntercept > rightIntercept) continue;

            // check scanline length
            if (rightIntercept - leftIntercept > 5) abort();

            // check scanline count
            if (topOrBottomScanline - middleScanline > 5) abort();
            
            for (int c = leftIntercept; c <= rightIntercept; c++)
            {
                ColorPixel(img, r, c, t->color[0], t->color[1], t->color[2]);
            }
        }
    }
}
