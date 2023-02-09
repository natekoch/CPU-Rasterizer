/* proj1F.c
 * Author: Nate Koch
 * Date: 2/6/23
 * LINK TO VIDEO: https://youtube.com/shorts/x2kGc2y-U8w
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NORMALS

#define MOVIE_MODE 1

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
    double          A[4][4];     // A[i][j] means row i, column j
} Matrix;

typedef struct
{
    double  X[3];
    double  Y[3];
    double  Z[3];
    double  colors[3][3];
    double  shadingValue[3];
#ifdef NORMALS
    double  normals[3][3]; // normals[2][0] is for V2, x-component
#endif
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

typedef struct
{
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
} Camera;

typedef struct 
{
    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
} LightingParameters;

double max(double a, double b);

double min(double a, double b);

void PrintMatrix(Matrix m);

Matrix ComposeMatrices(Matrix M1, Matrix M2);

void TransformPoint(Matrix m, const double *ptIn, double *ptOut);

double SineParameterize(int curFrame, int nFrames, int ramp);

Camera GetCamera(int frame, int nframes);

Image* AllocateScreen(int width, int height);

void InitializeScreen(Image *img);

void SaveImage(Image *img, char *filename);

void ColorPixel(Image *img, int row, int col, double R, double G, double B, double Z);

char * Read3Numbers(char *tmp, double *v1, double *v2, double *v3);

TriangleList* Get3DTriangles(void);

void RasterizeTriangle(Triangle *t, Image *img);

void TransformAndRenderTriangles(Camera c, TriangleList *tl, Image *img);

Matrix GetViewTransform(Camera c);

Matrix GetCameraTransform(Camera c);

Matrix GetDeviceTransform(Camera c, double n, double m);

LightingParameters GetLighting(Camera c);

double CalculateShading(LightingParameters *lp, double *viewDirection, double *normal);

int
main(void)
{
    char filename[28];
    
    TriangleList *tl = Get3DTriangles();
    Image *img = AllocateScreen(1000, 1000);
    for (int i = 0; i < 1000; i++)
    {
        if (i != 0) continue;

        InitializeScreen(img);
        Camera c = GetCamera(i, 1000);
        TransformAndRenderTriangles(c, tl, img);
        sprintf(filename, "proj1F_frame%04d.pnm", i);
        //sprintf(filename, "frames/proj1F_frame%04d.pnm", i);
        SaveImage(img, filename);
    }

    free(tl->triangles);
    free(tl);

    free(img->data);
    free(img->zbuffer);
    free(img);

    return 0; 
}

double
max(double a, double b)
{
    return (a > b) ? a : b;
}

double
min(double a, double b)
{
    return (a < b) ? a : b;
}

void
PrintMatrix(Matrix m)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        printf("(%.7f %.7f %.7f %.7f)\n", m.A[i][0], m.A[i][1], m.A[i][2], m.A[i][3]);
    }
}

Matrix
ComposeMatrices(Matrix M1, Matrix M2)
{
    Matrix m_out;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            m_out.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                m_out.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }
    return m_out;
}

void 
TransformPoint(Matrix m, const double *ptIn, double *ptOut)
{  
    ptOut[0] = ptIn[0]*m.A[0][0]
             + ptIn[1]*m.A[1][0]
             + ptIn[2]*m.A[2][0]
             + ptIn[3]*m.A[3][0];
    ptOut[1] = ptIn[0]*m.A[0][1]
             + ptIn[1]*m.A[1][1]
             + ptIn[2]*m.A[2][1]
             + ptIn[3]*m.A[3][1];
    ptOut[2] = ptIn[0]*m.A[0][2]
             + ptIn[1]*m.A[1][2]
             + ptIn[2]*m.A[2][2]
             + ptIn[3]*m.A[3][2];
    ptOut[3] = ptIn[0]*m.A[0][3]
             + ptIn[1]*m.A[1][3]
             + ptIn[2]*m.A[2][3]
             + ptIn[3]*m.A[3][3];
}

double SineParameterize(int curFrame, int nFrames, int ramp)
{  
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {        
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }        
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
} 

Camera       
GetCamera(int frame, int nframes)
{            
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0; 
    c.focus[1] = 0; 
    c.focus[2] = 0;
    c.up[0] = 0;    
    c.up[1] = 1;    
    c.up[2] = 0;    
    return c;       
}

Image* 
AllocateScreen(int width, int height)
{    
    Image *img;
    img = malloc(sizeof(Image));
    img->x = width;
    img->y = height;
    img->data = malloc(sizeof(Pixel) * img->x * img->y);
    img->zbuffer = malloc(sizeof(double) * img->x * img->y);
    return img;
}

void 
InitializeScreen(Image *img)
{
    // initialize img data to all zeros
    memset(img->data, 0, sizeof(Pixel) * img->x * img->y); 
    // initial z buffer to all -1
    for (int i = 0; i < img->x*img->y; i++) img->zbuffer[i] = -1;
}

void 
SaveImage(Image *img, char *filename) 
{
    FILE *pnm = fopen(filename, "wb");

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

void 
ColorPixel(Image *img, int row, int col, double R, double G, double B, double Z)
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
Read3Numbers(char *tmp, double *v1, double *v2, double *v3)
{
    *v1 = atof(tmp);
    while (*tmp != ' ')
       tmp++;
    tmp++; /* space */
    *v2 = atof(tmp);
    while (*tmp != ' ')
       tmp++;
    tmp++; /* space */
    *v3 = atof(tmp);
    while (*tmp != ' ' && *tmp != '\n')
       tmp++;
    return tmp;
}

TriangleList *
Get3DTriangles()
{
   FILE *f = fopen("ws_tris.txt", "r");
   if (f == NULL)
   {
       fprintf(stderr, "You must place the ws_tris.txt file in the current directory.\n");
       exit(EXIT_FAILURE);
   }
   fseek(f, 0, SEEK_END);
   int numBytes = ftell(f);
   fseek(f, 0, SEEK_SET);
   if (numBytes != 3892295)
   {
       fprintf(stderr, "Your ws_tris.txt file is corrupted.  It should be 3892295 bytes, but you have %d.\n", numBytes);
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
 
   if (numTriangles != 14702)
   {
       fprintf(stderr, "Issue with reading file -- can't establish number of triangles.\n");
       exit(EXIT_FAILURE);
   }

   TriangleList *tl = (TriangleList *) malloc(sizeof(TriangleList));
   tl->numTriangles = numTriangles;
   tl->triangles = (Triangle *) malloc(sizeof(Triangle)*tl->numTriangles);

   for (int i = 0 ; i < tl->numTriangles ; i++)
   {
       for (int j = 0 ; j < 3 ; j++)
       {
           double x, y, z;
           double r, g, b;
           double normals[3];
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
           tmp = Read3Numbers(tmp, &x, &y, &z);
           tmp += 3; /* space+slash+space */
           tmp = Read3Numbers(tmp, &r, &g, &b);
           tmp += 3; /* space+slash+space */
           tmp = Read3Numbers(tmp, normals+0, normals+1, normals+2);
           tmp++;    /* newline */

           tl->triangles[i].X[j] = x;
           tl->triangles[i].Y[j] = y;
           tl->triangles[i].Z[j] = z;
           tl->triangles[i].colors[j][0] = r;
           tl->triangles[i].colors[j][1] = g;
           tl->triangles[i].colors[j][2] = b;
#ifdef NORMALS
           tl->triangles[i].normals[j][0] = normals[0];
           tl->triangles[i].normals[j][1] = normals[1];
           tl->triangles[i].normals[j][2] = normals[2];
#endif
       }
   }
   fclose(f);
   free(buffer);
   return tl;
}

void
crossProduct(const double *vecA, const double *vecB, double *vecC)
{
    vecC[0] = vecA[1] * vecB[2] - vecA[2] * vecB[1];
    vecC[1] = vecA[2] * vecB[0] - vecA[0] * vecB[2];
    vecC[2] = vecA[0] * vecB[1] - vecA[1] * vecB[0];
}

double
dotProduct(const double *vecA, const double *vecB)
{
    double product = 0.0;
    for (int i = 0; i < 3; i++)
        product += vecA[i] * vecB[i];
    return product;
}

void
normalizeVector(const double *vecIn, double *vecOut)
{   
    double vecLength = sqrt((vecIn[0] * vecIn[0]) +
                            (vecIn[1] * vecIn[1]) +
                            (vecIn[2] * vecIn[2]));

    vecOut[0] = vecIn[0]/vecLength;
    vecOut[1] = vecIn[1]/vecLength;
    vecOut[2] = vecIn[2]/vecLength;
}

void
RasterizeTriangle(Triangle *t, Image *img)
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

            double shadingValueLeft, shadingValueRight, shadingValuePixel;

            // interpolate color at left and right ends 
            colorLeftEnd[0] = t->colors[left][0] + tPropLeft * (t->colors[topOrBottom][0] - t->colors[left][0]);
            colorLeftEnd[1] = t->colors[left][1] + tPropLeft * (t->colors[topOrBottom][1] - t->colors[left][1]);
            colorLeftEnd[2] = t->colors[left][2] + tPropLeft * (t->colors[topOrBottom][2] - t->colors[left][2]);

            colorRightEnd[0] = t->colors[right][0] + tPropRight * (t->colors[topOrBottom][0] - t->colors[right][0]);  
            colorRightEnd[1] = t->colors[right][1] + tPropRight * (t->colors[topOrBottom][1] - t->colors[right][1]);
            colorRightEnd[2] = t->colors[right][2] + tPropRight * (t->colors[topOrBottom][2] - t->colors[right][2]);

            shadingValueLeft = t->shadingValue[left] + tPropLeft * (t->shadingValue[topOrBottom] - t->shadingValue[left]);

            shadingValueRight = t->shadingValue[right] + tPropRight * (t->shadingValue[topOrBottom] - t->shadingValue[right]);
           
            for (int c = leftIntercept; c <= rightIntercept; c++)
            { 
                tPropPixel = (c - leftEnd) / (rightEnd - leftEnd);
                if (rightEnd != leftEnd)
                {
                    zPixel = zLeftEnd + tPropPixel * (zRightEnd - zLeftEnd);
                    colorPixel[0] = colorLeftEnd[0] + tPropPixel * (colorRightEnd[0] - colorLeftEnd[0]);
                    colorPixel[1] = colorLeftEnd[1] + tPropPixel * (colorRightEnd[1] - colorLeftEnd[1]);
                    colorPixel[2] = colorLeftEnd[2] + tPropPixel * (colorRightEnd[2] - colorLeftEnd[2]);
                    shadingValuePixel = shadingValueLeft + tPropPixel * (shadingValueRight - shadingValueLeft);
                }
                else
                {
                    zPixel = zRightEnd;
                    colorPixel[0] = colorRightEnd[0];
                    colorPixel[1] = colorRightEnd[1];
                    colorPixel[2] = colorRightEnd[2];
                    shadingValuePixel = shadingValueRight;
                }
                //printf("SVP: %f\n", shadingValuePixel);
                colorPixel[0] = min(1.0, colorPixel[0]*shadingValuePixel);
                colorPixel[1] = min(1.0, colorPixel[1]*shadingValuePixel);
                colorPixel[2] = min(1.0, colorPixel[2]*shadingValuePixel);

                ColorPixel(img, r, c, colorPixel[0], colorPixel[1], colorPixel[2], zPixel);
            }
        }
    }
}

void 
TransformAndRenderTriangles(Camera c, TriangleList *tl, Image *img)
{
    Matrix cameraTransform = GetCameraTransform(c);
    Matrix viewTransform = GetViewTransform(c);
    Matrix deviceTransform = GetDeviceTransform(c, img->x, img->y);

    Matrix M = ComposeMatrices(ComposeMatrices(cameraTransform, viewTransform), deviceTransform);
    
    LightingParameters lp = GetLighting(c);

    Triangle *newT = malloc(sizeof(Triangle));
    double pointIn[4];
    double newV[3][4];

    double viewDirection[3]; 
    double viewDirectionNormal[3];

    double shadingValue;

    for (int i = 0; i < tl->numTriangles; i++)
    {
        for (int vertexId = 0; vertexId < 3; vertexId++)
        {
            pointIn[0] = (tl->triangles+i)->X[vertexId];
            pointIn[1] = (tl->triangles+i)->Y[vertexId];
            pointIn[2] = (tl->triangles+i)->Z[vertexId];
            pointIn[3] = 1;
            TransformPoint(M, pointIn, newV[vertexId]);

            // divide newV[vertexId]'s x, y, & z by its w
            newT->X[vertexId] = newV[vertexId][0] / newV[vertexId][3];
            newT->Y[vertexId] = newV[vertexId][1] / newV[vertexId][3];
            newT->Z[vertexId] = newV[vertexId][2] / newV[vertexId][3];

            viewDirection[0] = c.position[0]-(tl->triangles+i)->X[vertexId];
            viewDirection[1] = c.position[1]-(tl->triangles+i)->Y[vertexId];
            viewDirection[2] = c.position[2]-(tl->triangles+i)->Z[vertexId];

            normalizeVector(viewDirection, viewDirectionNormal);

            (tl->triangles+i)->shadingValue[vertexId] = CalculateShading(&lp, viewDirectionNormal, (tl->triangles+i)->normals[vertexId]);
        
        }
        
        memcpy(newT->colors, (tl->triangles+i)->colors, sizeof(double)*3*3);
        memcpy(newT->shadingValue, (tl->triangles+i)->shadingValue, sizeof(double)*3);
        memcpy(newT->normals, (tl->triangles+i)->normals, sizeof(double)*3*3);

        RasterizeTriangle(newT, img);
    }
    free(newT);
}

double
cot(double a)
{
    return cos(a)/sin(a);
}

Matrix 
GetViewTransform(Camera c)
{
    Matrix m;

    double A[4][4] = 
    {
                    {cot(c.angle/2), 0, 0, 0},
                    {0, cot(c.angle/2), 0, 0},
                    {0, 0, (c.far+c.near)/(c.far-c.near), -1},
                    {0, 0, (2*c.far*c.near)/(c.far-c.near), 0}
    };

    memcpy(m.A, A, sizeof(double)*4*4);

    return m;
}

Matrix
GetCameraTransform(Camera c)
{   
    Matrix rv;
    
    /* Construct camera frame */
    double O[3];
    memcpy(O, c.position, sizeof(double) * 3);
    
    double w[3];
    w[0] = O[0] - c.focus[0];
    w[1] = O[1] - c.focus[1];
    w[2] = O[2] - c.focus[2];
    
    double u[3];
    crossProduct(c.up, w, u);

    double v[3];
    crossProduct(w, u, v);

    // Normalize vectors
    double uN[3], vN[3], wN[3], ON[3];
    normalizeVector(u, uN);
    normalizeVector(v, vN);
    normalizeVector(w, wN);
    normalizeVector(O, ON);

    /* Construct camera transform */
    double t[3];
    t[0] = 0 - ON[0];
    t[1] = 0 - ON[1];
    t[2] = 0 - ON[2];

    double A[4][4] = 
    { 
                    {uN[0], vN[0], wN[0], 0},
                    {uN[1], vN[1], wN[1], 0},
                    {uN[2], vN[2], wN[2], 0},
                    {dotProduct(u, t), dotProduct(v, t), dotProduct(w, t), 1}
    };

    memcpy(rv.A, A, sizeof(double)*4*4);
    
    return rv;
}

Matrix
GetDeviceTransform(Camera c, double n, double m)
{   
    Matrix rv;

    double A[4][4] = 
    { 
                    {n/2, 0, 0, 0},
                    {0, m/2, 0, 0},
                    {0, 0, 1, 0},
                    {n/2, m/2, 0, 1}
    };

    memcpy(rv.A, A, sizeof(double)*4*4);

    return rv;
}

LightingParameters 
GetLighting(Camera c)
{
    LightingParameters lp;
    lp.Ka = 0.3;
    lp.Kd = 0.7;
    lp.Ks = 2.8;
    lp.alpha = 50.5;
    lp.lightDir[0] = c.position[0]-c.focus[0];
    lp.lightDir[1] = c.position[1]-c.focus[1];
    lp.lightDir[2] = c.position[2]-c.focus[2];
    double mag = sqrt(lp.lightDir[0]*lp.lightDir[0]
                    + lp.lightDir[1]*lp.lightDir[1]
                    + lp.lightDir[2]*lp.lightDir[2]);
    if (mag > 0)
    {
        lp.lightDir[0] /= mag;
        lp.lightDir[1] /= mag;
        lp.lightDir[2] /= mag;
    }

    return lp;
}

double 
CalculateShading(LightingParameters *lp, double *viewDirection, double *normal)
{
    double LdotN, diffuse, RdotV, specular, shadingAmount;
    double reflection[3], reflectionNormal[3];

    LdotN = dotProduct(lp->lightDir, normal); 
    diffuse = lp->Kd * max(0.0, LdotN);

    reflection[0] = 2 * LdotN * normal[0] - lp->lightDir[0];
    reflection[1] = 2 * LdotN * normal[1] - lp->lightDir[1];
    reflection[2] = 2 * LdotN * normal[2] - lp->lightDir[2]; 

    normalizeVector(reflection, reflectionNormal);
    RdotV = dotProduct(reflectionNormal, viewDirection);
    specular = fabs(lp->Ks * (pow(max(0.0, RdotV), lp->alpha)));
    
    shadingAmount = lp->Ka + diffuse + specular;
    
    return shadingAmount;
}

