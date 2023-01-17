/* proj1A.c
 * Author: Nate Koch
 * Date: 1/16/23
 */

#include <stdio.h>
#include <stdlib.h>

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

int main(void)
{
    // Create the Image struct
    Image *img = NULL;

    img = malloc(sizeof(Image));
    img->x = 300;
    img->y = 300;
    img->data = malloc(sizeof(Pixel) * img->x * img->y);

    // Black Square
    AssignPixels(img, 0, 100, 0, 100, 0, 0, 0);
    // Gray Square
    AssignPixels(img, 100, 200, 0, 100, 128, 128, 128);
    // White Square
    AssignPixels(img, 200, 300, 0, 100, 255, 255, 255); 
    // Red Square
    AssignPixels(img, 0, 100, 100, 200, 255, 0, 0);
    // Green Square
    AssignPixels(img, 100, 200, 100, 200, 0, 255, 0);
    // Blue Square
    AssignPixels(img, 200, 300, 100, 200, 0, 0, 255);
    // Pink Square
    AssignPixels(img, 0, 100, 200, 300, 255, 0, 255);
    // Cyan Square
    AssignPixels(img, 100, 200, 200, 300, 0, 255, 255);
    // Yellow Square
    AssignPixels(img, 200, 300, 200, 300, 255, 255, 0);

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
