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

//void AssignPixels(Image *img, int colMin, int colMax, int rowMin, int rowMax, 
                   // unsigned char R, unsigned char G, unsigned char B);

int main(void)
{
    
    Image *img = NULL;

    img = malloc(sizeof(Image));
    img->x = 300;
    img->y = 300;
    img->data = malloc(sizeof(Pixel) * img->x * img->y);
    
    int index = 0;
    for (int row = 0; row < 300; row++)
    {
        for (int col = 0; col < 300; col++) 
        {
            img->data[index].R = 255;
            img->data[index].G = 255;
            img->data[index].B = 255;
            index++;
        }
    }

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