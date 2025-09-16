#include "ift.h"

iftImage *Watershed_Float(iftMImage *img, iftAdjRel *A, iftLabeledSet *S)
{
  iftFImage *cost   = iftCreateFImage(img->xsize, img->ysize, img->zsize);
  iftImage  *label  = iftCreateImage(img->xsize, img->ysize, img->zsize);
  iftFHeap  *Q      = iftCreateFHeap(img->n,cost->val);
  int i, p, q, tmp;
  iftVoxel    u, v;
  iftLabeledSet *seeds = S;

  // Initialization

  for (p = 0; p < cost->n; p++)
  {
    cost->val[p] = IFT_INFINITY_FLT;
  }

  while (seeds != NULL)
  {
    p = seeds->elem;
    label->val[p] = seeds->label;
    cost->val[p]  = 0;
    iftInsertFHeap(Q, p);
    seeds = seeds->next;
  }

  while (!iftEmptyFHeap(Q))
  {
    p = iftRemoveFHeap(Q);
    u = iftMGetVoxelCoord(img, p);

    for (i = 1; i < A->n; i++)
    {
      v = iftGetAdjacentVoxel(A, u, i);

      if (iftMValidVoxel(img, v))
      {
        q = iftMGetVoxelIndex(img, v);
        if (cost->val[q] > cost->val[p])
        {
	  float arcw = 0.0;
	  for (int b=0; b < img->m; b++)
	    arcw += powf(img->val[q][b]-img->val[p][b],2);

          tmp = iftMax(cost->val[p], arcw);
          if (tmp < cost->val[q]) {
            label->val[q] = label->val[p];
            cost->val[q]  = tmp;
	    if (Q->color[q] == IFT_WHITE)	      
	      iftInsertFHeap(Q, q);
	    else
	      iftGoUpFHeap(Q,Q->pos[q]);
          }
        }
      }
    }
  }
  
  iftDestroyFHeap(&Q);
  iftDestroyFImage(&cost);

  return(label);
}

  
int main(int argc, char *argv[])
{
  if (argc!=4){
    printf("Usage: %s <P1> <P2> <P3>\n",argv[0]);
    printf("P1: imagem original\n");
    printf("P2: arquivo de sementes\n");
    printf("P3: imagem de rótulos de saída\n");
    exit(0);
  }

  iftImage  *aux       = iftReadImageByExt(argv[1]);  
  iftMImage *img;

  if (iftIsColorImage(aux)){
    img = iftImageToMImage(aux,LAB_CSPACE);
  } else {
    img = iftImageToMImage(aux,GRAY_CSPACE);
  }
  iftLabeledSet *S     = iftReadSeeds(aux, argv[2]);
  iftAdjRel *A         = iftCircular(1.0);  
  iftImage *label      = Watershed_Float(img,A,S);

  iftImage *slabel     = iftSmoothRegionsByDiffusion(label,aux, 0.5, 2);

  
  iftWriteImageByExt(slabel,argv[3]);
  
  iftDestroyImage(&aux);  
  iftDestroyImage(&slabel);  
  iftDestroyMImage(&img);
  iftDestroyAdjRel(&A);
  iftDestroyLabeledSet(&S);
  iftDestroyImage(&label);
  
  return(0);
}

