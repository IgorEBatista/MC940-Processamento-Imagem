#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ift.h"

typedef int (*WeightFunc)(iftMImage *img, iftImage *root, iftImage *label, int **tree_sum, int *tree_size, int *tree_label, int **obj_sum, int *obj_size, int p, int q, int n_seeds);

int w1(iftMImage *img, iftImage *root, iftImage *label, int **tree_sum, int *tree_size, int *tree_label, int **obj_sum, int *obj_size, int p, int q, int n_seeds){
    double w1 = 0;
    for (int i = 0; i < img->m; i++)
    {
        w1 += pow((tree_sum[i][root->val[p]] / tree_size[root->val[p]]) - (int)img->val[q][i], 2);
    }
    return (int)sqrt(w1);
}

int w2(iftMImage *img, iftImage *root, iftImage *label, int **tree_sum, int *tree_size, int *tree_label, int **obj_sum, int *obj_size, int p, int q, int n_seeds){
    int w2 = 0, dist = 0, min_dist = IFT_INFINITY_INT, min_idx = -1;
    for (int j = 0; j < n_seeds; j++)
    {
        if (tree_label[j] == label->val[p]) // só compara com árvores do mesmo objeto
        {
            dist = 0;
            for (int i = 0; i < img->m; i++)
            {
                dist += pow((tree_sum[i][j] / tree_size[j]) - (int)img->val[q][i], 2);
            }
            if (dist < min_dist) {
                min_dist = dist;
                min_idx = j;
            }
        }
    }
    w2 = min_dist;
    return sqrt(w2);
}

int w3(iftMImage *img, iftImage *root, iftImage *label, int **tree_sum, int *tree_size, int *tree_label, int **obj_sum, int *obj_size, int p, int q, int n_seeds){
    int w3 = 0;
    for (int i = 0; i < img->m; i++)
    {
        w3 += pow((obj_sum[i][label->val[p]] / obj_size[label->val[p]]) - (int)img->val[q][i], 2);
    }
    return sqrt(w3);
}

int w4(iftMImage *img, iftImage *root, iftImage *label, int **tree_sum, int *tree_size, int *tree_label, int **obj_sum, int *obj_size, int p, int q, int n_seeds){
    double w4a = 0;
    for (int i = 0; i < img->m; i++)
    {
        w4a += pow((int)img->val[q][i] - (int)img->val[p][i], 2);
    }
    return (int)sqrt(w4a) + w1(img, root, label, tree_sum, tree_size, tree_label, obj_sum, obj_size, p, q, n_seeds);
}

int w5(iftMImage *img, iftImage *root, iftImage *label, int **tree_sum, int *tree_size, int *tree_label, int **obj_sum, int *obj_size, int p, int q, int n_seeds){
    int w5 = 0;
    for (int i = 0; i < img->m; i++)
    {
        w5 += pow((int)img->val[q][i] - (int)img->val[p][i], 2);
    }
    return sqrt(w5) + w2(img, root, label, tree_sum, tree_size, tree_label, obj_sum, obj_size, p, q, n_seeds);
}

int w6(iftMImage *img, iftImage *root, iftImage *label, int **tree_sum, int *tree_size, int *tree_label, int **obj_sum, int *obj_size, int p, int q, int n_seeds){
    int w6 = 0;
    for (int i = 0; i < img->m; i++)
    {
        w6 += pow((int)img->val[q][i] - (int)img->val[p][i], 2);
    }
    return sqrt(w6) + w3(img, root, label, tree_sum, tree_size, tree_label, obj_sum, obj_size, p, q, n_seeds);
}




iftImage *dynamic(iftMImage *img, iftAdjRel *A, iftLabeledSet *S, WeightFunc weight_func)
{
    iftGQueue *Q = NULL;
    iftImage *cost = NULL, *label = NULL, *marker = NULL, *root = NULL, *predecessor = NULL;
    int **tree_sum = NULL, *tree_size = NULL, *tree_label = NULL, **obj_sum = NULL, *obj_size = NULL;
    int n_seeds = 0, n_objects = 0;
    int i, j, p, q, tmp, w;
    iftVoxel u, v;
    iftLabeledSet *seeds = S;

    printf("número de bandas: %lu\n", img->m);

    label = iftCreateImage(img->xsize, img->ysize, img->zsize);
    cost  = iftCreateImage(img->xsize, img->ysize, img->zsize);
    // marker = iftCreateImage(img->xsize, img->ysize, img->zsize);
    root = iftCreateImage(img->xsize, img->ysize, img->zsize);
    predecessor = iftCreateImage(img->xsize, img->ysize, img->zsize);
    // Q     = iftCreateGQueue(iftMaximumValue(img) + 1, img->n, cost->val);
    Q     = iftCreateGQueue(2560, img->n, cost->val);


    // inicializa todos os pixels
    for (int p = 0; p < img->n; p++)
    {
        cost->val[p] = IFT_INFINITY_INT;
        label->val[p] = -1;
        root->val[p] = -1;
        predecessor->val[p] = -1;
        iftInsertGQueue(&Q, p);
    }

    n_seeds = 0;
    n_objects = 0;
    // conta as sementes
    while (seeds != NULL)
    {
        p = seeds->elem;
        n_seeds++;
        if (seeds->label > n_objects - 1)
            n_objects = seeds->label + 1;     // conta os objetos pelo maior rótulo
        seeds = seeds->next;
    }
    printf("Número de sementes: %d e número de objetos: %d\n", n_seeds, n_objects);
    
    // aloca estruturas para guardar soma e tamanho das árvores e objetos
    tree_sum = (int **)calloc(img->m, sizeof(int *));
    obj_sum = (int **)calloc(img->m, sizeof(int *));
    tree_label = (int *)calloc(n_seeds, sizeof(int));
    // Guarda soma para cada banda separadamente
    for (int i = 0; i < img->m; i++)
    {   
        tree_sum[i] = (int *)calloc(n_seeds, sizeof(int));
        obj_sum[i] = (int *)calloc(n_objects, sizeof(int));
    }
    
    tree_size = (int *)calloc(n_seeds, sizeof(int));
    obj_size = (int *)calloc(n_objects, sizeof(int));

    // inicializa soma e tamanho das árvores e objetos com os valores das sementes
    seeds = S;
    j = 0;
    while (seeds != NULL) {
        p = seeds->elem;
        iftRemoveGQueueElem(Q, p);
        cost->val[p]  = 0;
        label->val[p] = seeds->label;
        root->val[p]  = j;
        tree_label[j] = seeds->label;
        // marker->val[p] = seeds->marker;

        iftInsertGQueue(&Q, p);
        j++;
        seeds = seeds->next;
    }

    while (!iftEmptyGQueue(Q))
    {
        p = iftRemoveGQueue(Q);
        u = iftGetVoxelCoord(label, p);
        
        // Atualiza soma e tamanho da árvore e do objeto
        if (root->val[p] != -1 && label->val[p] != -1) {
            for (int i = 0; i < img->m; i++)
            {
                tree_sum[i][root->val[p]] += (int)img->val[p][i];
                obj_sum[i][label->val[p]] += (int)img->val[p][i];
            }
            tree_size[root->val[p]]++;
            obj_size[label->val[p]]++;
        }
        for (int j = 1; j < A->n; j++)
        {
            v = iftGetAdjacentVoxel(A, u, j);

            if (iftValidVoxel(img, v) && Q->L.elem[iftGetVoxelIndex(img, v)].color != IFT_BLACK)
            {
                q = iftGetVoxelIndex(img, v);

                w = weight_func(img, root, label, tree_sum, tree_size, tree_label, obj_sum, obj_size, p, q, n_seeds);

                tmp = iftMax(cost->val[p], w);
                if (tmp < cost->val[q])
                {
                    // Atualizar píxel conquistador
                    label->val[q] = label->val[p];
                    predecessor->val[q] = p;
                    root->val[q] = root->val[p];
                    cost->val[q]  = tmp;
                    if (Q->L.elem[q].color == IFT_GRAY) {
                        iftRemoveGQueueElem(Q, q);
                        iftInsertGQueue(&Q, q);
                    } else if (Q->L.elem[q].color == IFT_WHITE) {
                        iftInsertGQueue(&Q, q);
                        printf("insere %d\n", q);
                    }else{
                        printf("erro: elemento preto na fila\n");
                    }
                }
            }
        }
    }
    return label;

}

int main(int argc, char *argv[])
{
  // if (argc!=4){
  //   printf("Usage: %s <P1> <P2> <P3>\n",argv[0]);
  //   printf("P1: imagem original\n");
  //   printf("P2: arquivo de sementes\n");
  //   printf("P3: imagem de rótulos de saída\n");
  //   exit(0);
  // }

  
  
  
  // iftImage  *aux       = iftReadImageByExt(argv[1]);  
  iftImage  *aux       = iftReadImageByExt("data/bird.png");
  iftMImage *img;

  if (iftIsColorImage(aux)){
    img = iftImageToMImage(aux,LAB_CSPACE);
    // img = iftImageToMImage(aux,LABNorm2_CSPACE);
  } else {
    img = iftImageToMImage(aux,GRAY_CSPACE);
  }




  // iftLabeledSet *S     = iftReadSeeds(aux, argv[2]);
  iftLabeledSet *S     = iftReadSeeds(aux, "data/bird-seeds.txt");
  iftAdjRel *A         = iftCircular(1.0);  
  // iftImage *label      = Watershed(img,A,S);
  // iftImage *label      = w1(img,A,S);
  // iftImage *label      = w2(img,A,S);
  // iftImage *label      = w3(img,A,S);
    iftImage *label      = dynamic(img, A, S, w1);
  // iftWriteImageByExt(label,argv[3]);
  iftWriteImageByExt(label,"output/dynamic-bird-1.png");

  
  iftDestroyImage(&aux);  
  iftDestroyMImage(&img);
  iftDestroyAdjRel(&A);
  iftDestroyLabeledSet(&S);
  iftDestroyImage(&label);
  
  return(0);
}