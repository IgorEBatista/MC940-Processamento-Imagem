#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ift.h"

typedef int (*WeightFunc)(iftMImage *img, iftImage *root, iftImage *label, int **tree_sum, int *tree_size, int *tree_label, int **obj_sum, int *obj_size, int p, int q, int n_seeds);

int w1(iftMImage *img, iftImage *root, iftImage *label, int **tree_sum, int *tree_size, int *tree_label, int **obj_sum, int *obj_size, int p, int q, int n_seeds){
    double w1 = 0, x = 0;
    for (int i = 0; i < img->m; i++)
    {
        x = (tree_sum[i][root->val[p]] / tree_size[root->val[p]]) - (int)img->val[q][i];
        w1 += x * x;
    }
    return (int)sqrt(w1);
}

int w2(iftMImage *img, iftImage *root, iftImage *label, int **tree_sum, int *tree_size, int *tree_label, int **obj_sum, int *obj_size, int p, int q, int n_seeds){
    int w2 = 0, dist = 0, x = 0, min_dist = IFT_INFINITY_INT, min_idx = -1;
    for (int j = 0; j < n_seeds; j++)
    {
        if (tree_label[j] == label->val[p]) // só compara com árvores do mesmo objeto
        {
            dist = 0;
            for (int i = 0; i < img->m; i++)
            {
                x = (tree_sum[i][j] / tree_size[j]) - (int)img->val[q][i];
                dist += x * x;
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
    int w3 = 0, x = 0;
    for (int i = 0; i < img->m; i++)
    {
        x = (obj_sum[i][label->val[p]] / obj_size[label->val[p]]) - (int)img->val[q][i];
        w3 += x * x;
    }
    return sqrt(w3);
}

int w4(iftMImage *img, iftImage *root, iftImage *label, int **tree_sum, int *tree_size, int *tree_label, int **obj_sum, int *obj_size, int p, int q, int n_seeds){
    double w4a = 0, x = 0;
    for (int i = 0; i < img->m; i++)
    {
        x = (int)img->val[q][i] - (int)img->val[p][i];
        w4a += x * x;
    }
    return (int)sqrt(w4a) + w1(img, root, label, tree_sum, tree_size, tree_label, obj_sum, obj_size, p, q, n_seeds);
}

int w5(iftMImage *img, iftImage *root, iftImage *label, int **tree_sum, int *tree_size, int *tree_label, int **obj_sum, int *obj_size, int p, int q, int n_seeds){
    int w5 = 0, x = 0;
    for (int i = 0; i < img->m; i++)
    {
        x = (int)img->val[q][i] - (int)img->val[p][i];
        w5 += x * x;
    }
    return sqrt(w5) + w2(img, root, label, tree_sum, tree_size, tree_label, obj_sum, obj_size, p, q, n_seeds);
}

int w6(iftMImage *img, iftImage *root, iftImage *label, int **tree_sum, int *tree_size, int *tree_label, int **obj_sum, int *obj_size, int p, int q, int n_seeds){
    int w6 = 0, x = 0;
    for (int i = 0; i < img->m; i++)
    {
        x = (int)img->val[q][i] - (int)img->val[p][i];
        w6 += x * x;
    }
    return sqrt(w6) + w3(img, root, label, tree_sum, tree_size, tree_label, obj_sum, obj_size, p, q, n_seeds);
}




iftImage *dynamic(iftMImage *img, iftAdjRel *A, iftLabeledSet *S, WeightFunc weight_func)
{
    // iftGQueue *Q = NULL;
    iftFHeap *Q = NULL;
    iftFImage *cost = NULL;
    iftImage *label = NULL, *marker = NULL, *root = NULL, *predecessor = NULL;
    int **tree_sum = NULL, *tree_size = NULL, *tree_label = NULL, **obj_sum = NULL, *obj_size = NULL;
    int n_seeds = 0, n_objects = 0;
    int i, j, p, q, tmp, w;
    iftVoxel u, v;
    iftLabeledSet *seeds = S;

    printf("número de bandas: %lu\n", img->m);

    label = iftCreateImage(img->xsize, img->ysize, img->zsize);
    cost  = iftCreateFImage(img->xsize, img->ysize, img->zsize);
    // marker = iftCreateImage(img->xsize, img->ysize, img->zsize);
    root = iftCreateImage(img->xsize, img->ysize, img->zsize);
    predecessor = iftCreateImage(img->xsize, img->ysize, img->zsize);
    // Q     = iftCreateGQueue(iftMaximumValue(img) + 1, img->n, cost->val);
    // Q     = iftCreateGQueue(2560, img->n, cost->val);
    Q     = iftCreateFHeap(img->n, cost->val);


    // inicializa todos os pixels
    for (int p = 0; p < img->n; p++)
    {
        cost->val[p] = IFT_INFINITY_INT;
        root->val[p] = -1;
        predecessor->val[p] = -1;
        iftInsertFHeap(Q, p);
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
        cost->val[p]  = 0;
        label->val[p] = seeds->label;
        root->val[p]  = j;
        tree_label[j] = seeds->label;
        // marker->val[p] = seeds->marker;
        iftGoUpFHeap(Q, p);

        for (int i = 0; i < img->m; i++)
        {
            tree_sum[i][j] += (int)img->val[p][i];
            obj_sum[i][seeds->label] += (int)img->val[p][i];
        }
        tree_size[j] += 1;
        obj_size[seeds->label] += 1;
        j++;
        seeds = seeds->next;
    }

    while (!iftEmptyFHeap(Q))
    {
        p = iftRemoveFHeap(Q);
        u = iftGetVoxelCoord(label, p);
        
        // Atualiza soma e tamanho da árvore e do objeto
        for (int i = 0; i < img->m; i++)
        {
            tree_sum[i][root->val[p]] += (int)img->val[p][i];
            obj_sum[i][label->val[p]] += (int)img->val[p][i];
        }
        tree_size[root->val[p]]++;
        obj_size[label->val[p]]++;
        
        for (int j = 1; j < A->n; j++)
        {
            v = iftGetAdjacentVoxel(A, u, j);

            if (iftValidVoxel(img, v))
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
                    iftGoUpFHeap(Q, Q->pos[q]);
                    
                }
            }
        }
    }
    return label;

}

int main(int argc, char *argv[])
{
  if (argc < 4){
    printf("Usage: %s <P1> <P2> <P3> (<P4>)\n",argv[0]);
    printf("P1: imagem original\n");
    printf("P2: arquivo de sementes\n");
    printf("P3: imagem de rótulos de saída\n");
    printf("P4: função de custo w1, w2, w3, w4, w5 ou w6 (opcional)\n");
    exit(0);
  }
  
  iftImage  *aux       = iftReadImageByExt(argv[1]);  
//   iftImage  *aux       = iftReadImageByExt("data/bird.png");
  iftMImage *img;

  if (iftIsColorImage(aux)){
    img = iftImageToMImage(aux,LAB_CSPACE);
    // img = iftImageToMImage(aux,LABNorm2_CSPACE);
  } else {
    img = iftImageToMImage(aux,GRAY_CSPACE);
  }




    iftLabeledSet *S     = iftReadSeeds(aux, argv[2]);
    // iftLabeledSet *S     = iftReadSeeds(aux, "data/bird-seeds.txt");
    iftAdjRel *A         = iftCircular(1.0);  
    
    WeightFunc funcs[] = {w1, w2, w3, w4, w5, w6};
    const char *names[] = {"w1", "w2", "w3", "w4", "w5", "w6"};
    int idx = 0; // default w1

    if (argc == 5) {
        for (int i = 0; i < 6; i++) {
            if (strcmp(argv[4], names[i]) == 0) {
                idx = i;
                break;
            }
        }
    } else {
        idx = 0; // default w1
    }
    printf("Usando função de custo %s\n", names[idx]);
    iftImage *label = dynamic(img, A, S, funcs[idx]);

    iftWriteImageByExt(label,argv[3]);
    // iftWriteImageByExt(label,"output/bird-dynamic-w1.png");

  
  iftDestroyImage(&aux);  
  iftDestroyMImage(&img);
  iftDestroyAdjRel(&A);
  iftDestroyLabeledSet(&S);
  iftDestroyImage(&label);
  
  return(0);
}