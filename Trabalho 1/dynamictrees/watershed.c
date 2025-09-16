#include "ift.h"

iftImage *w1(iftMImage *img, iftAdjRel *A, iftLabeledSet *S)
{
  iftImage *basins = iftMImageBasins(img,A);
  iftImage *cost   = NULL, *label = NULL, *tree = NULL;
  iftGQueue     *Q = NULL;
  int i, p, q, tmp;
  iftVoxel    u, v;
  iftLabeledSet *seeds = S;
  int n_seeds, w1;
  int n_objects = 0;
  int *tree_sum = NULL;
  int *tree_size = NULL;

  // Initialization
  cost  = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
  label = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
  Q     = iftCreateGQueue(iftMaximumValue(basins) + 1, basins->n, cost->val);
  
  tree = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
  
  for (p = 0; p < basins->n; p++)
  {
    cost->val[p] = IFT_INFINITY_INT;
    tree->val[p] = -1; // nenhum índice de árvore ainda
  }
  n_seeds = 0;
  while (seeds != NULL)
  {
    p = seeds->elem;
    label->val[p] = seeds->label;
    cost->val[p]  = 0;
    tree->val[p]  = n_seeds; // índice da árvore
    // tree->val[p]  = p; // índice da árvore
    iftInsertGQueue(&Q, p);
    n_seeds++;
    if (seeds->label > n_objects - 1)
      n_objects = seeds->label + 1;     // garante o maior rótulo usado
    seeds = seeds->next;
  }
  tree_sum = (int *) calloc(n_seeds, sizeof(int));
  tree_size = (int *) calloc(n_seeds, sizeof(int));
  // tree_sum = (int *) calloc(basins->n, sizeof(int));
  // tree_size = (int *) calloc(basins->n, sizeof(int));

  printf("Número de sementes: %d\n", n_seeds);
  printf("Número de objetos: %d\n", n_objects);
  
  // Inicializar somas e tamanhos das árvores
  seeds = S;
  while (seeds != NULL) {
    p = seeds->elem;
    tree_sum[tree->val[p]] += (int)img->val[p][0]; // considerando apenas a primeira banda
    tree_size[tree->val[p]] += 1;
    seeds = seeds->next;
  }

  while (!iftEmptyGQueue(Q))
  {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(basins, p);

    for (i = 1; i < A->n; i++)
    {
      v = iftGetAdjacentVoxel(A, u, i);

      if (iftValidVoxel(basins, v))
      {
        q = iftGetVoxelIndex(basins, v);

        w1 = abs((int)img->val[q][0] - tree_sum[tree->val[p]]/tree_size[tree->val[p]]);

        tmp = iftMax(cost->val[p], w1);
        if (tmp < cost->val[q]) {
          tree_sum[tree->val[p]] += (int)img->val[q][0];
          tree_size[tree->val[p]]++;
          
          tree->val[q] = tree->val[p];
          label->val[q] = label->val[p];
          cost->val[q]  = tmp;
          iftInsertGQueue(&Q, q);
        }
        
      }
    }
  }

  printf("Somas e tamanhos finais e médias das árvores:\n");
  for (i = 0; i < n_seeds; i++) {
    printf("Árvore %d: soma = %d, tamanho = %d, média = %.2f\n", i, tree_sum[i], tree_size[i], (float)tree_sum[i] / tree_size[i]);
  }

  // for (i = 0; i < basins->n; i++) {
  //   if (tree_size[i] > 0) {
  //     printf("Árvore %d: soma = %d, tamanho = %d, média = %.2f\n", i, tree_sum[i], tree_size[i], (float)tree_sum[i] / tree_size[i]);
  //   }
  // }


  printf("Somas e tamanhos finais e médias dos objetos:\n");
  int object_sum[2] = {0,0};
  int object_size[2] = {0,0};
  for (p = 0; p < basins->n; p++) {
    int obj = label->val[p];
    if (obj == 0 || obj == 1) { // considerando rótulos 0 e 1
      object_sum[obj] += basins->val[p];
      object_size[obj] += 1;
    }
  }

  for (i = 0; i < n_objects; i++) {
    printf("Objeto %d: soma = %d, tamanho = %d, média = %.2f\n", i, object_sum[i], object_size[i], (float)object_sum[i] / object_size[i]);
  }

  iftDestroyGQueue(&Q);
  iftDestroyImage(&cost);

  return(label);
}

iftImage *w2(iftMImage *img, iftAdjRel *A, iftLabeledSet *S)
{
  iftImage *basins = iftMImageBasins(img,A);
  iftImage *cost   = NULL, *label = NULL;
  iftGQueue     *Q = NULL;
  int i, p, q, tmp;
  iftVoxel    u, v;
  iftLabeledSet *seeds = S;

  int n_seeds, w2, *tree_size, **tree;


  // Initialization
  cost  = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
  label = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
  Q     = iftCreateGQueue(iftMaximumValue(basins) + 1, basins->n, cost->val);

  for (p = 0; p < basins->n; p++)
  {
    cost->val[p] = IFT_INFINITY_INT;
  }
  
  n_seeds = 0;
  while (seeds != NULL)
  {
    p = seeds->elem;
    label->val[p] = seeds->label;
    cost->val[p]  = 0;
    iftInsertGQueue(&Q, p);
    n_seeds++;
    seeds = seeds->next;
  }
  tree = (int **) calloc(n_seeds + 1, sizeof(int *));
  tree_size = (int *) calloc(basins->n, sizeof(int));
  while (seeds != NULL)
  {
    p = seeds->elem;
    tree[seeds->label] = (int *) calloc(basins->n, sizeof(int));
    tree[seeds->label][tree_size[seeds->label]] = basins->val[p];
    tree_size[seeds->label]++;
    seeds = seeds->next;
  }

  while (!iftEmptyGQueue(Q))
  {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(basins, p);

    for (i = 1; i < A->n; i++)
    {
      v = iftGetAdjacentVoxel(A, u, i);

      if (iftValidVoxel(basins, v))
      {
        q = iftGetVoxelIndex(basins, v);
        
        w2 = IFT_INFINITY_INT;
        for (int j = 0; j < tree_size[label->val[p]]; j++)
        {
          int diff = abs(basins->val[q] - tree[label->val[p]][j]);
          if (diff < w2)
            w2 = diff;
        }
        
        tmp = iftMax(cost->val[p], w2);
        if (tmp < cost->val[q]) {
          
          tree[label->val[p]][tree_size[label->val[p]]] = basins->val[q];
          tree_size[label->val[p]]++;
          label->val[q] = label->val[p];
          cost->val[q]  = tmp;
          iftInsertGQueue(&Q, q);
        }
        
      }
    }
  }
  
  iftDestroyGQueue(&Q);
  iftDestroyImage(&cost);

  return(label);
}

iftImage *w3(iftMImage *img, iftAdjRel *A, iftLabeledSet *S)
{
  iftImage *basins = iftMImageBasins(img, A);
  iftImage *cost   = NULL, *label = NULL;
  iftGQueue *Q = NULL;
  int i, p, q, tmp;
  iftVoxel u, v;
  iftLabeledSet *seeds = S;
  
  // Para médias globais por objeto
  int n_objects = 0; 
  int *object_sum = NULL;
  // float *object_sum = NULL;
  int *object_size = NULL;

  // Inicialização
  cost  = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
  label = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
  Q     = iftCreateGQueue(iftMaximumValue(basins) + 1, basins->n, cost->val);

  for (p = 0; p < basins->n; p++) {
    cost->val[p] = IFT_INFINITY_INT;
  }

  while (seeds != NULL) {
    p = seeds->elem;
    label->val[p] = seeds->label;   // rótulo do objeto
    cost->val[p]  = 0;
    iftInsertGQueue(&Q, p);
    if (seeds->label > n_objects - 1)
      n_objects = seeds->label + 1;     // garante o maior rótulo usado
    seeds = seeds->next;
  }
  printf("Número de objetos: %d\n", n_objects);

  // Alocar médias globais por objeto
  object_sum = (int *) calloc(n_objects, sizeof(int));
  // object_sum = (float *) calloc(n_objects, sizeof(float));
  object_size = (int *) calloc(n_objects, sizeof(int));

  // Inicializar com os valores dos seeds
  seeds = S;
  while (seeds != NULL) {
    p = seeds->elem;
    int obj = seeds->label;
    // object_sum[obj] += basins->val[p];
    object_sum[obj] += (int)img->val[p][0];
    object_size[obj] += 1;
    seeds = seeds->next;
  }
 
  printf("Soma e tamanhos iniciais dos objetos:\n");
  for (i = 0; i < n_objects; i++) {
    printf("Objeto %d: soma = %d, tamanho = %d\n", i, object_sum[i], object_size[i]);
  }

  while (!iftEmptyGQueue(Q)) {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(basins, p);

    for (i = 1; i < A->n; i++) {
      v = iftGetAdjacentVoxel(A, u, i);

      if (iftValidVoxel(basins, v)) {
        q = iftGetVoxelIndex(basins, v);
        
        int obj = label->val[p]; // objeto atual
        int w3 = abs(basins->val[q] - (object_sum[obj] / object_size[obj]));
        // int w3 = abs((int)img->val[q][0] - (object_sum[obj] / object_size[obj]));

        tmp = iftMax(cost->val[p], w3);
        if (tmp < cost->val[q]) {
          // Atualizar soma e tamanho do objeto
          object_sum[obj] += basins->val[q];
          // object_sum[obj] += (int)img->val[q][0];
          object_size[obj]++;

          label->val[q] = obj;
          cost->val[q]  = tmp;
          iftInsertGQueue(&Q, q);
        }
      }
    }
  }

  printf("Somas e tamanhos finais e médias dos objetos:\n");
  for (i = 0; i < n_objects; i++) {
    printf("Objeto %d: soma = %d, tamanho = %d, média = %.2f\n", i, object_sum[i], object_size[i], (float)object_sum[i] / object_size[i]);
  }
  
  iftDestroyGQueue(&Q);
  iftDestroyImage(&cost);

  free(object_sum);
  free(object_size);

  return(label);
}

iftImage *Watershed(iftMImage *img, iftAdjRel *A, iftLabeledSet *S)
{
  iftImage *basins = iftMImageBasins(img,A);
  iftImage *cost   = NULL, *label = NULL;
  iftGQueue     *Q = NULL;
  int i, p, q, tmp;
  iftVoxel    u, v;
  iftLabeledSet *seeds = S;

  // Initialization
  cost  = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
  label = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
  Q     = iftCreateGQueue(iftMaximumValue(basins) + 1, basins->n, cost->val);

  for (p = 0; p < basins->n; p++)
  {
    cost->val[p] = IFT_INFINITY_INT;
  }

  
  while (seeds != NULL)
  {
    p = seeds->elem;
    label->val[p] = seeds->label;
    cost->val[p]  = 0;
    iftInsertGQueue(&Q, p);
    seeds = seeds->next;
  }

  while (!iftEmptyGQueue(Q))
  {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(basins, p);

    for (i = 1; i < A->n; i++)
    {
      v = iftGetAdjacentVoxel(A, u, i);

      if (iftValidVoxel(basins, v))
      {
        q = iftGetVoxelIndex(basins, v);
        if (cost->val[q] > cost->val[p])
        {
          tmp = iftMax(cost->val[p], basins->val[q]);
          if (tmp < cost->val[q]) {
            label->val[q] = label->val[p];
            cost->val[q]  = tmp;
            iftInsertGQueue(&Q, q);
          }
        }
      }
    }
  }

  printf("Somas e tamanhos finais e médias dos objetos:\n");
  int object_sum[2] = {0,0};
  int object_size[2] = {0,0};
  for (p = 0; p < basins->n; p++) {
    int obj = label->val[p];
    if (obj == 0 || obj == 1) { // considerando rótulos 0 e 1
      object_sum[obj] += basins->val[p];
      object_size[obj] += 1;
    }
  }
  int n_objects = 2; // considerando rótulos 0 e 1
  for (i = 0; i < n_objects; i++) {
    printf("Objeto %d: soma = %d, tamanho = %d, média = %.2f\n", i, object_sum[i], object_size[i], (float)object_sum[i] / object_size[i]);
  }
  
  iftDestroyGQueue(&Q);
  iftDestroyImage(&cost);

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
  // iftImage *label      = Watershed(img,A,S);
  // iftImage *label      = w1(img,A,S);
  // iftImage *label      = w2(img,A,S);
  iftImage *label      = w3(img,A,S);

  // iftWriteImageByExt(label,argv[3]);

  
  iftDestroyImage(&aux);  
  iftDestroyMImage(&img);
  iftDestroyAdjRel(&A);
  iftDestroyLabeledSet(&S);
  iftDestroyImage(&label);
  
  return(0);
}

