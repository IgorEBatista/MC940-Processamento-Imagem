import numpy as np
import subprocess
import os
import time
from skimage.segmentation import find_boundaries
from skimage.io import imread


def boundary_recall(seg, gt, tolerance=1):
    """
    Boundary Recall (BR).
    Compara fronteiras de seg com fronteiras de gt.
    tolerance = raio em pixels para considerar uma borda 'correta'.
    """
    seg_bound = find_boundaries(seg, mode='outer')
    gt_bound = find_boundaries(gt, mode='outer')

    from scipy.ndimage import binary_dilation
    gt_dilated = binary_dilation(gt_bound, iterations=tolerance)

    hits = np.logical_and(seg_bound, gt_dilated).sum()
    total = gt_bound.sum()
    return hits / total if total > 0 else 0.0


def undersegmentation_error(seg, gt):
    """
    Undersegmentation Error (UE).
    Mede quanto os superpixels invadem diferentes regiões do GT.
    """
    ue = 0
    for s in np.unique(seg):
        mask_s = (seg == s)
        overlap = [np.sum(np.logical_and(mask_s, gt == g)) for g in np.unique(gt)]
        ue += (mask_s.sum() - max(overlap))
    return ue / seg.size


def achievable_segmentation_accuracy(seg, gt):
    """
    Achievable Segmentation Accuracy (ASA).
    Mede fração de pixels rotuláveis corretamente se cada superpixel
    fosse rotulado com a classe majoritária.
    """
    correct = 0
    for s in np.unique(seg):
        mask_s = (seg == s)
        overlap = [np.sum(np.logical_and(mask_s, gt == g)) for g in np.unique(gt)]
        correct += max(overlap)
    return correct / seg.size


def dice_coefficient(mask1, mask2):
    """
    Calcula o DICE entre dois conjuntos binários.
    mask1, mask2: arrays booleanos (True = pixel pertencente ao conjunto)
    """
    intersection = np.logical_and(mask1, mask2).sum()
    size1 = mask1.sum()
    size2 = mask2.sum()
    return 2.0 * intersection / (size1 + size2) if (size1 + size2) > 0 else 0.0

def mean_dice(seg, gt):
    """
    Calcula o DICE médio entre superpixels de seg e regiões de gt.
    """
    dice_scores = []
    for s in np.unique(seg):
        mask_s = (seg == s)
        
        # melhor sobreposição no GT
        best = 0
        for g in np.unique(gt):
            mask_g = (gt == g)
            d = dice_coefficient(mask_s, mask_g)
            if d > best:
                best = d
        dice_scores.append(best)
    return np.mean(dice_scores)

def dataset_dice(seg_list, gt_list):
    """
    Calcula média e desvio-padrão do DICE em um conjunto de imagens.
    seg_list: lista de segmentações (arrays 2D)
    gt_list:  lista de referências (arrays 2D)
    """
    scores = []
    for seg, gt in zip(seg_list, gt_list):
        scores.append(mean_dice(seg, gt))
    return np.mean(scores), np.std(scores)

os.makedirs('output', exist_ok=True)

# Encontra o nome dos arquivos
nomes_arquivos = os.listdir('data/')
nomes_arquivos = [f for f in nomes_arquivos if f.endswith('.png') and not '-' in f]
nomes_arquivos = [f[:-4] for f in nomes_arquivos]
nomes_arquivos.sort()

# Processa cada imagem
metricas = []

for nome in nomes_arquivos:
    print(f'Processando {nome}...')

    # Executa o programa C para cada peso diferente
    tempos = []
    for i in range(1, 7):
        comando = f'./dynamic data/{nome}.png data/{nome}-seeds.txt output/{nome}_w{i}.png w{i}'
        tempo_inicial = time.time()
        resultado = subprocess.run(comando, shell=True, capture_output=True, text=True)
        tempo_final = time.time()
        tempos.append(tempo_final - tempo_inicial)
        if resultado.returncode != 0:
            print(f'Erro ao processar {nome}: {resultado.stderr}')
            continue
    # Lê a imagem de saída e o ground truth
    ground_truth = imread(f'data/{nome}-label.png')

    for i in range(1, 7):
        seg = imread(f'output/{nome}_w{i}.png')

        # Calcula métricas
        br = boundary_recall(seg, ground_truth)
        ue = undersegmentation_error(seg, ground_truth)
        asa = achievable_segmentation_accuracy(seg, ground_truth)

        metricas.append({
            'imagem': nome,
            'peso': f"w{i}",
            'tempo': tempos[i-1] if i-1 < len(tempos) else None,
            'BR': br,
            'UE': ue,
            'ASA': asa
        })
        print(f'  Peso w{i}: Tempo={tempos[i-1]:.2f}s, BR={br:.4f}, UE={ue:.4f}, ASA={asa:.4f}')
    print(f"Finalizado {nome}.\n")

# Salva as métricas em um formato para Latex, com cada imagem tendo uma tabela separada
with open('output/metricas.txt', 'w') as f:
    for nome in nomes_arquivos:
        f.write(f'\\begin{{table}}[h]\n\\centering\n')
        f.write(f'\\begin{{tabular}}{{|c|c|c|c|c|}}\n\\hline\n')
        f.write('Peso & Tempo (s) & BR & UE & ASA \\\\\n\\hline\n')
        for metrica in metricas:
            if metrica['imagem'] == nome:
                f.write(f"{metrica['peso']} & {metrica['tempo']:.2f} & {metrica['BR']:.4f} & {metrica['UE']:.4f} & {metrica['ASA']:.4f} \\\\\n")
        f.write('\\hline\n\\end{tabular}\n')
        f.write(f'\\caption{{Métricas para a imagem {nome}}}\n')
        f.write(f'\\label{{tab:{nome}}}\n')
        f.write('\\end{table}\n\n')
