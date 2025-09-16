import sys
"""
visualizador.py

Este script permite carregar e visualizar imagens utilizando OpenCV e Matplotlib, corrigindo automaticamente o formato de cor e o intervalo de valores dos pixels.

Funções:
-----------
load_and_correct_image(path):
    Carrega uma imagem do caminho especificado, converte de BGR para RGB se necessário, e garante que os valores dos pixels estejam no intervalo [0, 255] como uint8.

Como usar:
-----------
Execute o script diretamente pelo terminal, passando o caminho da imagem como argumento:
    python visualizador.py <imagem.png>

Exemplo:
    python visualizador.py minha_imagem.png

Isso abrirá uma janela exibindo a imagem fornecida.
"""
import numpy as np
import cv2
import matplotlib.pyplot as plt



def load_and_correct_image(path):
    img = cv2.imread(path, cv2.IMREAD_UNCHANGED)
    if img is None:
        raise FileNotFoundError(f"Imagem não encontrada: {path}")

    # Converte BGR para RGB se for colorida
    if len(img.shape) == 3 and img.shape[2] == 3:
        img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

    arr = img.astype(np.float32)

    # Se a imagem está em 0 e 1, converte para 0 e 255
    # if arr.max() <= 1.0:
    #     arr = (arr * 255).astype(np.uint8)
    # else:
    #     arr = arr.astype(np.uint8)

    return arr

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Uso: python {sys.argv[0]} <imagem.png>")
        sys.exit(1)

    img_path = sys.argv[1]
    img_arr = load_and_correct_image(img_path)

    plt.imshow(img_arr, cmap='gray' if len(img_arr.shape) == 2 else None)
    plt.axis('off')
    plt.show()

