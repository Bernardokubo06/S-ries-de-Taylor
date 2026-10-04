"""
Entrega 1: Série de Taylor. Função: y = x * e^x                                
""" 

import csv
import math
import os
import timeit

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

"""Função:y = x * e^x"""
def f(x):
    return x * math.exp(x)

"""Calculo Exponencial"""
def exp_taylor(x, n_termos):
    return sum(x**n / math.factorial(n) for n in range(n_termos))

"""Série de Taylor"""
def taylor(x, n_termos):
    return x * exp_taylor(x, n_termos)


# Cores Grafico
SUPERFICIE = "#fcfcfb"
TEXTO = "#0b0b0b"
TEXTO_2 = "#52514e"
GRADE = "#e3e2de"
AZUL, LARANJA, AQUA = "#2a78d6", "#eb6834", "#1baf7a"

plt.rcParams.update({
    "figure.facecolor": SUPERFICIE,
    "axes.facecolor": SUPERFICIE,
    "savefig.facecolor": SUPERFICIE,
    "text.color": TEXTO,
    "axes.labelcolor": TEXTO_2,
    "xtick.color": TEXTO_2,
    "ytick.color": TEXTO_2,
    "axes.edgecolor": GRADE,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 10,
})

#Aproximação
def polinomio_str(n):
    """Polinômio de Taylor com n termos, em texto."""
    termos = []
    for k in range(1, n + 1):
        fat = math.factorial(k - 1)
        if k == 1:
            termos.append("x")
        elif fat == 1:
            termos.append(f"x^{k}")
        else:
            termos.append(f"x^{k}/{fat}")
    return " + ".join(termos)


def limite_lagrange(x, n):
    """Limite do erro (resto de Lagrange) para x > 0 com n termos (grau n).

    f^(m)(x) = (x + m) * e^x  ->  |R_n(x)| <= (x + n + 1) * e^x * x^(n+1) / (n+1)!
    """
    return (x + n + 1) * math.exp(x) * x ** (n + 1) / math.factorial(n + 1)


def propriedades():
    print("=" * 64)
    print("APROXIMAÇÕES (polinômios de Taylor de x*e^x em torno de 0)")
    print("=" * 64)
    for n in (1, 2, 3, 4, 5):
        print(f"P{n}(x) = {polinomio_str(n)}")

    print()
    print("=" * 64)
    print("PROPRIEDADES")
    print("=" * 64)
    print("- Série: x*e^x = soma_{k>=1} x^k / (k-1)!  (coeficientes 1/(k-1)!)")
    print("- Raio de convergência: infinito (converge para todo x real)")
    print("- Derivada n-ésima: f^(n)(x) = (x + n) * e^x   ->   f^(n)(0) = n")
    print("- f(0) = 0 ;  f'(x) = (1 + x) e^x  ->  ponto crítico em x = -1")
    print(f"- Mínimo global em x = -1: f(-1) = -1/e = {-1/math.e:.6f}")
    print("- Limite: f(x) -> 0 quando x -> -infinito; f(x) -> +infinito quando x -> +infinito")
    print("- Primitiva: integral de x*e^x dx = (x - 1) e^x + C")
    print("- Erro de truncamento (grau n): R_n(x) = f^(n+1)(c) x^(n+1) / (n+1)!, c entre 0 e x")
    print("- Perto de 0 o erro é aproximadamente o 1º termo descartado: x^(n+1) / n!")

    print()
    print("Erro real vs. limite de Lagrange em x = 1:")
    print(f"{'termos':>7} | {'erro real':>11} | {'limite Lagrange':>15} | {'1º termo descartado':>19}")
    print("-" * 62)
    x = 1.0
    for n in (2, 4, 6, 8, 10):
        real = abs(f(x) - taylor(x, n))
        lim = limite_lagrange(x, n)
        prox = x ** (n + 1) / math.factorial(n)
        ok = "OK" if real <= lim else "ERRO"
        print(f"{n:>7} | {real:>11.2e} | {lim:>15.2e} | {prox:>19.2e}  {ok}")
    print()



#Tabela de valores
VALORES_X = [-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3]
ORDENS = (3, 5, 10)


def tabela_valores():
    cab = ["x", "exato x*e^x"] + [f"Taylor {n} termos" for n in ORDENS] + [
        f"erro {n} termos" for n in ORDENS
    ]
    linhas = []
    for x in VALORES_X:
        exato = f(x)
        aprox = [taylor(x, n) for n in ORDENS]
        erros = [abs(exato - a) for a in aprox]
        linhas.append([x, exato] + aprox + erros)

    print("=" * 100)
    print("TABELA DE VALORES MAIS USADOS")
    print("=" * 100)
    print(" | ".join(f"{c:>15}" for c in cab))
    print("-" * 100)
    for lin in linhas:
        print(" | ".join(
            [f"{lin[0]:>15g}"]
            + [f"{v:>15.8f}" for v in lin[1:2 + len(ORDENS)]]
            + [f"{v:>15.2e}" for v in lin[2 + len(ORDENS):]]
        ))
    print()

    with open("tabela_valores.csv", "w", newline="", encoding="utf-8") as arq:
        w = csv.writer(arq)
        w.writerow(cab)
        w.writerows(linhas)


#Gráficos
def grafico_aproximacoes():
    """Função exata x aproximações de Taylor (serie_taylor.png)."""
    import numpy as np

    xs = np.linspace(-3, 2, 400)
    fig, ax = plt.subplots(figsize=(8, 5), dpi=120)
    ax.plot(xs, xs * np.exp(xs), color=TEXTO, linewidth=2.5, label="x·eˣ (exata)")
    for n, cor in ((2, AZUL), (3, LARANJA), (5, AQUA)):
        ax.plot(xs, [taylor(v, n) for v in xs], "--", color=cor, linewidth=2,
                label=f"Taylor ({n} termos)")
    ax.set_ylim(-3, 15)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Série de Taylor de y = x·eˣ", loc="left", color=TEXTO, fontsize=12)
    ax.grid(True, color=GRADE, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig("serie_taylor.png")
    plt.close(fig)


def grafico_erros():
    ns = list(range(1, 16))
    series = [(-2.0, AZUL), (1.0, LARANJA), (3.0, AQUA)]

    fig, ax = plt.subplots(figsize=(8, 5), dpi=120)
    for x, cor in series:
        erros = [max(abs(f(x) - taylor(x, n)), 1e-17) for n in ns]
        ax.plot(ns, erros, color=cor, linewidth=2, marker="o", markersize=6,
                markeredgecolor=SUPERFICIE, markeredgewidth=1.5, label=f"x = {x:g}")
        ax.annotate(f"x = {x:g}", (ns[-1], erros[-1]), xytext=(8, 0),
                    textcoords="offset points", va="center", color=TEXTO_2, fontsize=9)

    ax.set_yscale("log")
    ax.set_xlim(0.5, ns[-1] + 2.2)
    ax.set_xticks(ns)
    ax.set_xlabel("Número de termos da série")
    ax.set_ylabel("Erro absoluto |exato − Taylor| (escala log)")
    ax.set_title("Erro da série de Taylor de x·eˣ", loc="left", color=TEXTO, fontsize=12)
    ax.grid(True, axis="y", color=GRADE, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig("grafico_erros.png")
    plt.close(fig)


def medir_tempos(ns, x=1.0, repeticoes=5, chamadas=2000):
    """Tempo médio por chamada (em microssegundos): melhor de 'repeticoes' rodadas."""
    t_taylor = []
    for n in ns:
        t = min(timeit.repeat(lambda: taylor(x, n), number=chamadas, repeat=repeticoes))
        t_taylor.append(t / chamadas * 1e6)
    t_exata = min(timeit.repeat(lambda: f(x), number=chamadas * 10, repeat=repeticoes))
    t_exata = t_exata / (chamadas * 10) * 1e6
    return t_taylor, t_exata


def grafico_tempo():
    ns = list(range(1, 31))
    t_taylor, t_exata = medir_tempos(ns)

    print("=" * 64)
    print("VELOCIDADE DE PROCESSAMENTO (x = 1, tempo por chamada)")
    print("=" * 64)
    print(f"x * math.exp(x) (biblioteca): {t_exata:.3f} µs")
    for n in (3, 5, 10, 20, 30):
        t = t_taylor[n - 1]
        print(f"Taylor com {n:>2} termos      : {t:.3f} µs  ({t / t_exata:.1f}x o tempo da exata)")
    print()

    fig, ax = plt.subplots(figsize=(8, 5), dpi=120)
    ax.plot(ns, t_taylor, color=AZUL, linewidth=2, marker="o", markersize=6,
            markeredgecolor=SUPERFICIE, markeredgewidth=1.5, label="Série de Taylor (x · série de eˣ)")
    ax.axhline(t_exata, color=LARANJA, linewidth=2, label="x · math.exp(x) (biblioteca)")
    ax.annotate("Taylor", (ns[-1], t_taylor[-1]), xytext=(8, 0),
                textcoords="offset points", va="center", color=TEXTO_2, fontsize=9)
    ax.annotate("math.exp", (ns[-1], t_exata), xytext=(8, 0),
                textcoords="offset points", va="center", color=TEXTO_2, fontsize=9)

    ax.set_xlim(0, ns[-1] + 3.5)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Número de termos da série")
    ax.set_ylabel("Tempo por chamada (µs)")
    ax.set_title("Velocidade de processamento", loc="left", color=TEXTO, fontsize=12)
    ax.grid(True, axis="y", color=GRADE, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig("grafico_tempo.png")
    plt.close(fig)


if __name__ == "__main__":
    propriedades()
    tabela_valores()
    grafico_aproximacoes()
    grafico_erros()
    grafico_tempo()
    print("Arquivos gerados: tabela_valores.csv, serie_taylor.png, "
          "grafico_erros.png, grafico_tempo.png")
    print(f"Pasta onde foram salvos: {os.getcwd()}")

 
    #input("\nPressione Enter para sair...")

