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

# Pasta de saída (gráficos e CSV)
PASTA = "resultados"
os.makedirs(PASTA, exist_ok=True)

# Série de Taylor (a = 0) de y = x * e^x:  x * e^x = x * soma(x^n / n!)
def f(x):
    """Função: y = x * e^x"""
    return x * math.exp(x)


def exp_taylor(x, n_termos):
    """Cálculo exponencial: soma de x^n / n!"""
    return sum(x**n / math.factorial(n) for n in range(n_termos))


def taylor(x, n_termos):
    """Série de Taylor de x * e^x"""
    return x * exp_taylor(x, n_termos)


# Derivadas feitas à mão (regra do produto), sem Sympy:
#   f' = (x+1)e^x,  f'' = (x+2)e^x,  f''' = (x+3)e^x   ->   f^(n)(x) = (x+n)e^x
#   Em x = 0: f^(n)(0) = n, então o coeficiente é n/n! = 1/(n-1)!
def derivada(n, x):
    """n-ésima derivada de x*e^x"""
    return (x + n) * math.exp(x)


def taylor_por_derivadas(x, n_max):
    """Taylor pela definição: soma de f^(n)(0)/n! * x^n"""
    return sum(derivada(n, 0) / math.factorial(n) * x**n for n in range(n_max + 1))


def mostrar_derivadas():
    print("=" * 64)
    print("DERIVADAS (feitas à mão, sem Sympy)")
    print("=" * 64)
    print("f(x)    = x * e^x")
    print("f'(x)   = e^x + x*e^x       = (x + 1) e^x   (regra do produto)")
    print("f''(x)  = e^x + (x+1)*e^x   = (x + 2) e^x")
    print("f'''(x) = e^x + (x+2)*e^x   = (x + 3) e^x")
    print("Padrão: f^(n)(x) = (x + n) * e^x")
    print()
    print(f"{'n':>3} | {'f^(n)(0)':>9} | {'f^(n)(0)/n!':>12}")
    print("-" * 31)
    for n in range(0, 8):
        d = derivada(n, 0)
        print(f"{n:>3} | {d:>9.0f} | {d / math.factorial(n):>12.6f}")
    print()
    por_derivadas = taylor_por_derivadas(1.0, 10)
    por_exp = taylor(1.0, 10)
    print("Conferência em x = 1 com 10 termos:")
    print(f"  pela definição (derivadas) : {por_derivadas:.10f}")
    print(f"  pela série (exp * x)       : {por_exp:.10f}")
    print()


# Limite N: menor N em que o 1º termo descartado, x^(N+1)/N!, fica menor que
# TOLERANCIA para todo |x| <= X_MAX
X_MAX = 3.0
TOLERANCIA = 1e-15


def descobrir_n_limite(x_max=X_MAX, tol=TOLERANCIA):
    n = 1
    while x_max ** (n + 1) / math.factorial(n) >= tol:
        n += 1
    return n


N_LIMITE = descobrir_n_limite()


def mostrar_n_limite():
    print("=" * 64)
    print("LIMITE N DA SÉRIE (e por quê)")
    print("=" * 64)
    print(f"N_LIMITE = {N_LIMITE}  (vale para |x| <= {X_MAX:g})")
    print("Por quê:")
    print(f"- Com N = {N_LIMITE}, o 1º termo descartado x^(N+1)/N! é menor que {TOLERANCIA:g}.")
    print("- O float (double) guarda ~15 a 16 dígitos: termos menores que isso")
    print("  não mudam mais o resultado.")
    print("- Depois desse N, cada termo extra só aumenta o tempo de execução.")
    print("- Quanto maior |x|, maior o N necessário (veja a tabela).")
    print()
    print(f"{'x':>6} | {'menor N com erro < 1e-12':>25}")
    print("-" * 34)
    for x in VALORES_X:
        n_min = None
        for n in range(1, N_LIMITE + 1):
            if abs(f(x) - taylor(x, n)) < 1e-12:
                n_min = n
                break
        texto = str(n_min) if n_min else f"> {N_LIMITE}"
        print(f"{x:>6g} | {texto:>25}")
    print()


# Função otimizada: cada termo vem do anterior (t_k = t_(k-1) * x / (k-1)),
# sem fatorial nem potência, e para quando o termo fica menor que a tolerância
def taylor_otimizada(x, n_max=N_LIMITE, tol=TOLERANCIA):
    termo = x
    soma = x
    for k in range(2, n_max + 1):
        termo = termo * x / (k - 1)
        soma += termo
        if abs(termo) < tol:
            break
    return soma


# Cores Gráfico
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


# Aproximações e propriedades
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
    """Resto de Lagrange (x > 0): (x + n + 1) * e^x * x^(n+1) / (n+1)!"""
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


# Tabela de valores
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

    with open(os.path.join(PASTA, "tabela_valores.csv"), "w", newline="", encoding="utf-8") as arq:
        w = csv.writer(arq)
        w.writerow(cab)
        w.writerows(linhas)


# Gráficos
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
    fig.savefig(os.path.join(PASTA, "serie_taylor.png"))
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
    fig.savefig(os.path.join(PASTA, "grafico_erros.png"))
    plt.close(fig)


def medir_tempo(funcao, repeticoes=5, chamadas=2000):
    """Tempo médio por chamada (em microssegundos): melhor de 'repeticoes' rodadas."""
    t = min(timeit.repeat(funcao, number=chamadas, repeat=repeticoes))
    return t / chamadas * 1e6


def medir_tempos(ns, x=1.0):
    """Mede o tempo da série padrão, da otimizada (sem parada antecipada) e do math.exp."""
    t_taylor = [medir_tempo(lambda: taylor(x, n)) for n in ns]
    t_otim = [medir_tempo(lambda: taylor_otimizada(x, n, 0)) for n in ns]
    t_exata = medir_tempo(lambda: f(x), chamadas=20000)
    return t_taylor, t_otim, t_exata


def grafico_tempo():
    ns = list(range(1, N_LIMITE + 1))
    t_taylor, t_otim, t_exata = medir_tempos(ns)

    print("=" * 64)
    print("VELOCIDADE DE PROCESSAMENTO (x = 1, tempo por chamada)")
    print("=" * 64)
    print(f"x * math.exp(x) (biblioteca): {t_exata:.3f} µs")
    print(f"{'termos':>7} | {'padrão (µs)':>11} | {'otimizada (µs)':>14} | {'ganho':>6}")
    print("-" * 49)
    for n in (3, 5, 10, 20, N_LIMITE):
        print(f"{n:>7} | {t_taylor[n - 1]:>11.3f} | {t_otim[n - 1]:>14.3f} | "
              f"{t_taylor[n - 1] / t_otim[n - 1]:>5.1f}x")
    print()

    fig, ax = plt.subplots(figsize=(8, 5), dpi=120)
    ax.plot(ns, t_taylor, color=AZUL, linewidth=2, marker="o", markersize=6,
            markeredgecolor=SUPERFICIE, markeredgewidth=1.5, label="Taylor padrão")
    ax.plot(ns, t_otim, color=AQUA, linewidth=2, marker="o", markersize=6,
            markeredgecolor=SUPERFICIE, markeredgewidth=1.5, label="Taylor otimizada")
    ax.axhline(t_exata, color=LARANJA, linewidth=2, label="x · math.exp(x) (biblioteca)")

    ax.set_ylim(bottom=0)
    ax.set_xlabel("Número de termos da série")
    ax.set_ylabel("Tempo por chamada (µs)")
    ax.set_title("Velocidade de processamento", loc="left", color=TEXTO, fontsize=12)
    ax.grid(True, axis="y", color=GRADE, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(PASTA, "grafico_tempo.png"))
    plt.close(fig)

    grafico_erro_vs_tempo(ns, t_taylor, t_otim, t_exata)


def grafico_erro_vs_tempo(ns, t_taylor, t_otim, t_exata, x=1.0):
    """Erro da função x tempo de execução (cada ponto é um valor de N)."""
    erros = [max(abs(f(x) - taylor(x, n)), 1e-17) for n in ns]

    fig, ax = plt.subplots(figsize=(8, 5), dpi=120)
    ax.plot(t_taylor, erros, color=AZUL, linewidth=2, marker="o", markersize=6,
            markeredgecolor=SUPERFICIE, markeredgewidth=1.5, label="Taylor padrão")
    ax.plot(t_otim, erros, color=AQUA, linewidth=2, marker="o", markersize=6,
            markeredgecolor=SUPERFICIE, markeredgewidth=1.5, label="Taylor otimizada")
    ax.axvline(t_exata, color=LARANJA, linestyle="--", linewidth=2,
               label="tempo do x · math.exp(x)")

    ax.set_yscale("log")
    ax.set_xlabel("Tempo de execução por chamada (µs)")
    ax.set_ylabel("Erro absoluto (escala log)")
    ax.set_title(f"Erro × tempo de execução (x = {x:g}, N = 1 a {ns[-1]})",
                 loc="left", color=TEXTO, fontsize=12)
    ax.grid(True, axis="y", color=GRADE, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(PASTA, "grafico_erro_vs_tempo.png"))
    plt.close(fig)


if __name__ == "__main__":
    mostrar_derivadas()
    mostrar_n_limite()
    propriedades()
    tabela_valores()
    grafico_aproximacoes()
    grafico_erros()
    grafico_tempo()
    print("Arquivos gerados: tabela_valores.csv, serie_taylor.png, "
          "grafico_erros.png, grafico_tempo.png, grafico_erro_vs_tempo.png")
    print(f"Pasta onde foram salvos: {os.path.abspath(PASTA)}")

    # Pausa opcional para a janela do terminal não fechar sozinha:
    # input("\nPressione Enter para sair...")
