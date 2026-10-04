## EP 1 – Cálculos Complexos com Séries de Taylor  
# Tema: Cálculos Complexos com Séries de Taylor – função y = x · eˣ (exponencial × x)  

Integrantes do grupo: Bernardo Pereira Kubo

## Instalação, execução e exemplos

**Instalação** (Python 3.8+):

```bash
pip install -r requirements.txt
```

**Execução:**

```bash
python serie_taylor_completo.py
```

O programa não pede entrada. Ele imprime os resultados no terminal e salva os gráficos e o CSV na pasta `resultados/`.

**Exemplo de entrada:**

```python
from serie_taylor_completo import f, taylor, taylor_otimizada

print(f(1.0))
print(taylor(1.0, 10))
print(taylor_otimizada(1.0))
```

**Saída:**

```
2.718281828459045
2.7182815255731922
2.7182818284590455
```
