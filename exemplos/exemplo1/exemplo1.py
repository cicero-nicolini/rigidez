"""Compatibilidade com o caminho anterior do Exemplo 1.

Executa o script original, agora organizado em analise_estrutural.
"""
from pathlib import Path
import runpy

CAMINHO_EXEMPLO = (
    Path(__file__).resolve().parents[1]
    / "analise_estrutural"
    / "exemplo1"
    / "exemplo1.py"
)
runpy.run_path(str(CAMINHO_EXEMPLO), run_name="__main__")
