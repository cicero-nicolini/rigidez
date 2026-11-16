# Exemplo 1 — Análise matricial de pórtico

**Origem:** planilha da disciplina de Análise Matricial do professor Quevedo, conforme o cabeçalho do script original.

**Objetivo:** determinar deslocamentos, solicitações e reações de um pórtico.

## Dados implementados

- Nós: 1 `(0, 75)`, 2 `(100, 75)` e 3 `(200, 0)`.
- Restrições: translações e rotação impedidas nos nós 1 e 3.
- Cargas nodais no nó 2: `Fy = -10` e `Mz = -1000`.
- Seção: `Secao("retangular", 0.29, 34.64, 10000.0)`.
- Barras: 1–2 e 2–3.
- Carregamento distribuído na barra 1: `Carregamento_distribuido(0, 100, 0.24, 0.24, barra1)`.
- Carregamento pontual na barra 2: `Carregamento_pontual(62.5, 12.0, 16.0, barra2)`.

Os valores e unidades seguem exatamente o código existente; nenhuma conversão foi aplicada.

## Procedimento computacional

1. Define nós, vínculos, seção, barras e cargas.
2. Monta as matrizes `k` e `k01`.
3. Calcula forças de engastamento perfeito e o vetor de cargas nodais.
4. Calcula deslocamentos e solicitações internas nas barras.
5. Calcula reações de apoio.

## Arquivos

- `exemplo1.py`: cálculo original, preservado sem alterações.
- `exemplo1.ftl`: modelo de referência para o FTOOL, preservado sem alterações.

## Execução

No terminal, a partir da raiz do projeto:

```powershell
python exemplos/analise_estrutural/exemplo1/exemplo1.py
```

**Dependências:** módulos do pacote `src` do próprio projeto e suas dependências.

## Validação

O exemplo imprime matrizes, vetores de deslocamentos, solicitações e reações. Os resultados numéricos de referência e as tolerâncias de comparação ainda não estão documentados; nenhuma equivalência com o FTOOL foi presumida.
