# Testes automatizados — Rigidez

## 1. Objetivo

A suíte de testes automatizados do Rigidez utiliza o `pytest` para verificar a análise matricial de estruturas planas, o pós-processamento das solicitações internas e dos deslocamentos e a consistência das combinações de carregamento. Os testes analíticos utilizam soluções da mecânica estrutural; os testes de regressão verificam propriedades que devem permanecer válidas após alterações no código.

São 14 funções de teste, correspondentes a 18 casos executados (uma função é parametrizada com quatro entradas).

## 2. Organização

```text
rigidez/
├── src/
├── tests/
│   ├── README.md
│   ├── conftest.py
│   ├── test_estrutura.py
│   ├── test_diagramas.py
│   └── test_combinacoes.py
├── pytest.ini
└── requirements-dev.txt
```

## 3. Instalação e execução

Abra a pasta raiz do projeto no VS Code. No terminal PowerShell, utilize o mesmo ambiente Python em que o Rigidez está instalado. Se necessário, crie e ative um ambiente virtual:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

Instale as dependências do projeto (caso ainda não estejam instaladas) e as dependências de desenvolvimento:

```powershell
python -m pip install -r requirements.txt
python -m pip install -r requirements-dev.txt
```

Execute a suíte a partir da raiz do projeto:

```powershell
python -m pytest
```

| Comando | Finalidade |
| --- | --- |
| `python -m pytest -v` | Exibir o nome e o resultado de cada caso. |
| `python -m pytest tests/test_estrutura.py -v` | Executar os testes de análise matricial. |
| `python -m pytest tests/test_diagramas.py -v` | Executar os testes de pós-processamento. |
| `python -m pytest tests/test_combinacoes.py -v` | Executar os testes de combinações. |
| `python -m pytest tests/test_estrutura.py::test_balanco -v` | Executar um teste específico. |
| `python -m pytest -k restauracao -v` | Selecionar testes pelo nome. |
| `python -m pytest -x` | Interromper na primeira falha. |
| `python -m pytest --collect-only -q` | Listar os casos sem executá-los. |
| `python -m pytest --lf -v` | Reexecutar os casos que falharam na última execução. |

Resultado esperado para a versão validada:

```text
.................. [100%]
18 passed
```

O tempo de execução depende do computador. `PASSED` indica aprovação, `FAILED` indica falha em uma verificação e `ERROR` indica erro durante a preparação ou execução do teste. Se houver falha, examine o traceback e os valores esperado e obtido antes de modificar as tolerâncias ou o código.

## 4. Modelos auxiliares (`conftest.py`)

As fixtures do pytest fornecem modelos reutilizáveis aos testes. As unidades dos modelos estruturais são cm e kN (tensões e módulo de elasticidade em kN/cm²).

- `balanco`: barra horizontal de 300 cm, engastada à esquerda, seção retangular de 20 × 40 cm e módulo de elasticidade de 3000 kN/cm².
- `biapoiada`: viga de 600 cm, dois elementos de 300 cm, mesma seção e módulo de elasticidade, carga uniforme de 0,2 kN/cm. O apoio esquerdo restringe as translações horizontal e vertical; o direito restringe a vertical.
- `modelo`: viga de 600 cm, largura 30 cm, altura 80 cm, uma abertura (`lf=250`, `b=100`, `a=30`, `c=25`, em cm), apoios de 40 cm, `fck=35`, `fyk=fywk=500` (MPa), espaçamento máximo de 30 cm. Cargas de 0,15 e 0,10 kN/cm com fatores 1,4 e 1,5, respectivamente. A fixture é compartilhada no módulo de testes de combinações.

A função auxiliar `analisar(portico)` monta o vetor de forças nodais, aplica as condições de contorno ao vetor, resolve os deslocamentos, calcula as solicitações nodais e as reações. A rigidez e as cargas de barra devem estar previamente definidas.

## 5. Análise estrutural (`test_estrutura.py`) — 7 casos

### 01 — `test_tracao`

**Objetivo:** verificar o comportamento axial de uma barra. Aplica força horizontal de 50 kN na extremidade livre do balanço. A solução analítica é

$$u_B=\frac{FL}{EA}.$$

Verifica o deslocamento horizontal `u[3]`, a reação horizontal de −50 kN e o vetor local de forças nodais `fl = [-50, 0, 0, 50, 0, 0]`. A tolerância absoluta da comparação do vetor é $10^{-8}$.

### 02 — `test_balanco`

**Objetivo:** verificar flexão e equilíbrio de uma viga em balanço. Aplica força vertical descendente de 20 kN na extremidade livre. Pela teoria de Euler–Bernoulli:

$$v_B=-\frac{FL^3}{3EI},\qquad \theta_B=-\frac{FL^2}{2EI}.$$

Verifica `u[4]`, `u[5]`, a reação vertical $R_A=F=20$ kN e o momento de reação $M_A=FL=6000$ kN·cm.

### 03 — `test_biapoiada`

**Objetivo:** verificar a resposta de uma viga biapoiada com carga uniforme. Com dois elementos finitos, o nó central permite comparar diretamente a flecha nodal com a solução analítica:

$$v_C=-\frac{5qL^4}{384EI},\qquad R_A=R_B=\frac{qL}{2}=60\ \mathrm{kN}.$$

Verifica `u[4]`, as reações `R[1]` e `R[7]` e o equilíbrio vertical $R_A+R_B=qL=120$ kN. Não verifica a flecha em todos os pontos entre nós.

### 04 — `test_matriz_original`

**Objetivo:** garantir que a aplicação das condições de contorno não altere a matriz global original $\mathbf K$. Guarda uma cópia de `k`, chama `monta_k01()` e verifica que `k` permanece igual à cópia, é simétrica e não compartilha memória com `k01`. Também verifica que as três primeiras linhas e colunas de `k01`, relativas ao engaste, correspondem às da matriz identidade de ordem seis.

### 05 — `test_remontagem`

**Objetivo:** impedir a acumulação indevida das contribuições elementares. Guarda a matriz global original e executa `monta_k()` três vezes. Após cada montagem, a matriz deve ser igual à original, pois as propriedades da estrutura não foram alteradas.

### 06 — `test_reacoes_repetidas`

**Objetivo:** impedir a acumulação de reações. Guarda o vetor de reações da viga biapoiada e chama `calcula_reacoes()` três vezes. Cada resultado deve coincidir com o vetor original (tolerância absoluta de $10^{-8}$), sem somar as reações das chamadas anteriores.

### 07 — `test_residuo`

**Objetivo:** verificar a solução do sistema linear com condições de contorno. Compara `k01 @ u` com `fnos`, utilizando tolerância absoluta de $10^{-8}$:

$$\mathbf K_{01}\mathbf u-\mathbf F_{01}\approx\mathbf 0.$$

Este teste verifica o resíduo algébrico do sistema resolvido, não substituindo a validação física por soluções analíticas.

## 6. Pós-processamento (`test_diagramas.py`) — 4 casos

### 08 — `test_deformada_balanco`

**Objetivo:** validar a interpolação cúbica de Hermite para o balanço com força descendente de 20 kN na ponta. Compara os deslocamentos transversais em 21 pontos com

$$v(x)=-\frac{Fx^2(3L-x)}{6EI}.$$

Verifica também que o deslocamento axial é nulo em todos os pontos. A tolerância absoluta é $10^{-10}$. O caso possui solução cúbica exata segundo Euler–Bernoulli; não demonstra exatidão da interpolação para qualquer carregamento.

### 09 — `test_diagrama_uniforme`

**Objetivo:** validar os diagramas de esforço normal, cortante e momento fletor da viga biapoiada sob carga uniforme. Para a coordenada global $X=x+x_i$, em que $x$ é a coordenada local da barra:

$$N(X)=0,$$
$$V(X)=q(X-L/2),$$
$$M(X)=\frac{qX(L-X)}{2}.$$

Compara todos os pontos retornados por `calcula_solicitacoes_internas(cargas)` em ambas as barras, com tolerância absoluta de $10^{-8}$. As expressões seguem a convenção de sinais dos vetores retornados pelo programa.

### 10 — `test_salto_cortante`

**Objetivo:** verificar a descontinuidade do cortante sob uma força concentrada de 20 kN aplicada a 120 cm do engaste. O algoritmo deve retornar dois registros na posição da carga, representando os limites à esquerda e à direita. Verifica que existem exatamente dois registros, que $V^+-V^-=20$ kN e que $M^+=M^-$ (tolerância absoluta de $10^{-8}$ para o momento).

### 11 — `test_deformada_pontos_invalidos`

**Objetivo:** validar a entrada do método `calcula_deformada()`. A chamada com apenas um ponto deve lançar `ValueError`, pois são necessários pelo menos dois pontos para representar as extremidades da barra.

## 7. Combinações de carregamento (`test_combinacoes.py`) — 7 casos

### 12 — `test_combinacoes_lineares`

**Objetivo:** verificar a proporcionalidade dos resultados de uma análise linear quando as cargas uniformes são majoradas. No modelo auxiliar:

$$w_k=0,15+0,10=0,25\ \mathrm{kN/cm},$$
$$w_d=1,4\times0,15+1,5\times0,10=0,36\ \mathrm{kN/cm}.$$

O fator de proporcionalidade é $w_d/w_k=1,44$. Verifica que os vetores globais `u`, `R` e `fnos`, assim como `fl` de cada barra, da combinação de cálculo correspondem aos característicos multiplicados por esse fator (tolerância absoluta de $10^{-7}$). A proporcionalidade é válida neste caso de carregamentos uniformes proporcionais e rigidez constante, não para combinações arbitrárias.

### 13 — `test_restauracao`

**Objetivo:** verificar a restauração completa dos estados de carregamento sem modificar o dimensionamento. Alterna entre `caracteristicos`, `calculo` e `caracteristicos`. Em cada seleção verifica:

- Identificador da combinação e lista de carregamentos ativos.
- Vetores globais `u`, `R`, `fnos` e `f` e vetores de cada barra `u`, `fl`, `f`, `fep` e `fepl` iguais aos armazenados, sem compartilhamento de memória NumPy.
- Preservação de `ASL`, `ASL2` e `ASW_S`.
- Após novo cálculo, igualdade de `R`, `u` e `fl` com os resultados armazenados (tolerâncias absolutas de $10^{-7}$, $10^{-9}$ e $10^{-7}$, respectivamente).

Ao final, um bloco `finally` restaura a combinação de cálculo, inclusive se ocorrer uma falha.

### 14 — `test_combinacao_invalida`

**Objetivo:** verificar que `selecionar_esforcos('invalido')` lança `ValueError`. Evita a seleção silenciosa de uma combinação desconhecida.

### 15 — `test_cargas_invalidas[cargas0]`

**Entrada:** `[]`. Verifica que `normalizar_cargas()` rejeita uma lista vazia de carregamentos, lançando `ValueError`.

### 16 — `test_cargas_invalidas[cargas1]`

**Entrada:** `[{'w': 1, 'gamaf': -1}]`. Verifica que um coeficiente de majoração negativo é rejeitado com `ValueError`.

### 17 — `test_cargas_invalidas[cargas2]`

**Entrada:** `[{'w': float('nan'), 'gamaf': 1}]`. Verifica que uma intensidade de carga não finita (`NaN`) é rejeitada com `ValueError`.

### 18 — `test_cargas_invalidas[cargas3]`

**Entrada:** `[{'w': 1}]`. Verifica que a ausência do coeficiente `gamaf` é rejeitada com `ValueError`.

Os quatro últimos casos são gerados por `@pytest.mark.parametrize` a partir da mesma função `test_cargas_invalidas()`.

## 8. Comparações numéricas e tolerâncias

Os testes utilizam `pytest.approx()` para comparações escalares, `numpy.testing.assert_allclose()` para vetores e matrizes, `numpy.shares_memory()` para verificar independência de armazenamento e `pytest.raises()` para verificar exceções esperadas.

Em `assert_allclose`, a comparação considera $|a-b|\leq atol+rtol\,|b|$, em que $b$ é o valor de referência. Quando o código informa apenas `atol`, a tolerância relativa padrão continua sendo aplicada. `pytest.approx()` utiliza suas próprias tolerâncias padrão quando não especificadas. Pequenas diferenças de arredondamento não devem ser confundidas com erros de formulação.

## 9. Escopo e limitações

| Grupo | Casos |
| --- | ---: |
| Análise estrutural | 7 |
| Pós-processamento | 4 |
| Combinações de carregamento | 7 |
| **Total** | **18** |

A aprovação dos 18 casos demonstra consistência nos cenários testados, mas não valida integralmente o software. Ainda são necessários testes específicos para dimensionamento conforme a NBR 6118, ancoragem, escalonamento, detalhamento, equivalência entre modelos com e sem aberturas e interface gráfica. A preservação das áreas de aço na troca de combinação não equivale à validação das fórmulas de dimensionamento.

Recomenda-se executar a suíte antes e depois de cada alteração relevante no código e antes de realizar um commit ou publicar uma nova versão.
