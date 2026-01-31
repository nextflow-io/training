# Parte 3: Organizando entradas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Na Parte 2, executamos o molkart com múltiplos parâmetros na linha de comando.
Agora vamos aprender duas abordagens melhores para gerenciar entradas: **arquivos de parâmetros** e **planilhas de amostras**.

## 1. Usando arquivos de parâmetros

### 1.1. O problema com linhas de comando longas

Lembre-se do nosso comando da Parte 2:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose" \
  --outdir results
```

Isso funciona, mas é difícil de reproduzir, compartilhar ou modificar.
E se você precisar executar a mesma análise novamente no próximo mês?
E se um colaborador quiser usar suas configurações exatas?

### 1.2. Solução: Use um arquivo de parâmetros

Crie um arquivo chamado `params.yaml`:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

Agora seu comando se torna:

```bash
nextflow run ./molkart -params-file params.yaml -resume
```

É isso! O arquivo de parâmetros documenta sua configuração exata e facilita a reexecução ou compartilhamento.

### 1.3. Sobrescrevendo parâmetros

Você ainda pode sobrescrever parâmetros específicos da linha de comando:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "stardist" --outdir stardist_results -resume
```

A linha acima altera o `segmentation_method` para `stardist` e o nome do `--outdir` para `stardist_results` em vez dos parâmetros no arquivo `params.yaml`.
Além disso, você pode ver que a flag `-resume` nos permitiu reutilizar os resultados de pré-processamento da execução anterior, economizando tempo.
Você pode usar esse padrão para testar rapidamente diferentes variações do fluxo de trabalho.

### Conclusão

Arquivos de parâmetros tornam suas análises reproduzíveis e fáceis de compartilhar.
Use-os para qualquer trabalho de análise real.

### O que vem a seguir?

Aprenda como as planilhas de amostras organizam informações sobre múltiplas amostras.

---

## 2. O padrão de planilha de amostras

### 2.1. O que é uma planilha de amostras?

Uma planilha de amostras é um arquivo CSV que descreve suas amostras de entrada.
Cada linha é uma amostra, e as colunas especificam os arquivos e metadados para aquela amostra.

Vamos olhar a planilha de amostras que estamos usando:

```bash
cat data/samplesheet.csv
```

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
```

As colunas são:

- `sample`: Identificador único da amostra
- `nuclear_image`: Imagem de coloração nuclear (TIFF)
- `spot_table`: Pontos de transcrição (TXT)
- `membrane_image`: Imagem de coloração de membrana (TIFF, opcional)

### 2.2. Caminhos de arquivos

Planilhas de amostras aceitam múltiplos tipos de caminho:

- **URLs**: Nextflow baixa automaticamente (como mostrado acima)
- **Caminhos locais**: `data/nuclear.tiff` ou `/absolute/path/to/nuclear.tiff`
- **Armazenamento em nuvem**: `s3://bucket/nuclear.tiff`, `gs://bucket/nuclear.tiff`, `az://container/nuclear.tiff`

Você pode misturar tipos de caminho na mesma planilha de amostras.

### 2.3. Criando sua própria planilha de amostras

Primeiro, vamos baixar os arquivos de dados de teste localmente:

```bash
cd /workspaces/training/nf4-science/imaging/data
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
cd ..
```

Agora vamos modificar a planilha de amostras para referenciar esses arquivos locais:

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,data/nuclear.tiff,data/spots.txt,data/membrane.tiff
```

!!! warning "Aviso"

    Observe que os caminhos na planilha de amostras são relativos a onde você **executa** o Nextflow, não onde a planilha de amostras está localizada.

Finalmente, vamos executar o nf-core/molkart mais uma vez com a planilha de amostras com caminhos de arquivos locais:

`nextflow run ./molkart -params-file params.yaml -resume`

Como você pode ver, o Nextflow executa essa execução de forma similar a quando os arquivos foram baixados do Github. Esta é uma das grandes funcionalidades do Nextflow, ele prepara os dados adequadamente para você, independentemente de onde estão localizados.

### Conclusão

Planilhas de amostras organizam conjuntos de dados com múltiplas amostras de uma forma que permite definir explicitamente seus metadados junto com os caminhos dos arquivos.
A maioria dos fluxos de trabalho do nf-core usa esse padrão.

### O que vem a seguir?

Agora que cobrimos as entradas, vamos explorar como configurar fluxos de trabalho Nextflow para diferentes ambientes computacionais.
