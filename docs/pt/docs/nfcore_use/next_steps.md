# Resumo do curso

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Parabéns por concluir o curso de treinamento Use nf-core! 🎉

<!-- placeholder for video -->

## Sua jornada

Você começou encontrando e recuperando o pipeline `nf-core/demo`, depois aprendeu a executá-lo usando seu perfil de teste e a examinar suas saídas.
Em seguida, você configurou sua execução por meio de parâmetros de pipeline e arquivos de configuração, e viu como os pipelines nf-core validam parâmetros e dados de entrada.
Por fim, você aplicou essas mesmas habilidades ao `nf-core/rnaseq`, um pipeline em escala de produção, e aprendeu como substituir suas alocações de recursos padrão para se adequar ao hardware disponível para você.

### O que você aprendeu

Agora você é capaz de encontrar, recuperar, executar e configurar pipelines nf-core.

- Os pipelines nf-core são recuperados com `nextflow pull` e seguem uma organização de código padronizada.
- Todo pipeline nf-core vem com um perfil `test` para validação rápida em um conjunto de dados pequeno.
- Os parâmetros de pipeline (definidos via `--param_name` ou `-params-file`) e a configuração (definida via `-c`) têm finalidades diferentes: entradas e opções de análise versus logística de execução, como alocação de recursos.
- Os pipelines nf-core validam parâmetros e arquivos de entrada automaticamente, detectando erros antes que qualquer trabalho seja realizado.
- Os recursos padrão são atribuídos por meio de labels (`process_low`, `process_medium`, `process_high`) definidos em `conf/base.config`, que você pode substituir com um arquivo de configuração personalizado.

### Habilidades adquiridas

Ao longo deste curso prático, você aprendeu como:

- Encontrar um pipeline nf-core no site nf-co.re e recuperar seu código-fonte
- Executar um pipeline usando seu perfil de teste integrado e examinar suas saídas
- Obter ajuda, definir parâmetros e entender a validação de parâmetros e entradas
- Personalizar a alocação de recursos e os argumentos de ferramentas por meio de arquivos de configuração
- Recuperar e executar um pipeline em escala de produção, e substituir seus labels de recursos padrão

Agora você está equipado com o conhecimento fundamental para começar a executar pipelines nf-core em suas próprias análises.

## Próximos passos para desenvolver suas habilidades

Aqui estão nossas principais sugestões sobre o que fazer a seguir:

- Lance e monitore esses pipelines em escala com [Scale with Seqera](../seqera_scale/index.md)
- Não apenas execute pipelines nf-core, desenvolva-os! Aprenda as boas práticas do nf-core com [Build with nf-core](../nfcore_build/index.md)
- Novo no Nextflow? Comece com [Nextflow Run](../nextflow_run/index.md)
- Aplique o Nextflow a um caso de uso de análise científica com [Nextflow for Science](../nf4_science/index.md)
- Explore recursos mais avançados do Nextflow com os [Side Quests](../side_quests/index.md)

## Obtendo ajuda

Para recursos de ajuda e suporte da comunidade, consulte a [página de Ajuda](../help.md).

## Pesquisa de feedback

Antes de continuar, reserve um minuto para responder à pesquisa do curso! Seu feedback nos ajuda a melhorar nossos materiais de treinamento para todos.

[Responder à pesquisa :material-arrow-right:](survey.md){ .md-button .md-button--primary }
