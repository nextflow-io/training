# Resumo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Você concluiu o treinamento de Desenvolvimento de Plugins.
Esta página recapitula o que você construiu em cada parte, aborda a distribuição e fornece orientações sobre os próximos passos.

---

## O que você aprendeu

### Parte 1: Usando plugins

Você descobriu como os plugins do Nextflow funcionam do ponto de vista do usuário.
Você instalou nf-schema e nf-co2footprint, os configurou via `nextflow.config`, e viu como os plugins podem validar entradas, adicionar funções e se conectar aos eventos do ciclo de vida do pipeline.

### Parte 2: Configurando o ambiente

Você configurou um ambiente de desenvolvimento de plugins com Java 21+, criou um novo projeto de plugin usando o comando `nextflow plugin create`, e aprendeu a estrutura de projeto que o Nextflow espera: arquivos de código-fonte, configuração de build e o fluxo de trabalho com o Makefile.

### Parte 3: Funções personalizadas

Você implementou seu primeiro ponto de extensão criando métodos anotados com `@Function` em uma classe `PluginExtensionPoint`.
Você construiu `reverseGreeting` e `decorateGreeting`, depois os importou e chamou a partir de um script de pipeline.

### Parte 4: Testes

Você escreveu testes unitários para suas funções personalizadas usando o framework de testes do Groovy.
Você aprendeu a executar os testes com `make test` e verificar que seu plugin se comporta corretamente antes de instalá-lo.

### Parte 5: Observers

Você implementou a interface `TraceObserver` para se conectar aos eventos do ciclo de vida do pipeline.
Você construiu `GreetingObserver` (reagindo ao início e à conclusão do pipeline) e `TaskCounterObserver` (contando tarefas concluídas), depois os registrou por meio de uma `TraceObserverFactory`.

### Parte 6: Configuração

Você tornou seu plugin configurável via `nextflow.config` usando `session.config.navigate()` para ler valores em tempo de execução.
Você adicionou uma classe `@ConfigScope` para declarar formalmente as opções do seu plugin, eliminando os avisos de "Unrecognized config option" e habilitando o suporte a IDEs.

---

## Distribuição

Assim que seu plugin estiver funcionando localmente, você pode compartilhá-lo com outras pessoas por meio do registro de plugins do Nextflow.

### Versionamento

Siga o [versionamento semântico](https://semver.org/) para seus lançamentos:

| Mudança de versão         | Quando usar                              | Exemplo                                              |
| ------------------------- | ---------------------------------------- | ---------------------------------------------------- |
| **MAJOR** (1.0.0 → 2.0.0) | Mudanças incompatíveis                   | Remover uma função, alterar tipos de retorno         |
| **MINOR** (1.0.0 → 1.1.0) | Novos recursos, compatíveis com versões anteriores | Adicionar uma nova função               |
| **PATCH** (1.0.0 → 1.0.1) | Correções de bugs, compatíveis com versões anteriores | Corrigir um bug em uma função existente |

Atualize a versão em `build.gradle` antes de cada lançamento:

```groovy title="build.gradle"
version = '1.0.0'  // Use versionamento semântico: MAJOR.MINOR.PATCH
```

### Publicando no registro

O [registro de plugins do Nextflow](https://registry.nextflow.io/) é a forma oficial de compartilhar plugins com a comunidade.

O fluxo de trabalho de publicação:

1. Reivindique o nome do seu plugin no [registro](https://registry.nextflow.io/) (faça login com sua conta do GitHub)
2. Configure suas credenciais de API em `~/.gradle/gradle.properties`
3. Execute os testes para verificar que tudo funciona: `make test`
4. Publique com `make release`

Para instruções passo a passo, consulte a [documentação oficial de publicação](https://www.nextflow.io/docs/latest/guides/gradle-plugin.html#publishing-a-plugin).

Após a publicação, os usuários instalam seu plugin sem nenhuma configuração local:

```groovy title="nextflow.config"
plugins {
    id 'nf-greeting@1.0.0'
}
```

O Nextflow baixa automaticamente o plugin do registro no primeiro uso.

---

## Lista de verificação para desenvolvimento de plugins

- [ ] Java 21+ instalado
- [ ] Criar projeto com `nextflow plugin create <name> <org>`
- [ ] Implementar classe de extensão com métodos `@Function`
- [ ] Escrever testes unitários e executá-los com `make test`
- [ ] Fazer build e instalar com `make install`
- [ ] Opcionalmente, adicionar implementações de `TraceObserver` para eventos do fluxo de trabalho
- [ ] Opcionalmente, adicionar `ConfigScope` para configuração do plugin
- [ ] Habilitar em `nextflow.config` com `plugins { id 'plugin-id' }`
- [ ] Importar funções com `include { fn } from 'plugin/plugin-id'`
- [ ] Versionar e publicar no registro

---

## Padrões de código principais

**Definição de função:**

```groovy
@Function
String myFunction(String input, String optional = 'default') {
    return input.transform()
}
```

**Configuração do plugin:**

```groovy
nextflowPlugin {
    provider = 'my-org'
    className = 'my.org.MyPlugin'
    extensionPoints = ['my.org.MyExtension']
}
```

**Usando em fluxos de trabalho:**

```groovy
include { myFunction } from 'plugin/my-plugin'

workflow {
    channel.of('a', 'b', 'c')
        .map { item -> myFunction(item) }
        .view()
}
```

---

## Resumo dos pontos de extensão

| Tipo                | Classe/Anotação  | Finalidade                                                    |
| ------------------- | ---------------- | ------------------------------------------------------------- |
| Função              | `@Function`      | Pode ser chamada a partir de fluxos de trabalho               |
| Trace Observer      | `TraceObserver`  | Conecta-se aos eventos do ciclo de vida do fluxo de trabalho  |
| Escopo de configuração | `@ScopeName`  | Define a configuração do plugin no nextflow.config            |

---

## O que fazer a seguir

Aqui estão alguns próximos passos práticos para continuar sua jornada de desenvolvimento de plugins.

**Construa algo real.**
Escolha um caso de uso do seu próprio trabalho: uma função personalizada que sua equipe usa repetidamente, um observer que envia notificações no Slack quando o pipeline é concluído, ou um escopo de configuração que padroniza opções nos pipelines da sua organização.
Começar a partir de um problema real é a forma mais rápida de aprofundar seu entendimento.

**Use nf-hello como referência.**
O repositório [nf-hello](https://github.com/nextflow-io/nf-hello) é o exemplo oficial mínimo de plugin.
É um bom ponto de partida para novos projetos e uma referência útil quando você precisa verificar como algo está estruturado.

**Leia a documentação oficial.**
A documentação do Nextflow abrange tópicos além deste treinamento, incluindo fábricas de canais, sobrecarga de operadores e padrões avançados de observers.
O guia [developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html) é a referência mais abrangente.

**Estude plugins existentes.**
O [repositório de plugins do Nextflow](https://github.com/nextflow-io/plugins) contém o código-fonte de plugins oficiais como nf-schema, nf-wave e nf-tower.
Ler o código de plugins em produção é uma das melhores formas de aprender padrões e convenções que vão além dos exemplos introdutórios.

---

## Recursos adicionais

**Documentação oficial:**

- [Using plugins](https://www.nextflow.io/docs/latest/plugins/plugins.html): guia abrangente para instalar e configurar plugins
- [Developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html): referência detalhada para desenvolvimento de plugins
- [Config scopes](https://nextflow.io/docs/latest/developer/config-scopes.html): criando escopos de configuração para plugins

**Descoberta de plugins:**

- [Nextflow Plugin Registry](https://registry.nextflow.io/): navegue e descubra plugins disponíveis
- [Plugin registry docs](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html): documentação do registro

**Exemplos e referências:**

- [nf-hello](https://github.com/nextflow-io/nf-hello): plugin de exemplo simples (ótimo ponto de partida)
- [Nextflow plugins repository](https://github.com/nextflow-io/plugins): coleção de plugins oficiais para referência
