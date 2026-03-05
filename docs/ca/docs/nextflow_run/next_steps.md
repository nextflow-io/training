# Resum del curs

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Felicitats per haver completat el curs de formació Nextflow Run! 🎉

<!-- placeholder for video -->

## El teu recorregut

Has començat amb un workflow molt bàsic, i has après a executar-lo, trobar les sortides i gestionar la seva execució.
Després, has treballat amb versions cada vegada més complexes d'aquest workflow i has après a reconèixer els conceptes i mecanismes essencials que impulsen els pipelines de Nextflow, incloent canals i operadors, modularització de codi i contenidors.
Finalment, has après com personalitzar la configuració d'un pipeline per adaptar-lo a les teves preferències i la teva infraestructura computacional.

### Què has après

Ara ets capaç de gestionar l'execució del pipeline Hello, descriure com està estructurat i identificar les peces principals de codi implicades.

- La forma final del workflow Hello pren com a entrada un fitxer CSV que conté salutacions de text.
- Els quatre passos s'implementen com a processos de Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` i `cowpy`) emmagatzemats en fitxers de mòdul separats.
- Els resultats es publiquen en un directori anomenat `results/`.
- La sortida final del pipeline és un fitxer de text pla que conté art ASCII d'un personatge dient les salutacions en majúscules.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Escriu cada salutació al seu propi fitxer de sortida (_p. ex._ "Hello-output.txt")
2. **`convertToUpper`:** Converteix cada salutació a majúscules (_p. ex._ "HELLO")
3. **`collectGreetings`:** Recull totes les salutacions en majúscules en un únic fitxer per lots
4. **`cowpy`:** Genera art ASCII utilitzant l'eina `cowpy`

La configuració del workflow permet proporcionar entrades i paràmetres de manera flexible i reproduïble.

### Habilitats adquirides

A través d'aquest curs pràctic, has après com:

- Executar un workflow de Nextflow localment
- Trobar i interpretar sortides (resultats) i fitxers de registre generats per Nextflow
- Reconèixer els components bàsics de Nextflow que constitueixen un workflow simple de múltiples passos
- Descriure conceptes del següent nivell com operadors i factories de canals
- Configurar pipelines per a diferents entorns computacionals

Ara estàs equipat amb el coneixement fonamental per començar a integrar pipelines de Nextflow existents en el teu propi treball.

## Següents passos per desenvolupar les teves habilitats

Aquí tens les nostres principals suggerències sobre què fer a continuació:

- No només executis Nextflow, escriu-lo! Converteix-te en desenvolupador de Nextflow amb [Hello Nextflow](../hello_nextflow/index.md)
- Aplica Nextflow a un cas d'ús d'anàlisi científica amb [Nextflow for Science](../nf4_science/index.md)
- Comença amb nf-core amb [Hello nf-core](../hello_nf-core/index.md)
- Aprèn tècniques de resolució de problemes amb la [Missió Secundària de Depuració](../side_quests/debugging.md)

Finalment, et recomanem que donis una ullada a [**Seqera Platform**](https://seqera.io/), una plataforma basada en el núvol desenvolupada pels creadors de Nextflow que fa encara més fàcil executar i gestionar els teus workflows, així com gestionar les teves dades i executar anàlisis de manera interactiva en qualsevol entorn.

## Obtenir ajuda

Per a recursos d'ajuda i suport de la comunitat, consulta la [pàgina d'Ajuda](../help.md).

## Enquesta de valoració

Abans de continuar, si us plau, dedica un minut a completar l'enquesta del curs! La teva valoració ens ajuda a millorar els nostres materials de formació per a tothom.

[Fes l'enquesta :material-arrow-right:](survey.md){ .md-button .md-button--primary }
