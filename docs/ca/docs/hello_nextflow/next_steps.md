# Resum del curs

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Felicitats per haver completat el curs de formació Hello Nextflow! 🎉

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Consulteu la [llista de reproducció completa al canal de YouTube de Nextflow](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n).

:green_book: Podeu llegir la [transcripció del vídeo](./transcripts/07_next_steps.md) juntament amb el vídeo.
///

## El vostre recorregut

Vau començar amb un workflow molt bàsic que executava una comanda codificada de manera fixa.
Al llarg de sis parts, heu transformat aquest workflow bàsic en un pipeline modular de múltiples passos que exercita característiques clau de Nextflow, incloent canals, operadors, suport integrat per a contenidors i opcions de configuració.

### El que heu construït

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

A través d'aquest curs pràctic, heu après a:

- Descriure i utilitzar components bàsics de Nextflow suficients per construir un workflow senzill de múltiples passos
- Descriure conceptes del següent nivell com ara operadors i factories de canals
- Llançar un workflow de Nextflow localment
- Trobar i interpretar sortides (resultats) i fitxers de registre generats per Nextflow
- Resoldre problemes bàsics

Ara esteu equipats amb el coneixement fonamental per començar a desenvolupar els vostres propis pipelines amb Nextflow.

## Passos següents per desenvolupar les vostres habilitats

Aquí teniu les nostres 3 principals suggerències sobre què fer a continuació:

- Apliqueu Nextflow a un cas d'ús d'anàlisi científica amb [Nextflow for Science](../nf4_science/index.md)
- Inicieu-vos amb nf-core amb [Hello nf-core](../hello_nf-core/index.md)
- Exploreu característiques més avançades de Nextflow amb les [Side Quests](../side_quests/index.md)

Finalment, us recomanem que doneu una ullada a [**Seqera Platform**](https://seqera.io/), una plataforma basada en el núvol desenvolupada pels creadors de Nextflow que fa encara més fàcil llançar i gestionar els vostres workflows, així com gestionar les vostres dades i executar anàlisis de manera interactiva en qualsevol entorn.

## Enquesta de valoració

Abans de continuar, si us plau, dediqueu un minut a completar l'enquesta del curs! La vostra valoració ens ajuda a millorar els nostres materials de formació per a tothom.

[Feu l'enquesta :material-arrow-right:](survey.md){ .md-button .md-button--primary }
