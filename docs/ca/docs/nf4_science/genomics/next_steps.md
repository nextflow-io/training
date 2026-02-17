# Resum del curs

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Felicitats per haver completat el curs de formació de Nextflow per a Genòmica! 🎉

## El teu recorregut

Has començat executant eines de detecció de variants manualment al terminal per entendre la metodologia.
Després has construït un pipeline de Nextflow per a una sola mostra per automatitzar el procés, l'has escalat per gestionar múltiples mostres en paral·lel, i has afegit genotipat conjunt de múltiples mostres utilitzant operadors de canal.

### El que has construït

- Un pipeline de detecció de variants que pren fitxers BAM com a entrada i produeix VCF amb detecció conjunta com a sortida.
- Tres processos (`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER` i `GATK_JOINTGENOTYPING`) emmagatzemats en fitxers de mòdul separats.
- El pipeline paral·lelitza automàticament el processament de les mostres d'entrada utilitzant el paradigma de flux de dades de Nextflow.
- Els resultats es publiquen en un directori anomenat `results/`.

### Habilitats adquirides

A través d'aquest curs pràctic, has après com:

- Escriure un workflow lineal per aplicar detecció de variants a una sola mostra
- Gestionar fitxers accessoris com fitxers d'índex i recursos del genoma de referència adequadament
- Aprofitar el paradigma de flux de dades de Nextflow per paral·lelitzar la detecció de variants per mostra
- Implementar detecció conjunta de múltiples mostres utilitzant operadors de canal rellevants

Ara estàs preparat/da per començar a aplicar Nextflow a workflows d'anàlisi genòmica en el teu propi treball.

## Propers passos per desenvolupar les teves habilitats

Aquí tens les nostres principals recomanacions sobre què fer a continuació:

- Aplica Nextflow a altres casos d'ús d'anàlisi científica amb [Nextflow for Science](../index.md)
- Comença amb nf-core amb [Hello nf-core](../../hello_nf-core/index.md)
- Explora funcionalitats més avançades de Nextflow amb les [Side Quests](../../side_quests/index.md)

Finalment, et recomanem que donis una ullada a [**Seqera Platform**](https://seqera.io/), una plataforma basada en el núvol desenvolupada pels creadors de Nextflow que fa encara més fàcil llançar i gestionar els teus workflows, així com gestionar les teves dades i executar anàlisis de manera interactiva en qualsevol entorn.

## Obtenir ajuda

Per a recursos d'ajuda i suport de la comunitat, consulta la [pàgina d'Ajuda](../../help.md).

## Enquesta de valoració

Abans de continuar, si us plau, dedica un minut a completar l'enquesta del curs! La teva valoració ens ajuda a millorar els nostres materials de formació per a tothom.

[Fes l'enquesta :material-arrow-right:](survey.md){ .md-button .md-button--primary }
