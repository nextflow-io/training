# Resum del curs

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Felicitats per completar el curs de formació Use nf-core! 🎉

<!-- placeholder for video -->

## El vostre recorregut

Heu començat trobant i recuperant el pipeline `nf-core/demo`, i després heu après a executar-lo utilitzant el seu perfil de prova i a examinar les seves sortides.
A continuació, heu configurat la seva execució mitjançant paràmetres de pipeline i fitxers de configuració, i heu vist com els pipelines d'nf-core validen els paràmetres i les dades d'entrada.
Finalment, heu aplicat aquestes mateixes habilitats a `nf-core/rnaseq`, un pipeline a escala de producció, i heu après a sobreescriure les seves assignacions de recursos per defecte per adaptar-les al maquinari disponible.

### Què heu après

Ara sou capaços de trobar, recuperar, executar i configurar pipelines d'nf-core.

- Els pipelines d'nf-core es recuperen amb `nextflow pull` i segueixen una organització de codi estàndard.
- Cada pipeline d'nf-core inclou un perfil `test` per a una validació ràpida amb un conjunt de dades petit.
- Els paràmetres de pipeline (definits amb `--param_name` o `-params-file`) i la configuració (definida amb `-c`) tenen propòsits diferents: entrades i opcions d'anàlisi versus logística d'execució com l'assignació de recursos.
- Els pipelines d'nf-core validen els paràmetres i els fitxers d'entrada automàticament, detectant errors abans que es faci cap feina.
- Els recursos per defecte s'assignen mitjançant etiquetes (`process_low`, `process_medium`, `process_high`) definides a `conf/base.config`, que podeu sobreescriure amb un fitxer de configuració personalitzat.

### Habilitats adquirides

A través d'aquest curs pràctic, heu après a:

- Trobar un pipeline d'nf-core al lloc web nf-co.re i recuperar el seu codi font
- Executar un pipeline utilitzant el seu perfil de prova integrat i examinar les seves sortides
- Obtenir ajuda, definir paràmetres i entendre la validació de paràmetres i dades d'entrada
- Personalitzar l'assignació de recursos i els arguments de les eines mitjançant fitxers de configuració
- Recuperar i executar un pipeline a escala de producció, i sobreescriure les seves etiquetes de recursos per defecte

Ara disposeu dels coneixements fonamentals per començar a executar pipelines d'nf-core per a les vostres pròpies anàlisis.

## Passos següents per millorar les vostres habilitats

Aquí teniu les nostres principals recomanacions sobre què fer a continuació:

- Llanceu i monitoritzeu aquests pipelines a escala amb [Scale with Seqera](../seqera_scale/index.md)
- No us limiteu a executar pipelines d'nf-core, desenvolupeu-los! Apreneu les bones pràctiques d'nf-core amb [Build with nf-core](../nfcore_build/index.md)
- Sou nous a Nextflow? Comenceu amb [Nextflow Run](../nextflow_run/index.md)
- Apliqueu Nextflow a un cas d'ús d'anàlisi científica amb [Nextflow for Science](../nf4_science/index.md)
- Exploreu funcionalitats més avançades de Nextflow amb els [Side Quests](../side_quests/index.md)

## Obtenir ajuda

Per a recursos d'ajuda i suport de la comunitat, consulteu la [pàgina d'ajuda](../help.md).

## Enquesta de valoració

Abans de continuar, dediqueu un minut a completar l'enquesta del curs! Els vostres comentaris ens ajuden a millorar els materials de formació per a tothom.

[Feu l'enquesta :material-arrow-right:](survey.md){ .md-button .md-button--primary }
