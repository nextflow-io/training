# GitHub Codespaces

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

GitHub Codespaces és una plataforma basada en web que ens permet proporcionar un entorn preconfigurat per a la formació, recolzat per màquines virtuals al núvol.
La plataforma està operada per Github (que és propietat de Microsoft), i és accessible gratuïtament (amb quotes d'ús) per a qualsevol persona amb un compte de Github.

!!! warning "Advertència"

    Els comptes vinculats a organitzacions poden estar subjectes a certes restriccions addicionals.
    Si aquest és el vostre cas, potser haureu d'utilitzar un compte personal independent, o utilitzar una instal·lació local en el seu lloc.

## Crear un compte de GitHub

Podeu crear un compte gratuït de GitHub des de la [pàgina d'inici de GitHub](https://github.com/).

## Llançar el vostre GitHub Codespace

Un cop hàgiu iniciat sessió a GitHub, obriu aquest enllaç al vostre navegador per obrir l'entorn de formació de Nextflow: <https://codespaces.new/nextflow-io/training?quickstart=1&ref=master>

Alternativament, podeu fer clic al botó que es mostra a continuació, que es repeteix a cada curs de formació (normalment a la pàgina d'Orientació).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Se us hauria de presentar una pàgina on podeu crear un nou GitHub Codespace:

![Create a GitHub Codespace](img/codespaces_create.png)

### Configuració

Per a l'ús general, no hauríeu de necessitar configurar res.
Tret que s'especifiqui el contrari al curs que esteu començant, simplement podeu fer clic al botó principal per continuar.

No obstant això, és possible personalitzar l'entorn fent clic al botó "Change options".

??? info "Opcions de configuració"

    Si feu clic al botó "Change options", se us donarà l'opció de personalitzar el següent:

    #### Branch

    Això us permet seleccionar una versió diferent dels materials de formació.
    La branca `master` generalment conté correccions d'errors i materials que s'han desenvolupat i aprovat recentment però que encara no s'han publicat al lloc web.
    Altres branques contenen treballs en curs que poden no ser completament funcionals.

    #### Machine type

    Això us permet personalitzar la màquina virtual que utilitzareu per treballar amb la formació.

    Utilitzar una màquina amb més nuclis us permet aprofitar millor la capacitat de Nextflow per paral·lelitzar l'execució del workflow.
    No obstant això, consumirà la vostra assignació de quota gratuïta més ràpidament, per la qual cosa no recomanem canviar aquesta configuració tret que s'aconselli a les instruccions del curs que teniu previst fer.

    Consulteu 'Quotes de GitHub Codespaces' a continuació per a més detalls sobre les quotes.

### Temps d'inici

Obrir un nou entorn de GitHub Codespaces per primera vegada pot trigar diversos minuts, perquè el sistema ha de configurar la vostra màquina virtual, així que no us preocupeu si hi ha un temps d'espera.
No obstant això, no hauria de trigar més de cinc minuts.

## Navegar per la interfície de formació

Un cop s'hagi carregat el vostre GitHub Codespaces, hauríeu de veure alguna cosa similar al següent (que pot obrir-se en mode clar depenent de les preferències del vostre compte):

![GitHub Codespaces welcome](img/codespaces_welcome.png)

Aquesta és la interfície de l'IDE VSCode, una aplicació popular de desenvolupament de codi que recomanem utilitzar per al desenvolupament amb Nextflow.

- **L'editor principal** és on s'obriran el codi de Nextflow i altres fitxers de text. Aquí és on editareu el codi. Quan obriu el codespace, això us mostrarà una vista prèvia del fitxer `README.md`.
- **El terminal** sota l'editor principal us permet executar comandes. Aquí és on executareu totes les línies de comandes donades a les instruccions del curs.
- **La barra lateral** us permet personalitzar el vostre entorn i realitzar tasques bàsiques (copiar, enganxar, obrir fitxers, cercar, git, etc.). Per defecte està oberta a l'explorador de fitxers, que us permet navegar pel contingut del repositori. Fer clic en un fitxer a l'explorador l'obrirà dins de la finestra de l'editor principal.

Podeu ajustar les proporcions relatives dels panells de la finestra com vulgueu.

<!-- TODO (future) Link to development best practices side quest? -->

## Altres notes sobre l'ús de GitHub Codespaces

### Reprendre una sessió

Un cop hàgiu creat un entorn, podeu reprendre'l o reiniciar-lo fàcilment i continuar des d'on ho vau deixar.
El vostre entorn expirarà després de 30 minuts d'inactivitat i desarà els vostres canvis durant un màxim de 2 setmanes.

Podeu reobrir un entorn des de <https://github.com/codespaces/>.
Es llistaran els entorns anteriors.
Feu clic a una sessió per reprendre-la.

![List GitHub Codespace sessions](img/codespaces_list.png)

Si heu desat l'URL del vostre entorn anterior de GitHub Codespaces, simplement podeu obrir-lo al vostre navegador.
Alternativament, feu clic al mateix botó que vau utilitzar per crear-lo en primer lloc:

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Hauríeu de veure la sessió anterior, l'opció per defecte és reprendre-la:

![Resume a GitHub Codespace](img/codespaces_resume.png)

### Desar fitxers a la vostra màquina local

Per desar qualsevol fitxer des del panell de l'explorador, feu clic amb el botó dret al fitxer i seleccioneu `Download`.

### Gestionar les quotes de GitHub Codespaces

GitHub Codespaces us proporciona fins a 15 GB-mes d'emmagatzematge per mes, i 120 hores-nucli per mes.
Això equival a aproximadament 60 hores de temps d'execució de l'entorn per defecte utilitzant l'espai de treball estàndard (2 nuclis, 8 GB de RAM i 32 GB d'emmagatzematge).

Podeu crear-los amb més recursos (vegeu l'explicació anterior), però això consumirà el vostre ús gratuït més ràpidament i tindreu menys hores d'accés a aquest espai.
Per exemple, si seleccioneu una màquina de 4 nuclis en lloc dels 2 nuclis per defecte, la vostra quota s'esgotarà en la meitat del temps.

Opcionalment, podeu comprar accés a més recursos.

Per a més informació, consulteu la documentació de GitHub:
[About billing for GitHub Codespaces](https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces)
