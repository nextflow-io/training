# Part 1: Hello World - Transcripció del Vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importants"

    Aquesta pàgina mostra només la transcripció. Per a instruccions completes pas a pas, torneu al [material del curs](../01_hello_world.md).

    Els números de secció mostrats a la transcripció es proporcionen només amb finalitats indicatives i poden no incloure tots els números de secció dels materials.

## Benvinguda

Hola, i benvinguts de nou.

Ara esteu a la Part 1 del curs "Hello Nextflow" anomenada "Hello World". En aquest capítol, començarem a construir una comprensió dels conceptes més bàsics de Nextflow.

Així que esperem que ara estigueu configurats a Codespaces o en algun lloc equivalent amb VS Code en execució, i tingueu la vostra carpeta Hello Nextflow a l'espai de treball a l'Explorador amb tots aquests fitxers diferents aquí.

Començarem fent algunes coses molt bàsiques al terminal utilitzant Bash, i després veurem si podem fer les mateixes coses dins de Nextflow perquè tingueu una idea de com es veu la sintaxi.

## 0. Escalfament

Així que comencem de manera molt simple. Comencem amb "echo", per imprimir alguna cosa a un terminal. "Hello World". Premo enter i això va a un terminal. Hello World. Esperem que això no sigui una sorpresa per a ningú que estigui veient aquest curs.

D'acord, fem alguna cosa amb això. En lloc de només imprimir-ho al terminal, escrivim-ho a un fitxer. Premeré la tecla de cursor amunt del meu teclat, que recorre l'historial de Bash, així que em dóna la meva última comanda, i afegiré al final d'aquesta, un petit símbol de major que, que redirigeix la sortida d'aquesta comanda a un fitxer, i l'anomenaré output.txt.

Enter de nou, per executar aquesta comanda, res al terminal aquesta vegada, però podem veure al costat esquerre, el nou fitxer ha aparegut aquí, anomenat output.txt.

Podem veure això en un terminal amb alguna cosa com cat. Així que cat output.txt i efectivament diu "Hello World". També podem fer doble clic i s'obre a l'editor de codi a VS Code.

## 1.1. Examinar el codi

D'acord. Us vaig dir que era simple. Què segueix? Intentem agafar aquest procés i fer-ho de nou, però aquesta vegada, fem-ho dins de Nextflow.

Com he dit, tots els diferents capítols d'aquest curs comencen amb un script i aquest s'anomena Hello World. Així que trobaré Hello World. El previsualitza si faig un sol clic, faré doble clic per obrir-lo a l'editor aquí. I eliminaré ràpidament el terminal.

Ara aquest és un script molt simple, tan simple com es pot aconseguir. Només té 22 línies de llarg, i fa bàsicament el mateix. De fet, algunes d'aquestes coses haurien de semblar familiars. És el que acabem d'escriure. Podem veure la nostra comanda bash redirigint a un fitxer allà.

D'acord. Què més? També, en aquest fitxer, podem començar a veure alguns dels conceptes bàsics de Nextflow. Tenim un procés en vermell aquí i un workflow. Aquestes són les paraules clau especials i la terminologia especial a Nextflow.

## 1.1.1. La definició del procés

Diferents processos dins d'un workflow embolcallen diferents unitats lògiques del vostre workflow. Cada procés fa una cosa.

Quan l'executem, genera una tasca o múltiples tasques, que són els passos reals d'execució d'un pipeline. Tots els processos són després orquestrats dins d'un bloc workflow, que veiem a la part inferior, i en aquest cas, només executa aquest únic procés.

El nom del procés segueix aquesta paraula clau aquí, i això pot ser bàsicament qualsevol cosa. I després els continguts del procés estan dins d'aquestes claus.

Només hi ha realment un requisit per al procés, que és que inclogui algun tipus de bloc script o exec. Això està entre les triples cometes aquí, i aquest és el script bash que s'escriu al directori de treball quan executem el pipeline i és la cosa que realment s'executa al vostre ordinador o servidor.

Això és bash típicament, però també podeu posar un shebang diferent aquí a la part superior, i podria ser un script Python o un script R. No importa. El que sigui que estigui en aquest script s'executarà.

Hi ha una altra cosa que hem afegit a aquest procés aquí, que és la declaració d'output. Això diu a Nextflow que aquest procés espera un fitxer de sortida anomenat output.txt. Diu que és un path, així que s'ha de gestionar com un fitxer, no diguem, si això fos val, diria que és com una variable o valor.

Tingueu en compte que això no està creant aquest fitxer. No el genera realment. Això es fa pel script aquí a baix. Només està dient a Nextflow que esperi un fitxer de sortida amb aquest nom de fitxer.

## 1.1.2. La definició del workflow

D'acord. I després a la part inferior tenim un workflow aquí, i de nou, tenim una declaració. Aquest s'anomena Main. Aquest és l'equivalent del workflow d'un bloc script, si voleu. És la part del workflow que fa alguna cosa. I en aquest cas, estem dient, crida el procés anomenat sayHello.

Normalment, per descomptat, el vostre pipeline es veurà molt més complex que això. Probablement tindreu més d'un procés, i utilitzareu canals per orquestrar el flux de dades entre ells. Arribarem a això a les següents parts d'aquest curs, però per ara, això és suficient. Aquest és un pipeline vàlid, que hauria de funcionar.

Fins i tot puc fer clic a previsualitzar DAG aquí a VS Code. El DAG o DAG és una representació d'una estructura de flux de dades al pipeline, i podem veure-ho representat al costat com un diagrama mermaid. En aquest cas és molt simple. Hi ha una caixa, que és el workflow i un procés, que s'anomena sayHello, però això podria semblar més interessant a mesura que avancem.

## 1.2. Executar el workflow

D'acord, intentem executar aquest workflow i veure què passa.

Tornaré a obrir el terminal a la part inferior, netejaré la sortida, i escriuré Nextflow Run. I després només escriuré el nom de l'script, que és hello-world.nf. I premeré enter.

D'acord, té algunes coses estàndard a la part superior, que ens diu que Nextflow s'ha executat i quina versió estava executant-se i quin era el nom de l'script i tot això.

I realment la cosa important que estem buscant aquí és _aquí_, que és un resum de les diferents tasques que es van executar.

Si el vostre es veu així amb una petita marca verda, doncs molt bé. Acabeu d'executar el vostre primer pipeline. Fantàstic.

Ens diu aquí el nom del procés, que es va executar, que s'anomenava Say Hello, i ens va dir que es va executar una vegada i que va tenir èxit. Això s'actualitza a mesura que avanceu, així que quan esteu executant un pipeline més gran, veureu el progrés representat aquí. Però com que això és tan petit, s'executa bàsicament immediatament.

## 1.2.2. Trobar la sortida i els registres al directori work

Ara quan executeu un pipeline de Nextflow, cadascun d'aquests processos està unit, i cada procés, com he dit abans, pot generar tasques una o múltiples. Així que en aquest cas, teníem una sola tasca d'aquest procés. Només es va executar una vegada i això es va fer sota aquest _hash_ de tasca.

Nextflow no tracta amb els fitxers del vostre directori de treball directament, crea una carpeta especial anomenada work. I si faig "ls", veurem que ha aparegut aquí: _work_, i dins d'aquí hi ha subdirectoris per a cada tasca individual que s'executa. I això coincideix amb aquest hash. Així que podeu veure si vaig a "ls work/c4", i després està truncat, però comença 203, i aquest és el directori de treball, que va ser creat per aquest procés quan vam executar el pipeline. I podeu veure-ho al costat també.

Quan llisto aquests fitxers, podeu veure que el fitxer output.txt es va generar. Podeu veure-ho aquí també. I hi ha un munt de fitxers ocults, que no es mostren amb el meu "ls" regular.

Si faig clic a output.txt, efectivament, tenim la nostra sortida. Fantàstic. Així que el pipeline va funcionar.

Pot semblar molt codi repetitiu per executar el que era essencialment un script bash d'una línia, però tindrà més sentit a mesura que els nostres processos es facin més complicats. I aquest directori work amb Nextflow i aquests fitxers, que es creen és realment la columna vertebral del que fa Nextflow tan potent.

Cada tasca, cada element d'un pipeline està aïllat de qualsevol altra tasca. És reproduïble. No entren en conflicte entre si, i tot pot executar-se en paral·lel. És realment una manera agradable quan t'hi acostumes a causa d'aquest aïllament que pots entrar i veure exactament què va passar per a una sola tasca i depurar.

Fem una ullada ràpida a aquests altres fitxers al directori work. De dalt a baix, tenim un fitxer anomenat _.command.begin_. Això està buit. És només el que s'anomena un fitxer sentinella, creat per Nextflow dient, d'acord, estic començant la tasca. Res interessant allà.

Després hi ha _.command.error_, _.command.log_ i _.command.out_. Aquestes són totes sortides de la comanda bash o aquest script que es va executar. Això és error estàndard. Això és sortida estàndard, i això són les dues combinades tal com van sortir. Així que obteniu l'ordre lògic.

D'acord, aquests també estaven tots buits per a això, així que no gaire interessant, però les coses es tornen més interessants quan arribeu a _.command.run_.

Aquest és típicament un script molt llarg. I això és el que Nextflow realment executa. Si entreu aquí, començareu a veure tota la lògica interna de Nextflow i veure què està fent i com està executant el vostre procés. Això dependrà d'on estigueu executant, si estem executant localment o enviant-ho com un treball a SLURM, en aquest cas tindrem capçaleres SLURM a la part superior. Totes aquestes configuracions diferents.

Generalment, no necessiteu mai mirar aquest fitxer, però. És autogenerat per Nextflow i no hi ha res realment particularment únic del vostre pipeline, que hi sigui. Però això és realment el nucli del que s'està executant.

El següent és molt més interessant. _.command.sh_ és l'script generat, que va venir del vostre procés, i aquí podeu veure que Nextflow va afegir la capçalera Bash, i després va executar la nostra comanda, que estava al nostre bloc script.

I això és tot el que fa el fitxer _.command.run_ és que només executa aquest fitxer _.command.sh_.

Aquest és realment útil, que és el que normalment acabeu mirant més quan esteu intentant depurar alguna cosa i comprovar que la lògica del vostre pipeline de Nextflow està fent el que espereu que faci.

Finalment, tenim un fitxer anomenat _.exitcode_, i això només captura el codi de sortida d'una tasca, que en aquest cas va tenir èxit. Així que el codi de sortida era zero.

Si alguna cosa va malament, us quedeu sense memòria o alguna altra cosa i falla, llavors això és molt útil per entendre què va anar malament.

## 1.3. Executar el workflow de nou

Una cosa més a entendre sobre els directoris work és que si continuo executant aquest pipeline repetidament, així que si faig _"nextflow run hello-world.nf"_, farà exactament el mateix, però aquesta vegada tindrà un nou id de tasca. Podeu veure que aquest hash aquí és diferent, i ara si miro a work, hi ha dos directoris hash. I aquests són, de nou, separats l'un de l'altre.

Així que cada vegada que executeu un workflow de Nextflow, tret que utilitzeu el resume, que utilitza la memòria cau, tocarem més endavant, tornarà a executar aquests processos en nous directoris work, que estan separats l'un de l'altre. No tindreu cap col·lisió de noms de fitxers, no tindreu cap problema com aquest. Tot està aïllat i net.

I si entrem en aquest directori, podeu veure tots els mateixos fitxers i el mateix _output.txt_, que s'ha recreat des de zero.

## 2. Publicar sortides

D'acord, això és genial per a Nextflow per si mateix, mentre està executant el vostre pipeline perquè totes les coses estan separades l'una de l'altra i netes i es poden gestionar.

Però no és súper útil si sou una persona intentant explorar els vostres resultats. Realment no voleu estar excavant a través de milers i milers de directoris work diferents intentant trobar els vostres fitxers de resultats. I realment no se suposa que ho feu. Els directoris work no estan destinats a ser l'estat final d'on es creen els vostres fitxers.

Ho fem publicant els nostres fitxers.

## 2.1.1. Declarar la sortida del procés sayHello

Així que si torno al nostre script, treballarem al nostre bloc workflow aquí. Li direm quins fitxers esperar, quins fitxers ens importen, i després crearem un nou bloc a sota anomenat bloc output.

Aquesta és la nova sintaxi, que va venir amb l'analitzador de sintaxi i per defecte a la versió 26.04 de Nextflow. Així que si heu utilitzat Nextflow una mica abans, aquesta és una de les coses que és nova.

Així que tenim el bloc main, i després diré publish i diré a Nextflow què esperar de la publicació. L'anomenarem _first_output_, i l'anomenarem, _sayHello.out_.

Accidentalment vaig fer un error tipogràfic allà, però aquesta és una bona oportunitat per assenyalar també algunes de les característiques de l'extensió Nextflow VS Code. Podeu veure que immediatament em va donar una petita línia vermella ondulada sota això dient que alguna cosa està malament. I si passo el cursor per sobre, em dirà que aquesta variable no està definida. No sé què és.

És bastant obvi en aquest cas, vaig fer un error tipogràfic. Volia escriure, sayHello, i després la línia ondulada desapareix.

Ara és porpra. L'analitzador de sintaxi de Nextflow sap que això és un procés i quan passo el cursor per sobre, em dóna una representació reduïda de com es veu aquest procés. Així que puc veure molt ràpidament d'una ullada que no pren cap entrada i ens dóna aquesta sortida. Així que treballar a VS Code amb aquesta extensió us dóna molta informació contextual mentre esteu escrivint codi.

Tingueu en compte que podem referir-nos a la sortida d'aquest procés amb la sintaxi _.out_. I de moment podem anomenar això com vulguem, és només un nom de variable arbitrari.

## 2.1.2. Afegir un bloc output: a l'script

On es torna important és quan fem el nostre nou bloc aquí, i això està per sota del bloc workflow ara, ja no estem dins del workflow. Claus de nou. I aquí és on només diem a Nextflow on posar tots els fitxers, que són creats pel workflow.

Ara agafaré aquest nom de variable, que vaig crear aquí, i el posaré allà i posaré algunes claus per a això. I diré a Nextflow que utilitzi un path. Ups. Path, entre cometes. I utilitzaré punt. Això només diu a Nextflow que posi el fitxer a l'arrel del directori results. Així que no cap subdirectori ni res.

Intentem executar el nostre workflow de nou. Si faig _"nextflow run hello-world.nf"_, llavors esperem que hauria de semblar bàsicament exactament el mateix. Res ha canviat realment amb Nextflow aquí. Està executant les mateixes coses. Només les està fent en directoris work de nou.

Però ara si faig _"ls results/"_, veureu que hi ha un nou directori aquí que s'ha creat anomenat results, que és el directori base per defecte per a la publicació del workflow. I allà hi ha un fitxer anomenat _output.txt_.

Si faig _"ls -l results"_, veureu que això està realment enllaçat de manera suau al directori work. Així que això no és un fitxer real, està enllaçat al directori work i ha recollit tots els fitxers allà per a nosaltres.

## 2.2. Establir una ubicació personalitzada

"Results" és el nom per defecte per a aquest path. Si executo el workflow de nou, i aquesta vegada faig _guió_ guió únic, això és perquè és una opció bàsica de Nextflow. _"-output-dir **my**results"._ També podria fer només _"-o"_ per abreujar. Llavors establirà un directori base diferent per a on s'emmagatzemen els fitxers i una vegada més, aquí a _myresults/_, ara tenim un _output.txt_.

Això és genial, però probablement no volem tots els fitxers només a l'arrel. Volem una mica d'organització, així que també podem crear un subdirectori aquí anomenat com vulguem. Diguem _"path 'hello_world'"_, i només executo això de nou. _"nextflow run hello-world.nf"_. Hauria d'anar al directori results a un subdirectori i efectivament, ara sota results aquí a la part superior tenim _hello_world/_ i tenim _output.txt_.

Cosa important a notar, l'antic fitxer _output.txt_ encara hi és. El directori results no s'esborra quan feu això. Només es copien nous fitxers allà. Sobreescriuran fitxers que ja hi són si tenen el mateix nom de fitxer, però no netejaran els antics. Així que heu de tenir una mica de cura sobre quan torneu a executar pipelines. Si no voleu que estiguin a sobre dels fitxers que ja hi són. Assegureu-vos d'utilitzar un directori buit.

## 2.3. Establir el mode de publicació a copy

D'acord. He esmentat que aquests fitxers són enllaços suaus, així que si faig _"ls -l results/hello_world/"_, podeu veure que està enllaçant de manera suau al directori work. Això és generalment una cosa bona si esteu treballant en alguna cosa com HPC, i aquests són fitxers realment enormes i no voleu duplicar-los, perquè significa que els fitxers només s'emmagatzemen una vegada al sistema de fitxers.

No obstant això, sí que significa que si elimineu el directori work: si faig _"rm -r work"_ i netejo tots aquests fitxers intermedis que es van crear. Ara, si intento llegir aquest fitxer _"results/hello_world/"_. Estarà apuntant com un enllaç suau a un fitxer que ja no existeix i les dades s'han perdut per sempre i són irrecuperables, cosa que potser no és genial.

Així que generalment nosaltres, dic que és una bona pràctica copiar els fitxers en lloc d'enllaçar-los de manera suau si podeu, perquè és més segur. Només tingueu en compte que utilitzarà el doble d'espai de disc tret que elimineu aquests directoris work.

Per fer això amb el bloc output, aniré al primer output aquí. Vaig establir el path abans i ara establiré el mode i podeu veure mentre escric, l'extensió VS Code està, suggerint coses que sap que és una directiva d'output aquí. I diré copy. Premo desar.

Tornem a executar el workflow. Crearà els fitxers de nou, nou directori work.

Ara, si vaig a _"ls -l results/hello_world/"_ podeu veure que això és un fitxer real i ja no és un enllaç suau, i Nextflow va copiar això. Bo de saber. Així que path i mode són coses que us trobareu escrivint força.

Ara, per descomptat, això és molt simple. Farem això més complex i potent a mesura que avancem, i veureu com fer aquestes coses dinàmiques i no massa verboses.

## 2.4. Nota sobre les directives publishDir a nivell de procés

Ara, vaig dir quan vam començar amb això, que aquesta és una forma de sintaxi bastant nova. Només està disponible a les últimes versions de Nextflow mentre gravo això, i s'anomena Workflow Outputs.

Si utilitzeu això, és genial. Desbloqueja moltes altres característiques interessants dins de Nextflow, com ara, Nextflow Lineage per ajudar a rastrejar l'herència d'aquests fitxers a mesura que es creen, i aviat serà el valor per defecte a 26.04. I en una data posterior al futur, aquesta serà l'única manera d'escriure els vostres workflows.

No obstant això, com estem en aquesta fase de transició ara mateix, és molt possible que vegeu pipelines en estat salvatge, que utilitzeu alguna cosa anomenada publishDir, que és l'antiga manera de fer-ho, i això es defineix no al nivell de workflow i output, sinó que això es defineix a nivell de procés.

I aquesta declaració diu bàsicament el mateix. Diu, publica els fitxers de resultats en un directori anomenat results, i utilitza un mode copy. Així que podeu veure que la sintaxi és molt similar. Però quan esteu escrivint nous pipelines ara, intenteu no utilitzar aquesta directiva publishDir, fins i tot si la veieu, en resultats d'IA o en documentació o altres pipelines, perquè aquesta és l'antiga manera de fer-ho.

El 2026 tots hauríem d'estar utilitzant workflow outputs.

Tot això està documentat, si esteu fent això i heu utilitzat Nextflow abans, podeu anar als documents de Nextflow aquí, nextflow.io/docs/. I si faig scroll cap avall a tutorials, hi ha un tutorial anomenat _Migrating to Workflow Outputs_.

És realment bo. Passa per tota la sintaxi, com és equivalent a l'antiga sintaxi, per què la vam canviar, i, té una línia de temps i tot. I passa per tots els diferents escenaris amb càrregues i càrregues d'exemples. Així que podeu convertir fàcilment el codi Nextflow existent a la nova sintaxi.

## 3.1. Canviar el procés sayHello per esperar una entrada variable

D'acord, així que tenim el nostre script simple, que està executant un procés, creant un fitxer, dient a Nextflow que és una sortida, i després estem dient a Nextflow on desar aquest fitxer. Això és un bon començament.

Però seria més interessant si no estigués tot codificat de manera fixa. Així que a continuació, pensem en com dir a Nextflow que aquest procés pot prendre una entrada variable, que és alguna cosa que podem controlar en temps d'execució quan llancem un workflow.

Necessitem fer algunes coses diferents per fer que això passi.

Primer, necessitem dir a aquest procés que pot acceptar una variable d'entrada i escrivim _input_ aquí com un nou bloc de declaració. I l'anomenarem _"val greeting"_.

La part val és l'equivalent d'un path aquí a baix. Diu a Nextflow que això és una variable, com una cadena en aquest cas. I si passeu el cursor per sobre de nou, us diu de l'extensió del que això significa.

A continuació direm a Nextflow què fer amb això. No és suficient només dir que hi ha una variable. Heu de dir a l'script com utilitzar aquesta variable. I així que eliminaré aquesta cadena codificada de manera fixa aquí, i posaré una variable.

Ho faré ràpidament sense claus només per mostrar-vos que això està permès, i aquesta és l'antiga manera d'estil de fer-ho. Però ara amb la nova sintaxi, realment recomanem posar-ho dins de claus així, i deixa molt clar que això està sent interpolat per Nextflow aquí.

Genial. Així que _"input greeting"_ va a _$\{greeting\}._ L'última cosa és que necessitem dir a Nextflow al nivell de workflow que aquest procés ara pren una entrada. I per fer això, bàsicament li donarem una variable.

## 3.2. Configurar un paràmetre de línia de comandes per capturar l'entrada de l'usuari

Podríem codificar-ho de manera fixa de nou, com Hello World, i això funcionaria bé, però òbviament no ens dóna realment cap avantatge. Volíem poder configurar això en temps d'execució, així que volem poder fer-ho a la CLI, quan llancem Nextflow.

I la manera de fer-ho és un concepte especial de Nextflow anomenat _params_. L'anomenarem _params.input_.

El que això fa és que exposa aquesta variable input a la CLI i és on utilitzem un guió doble quan llancem Nextflow.

Puc anomenar això com vulgui, puc anomenar-ho _hello, greeting_. No importa. El que faci allà s'exposarà com una opció CLI quan llancem un pipeline. I aquest és un veritable truc de màgia per Nextflow perquè significa que podeu construir el vostre script de workflow molt ràpidament amb aquests paràmetres, i essencialment esteu construint una CLI personalitzada per al vostre pipeline, fent-lo realment fàcil de personalitzar diferents opcions sobre la marxa quan llancem.

Així que. Provem-ho. Tornem al nostre terminal. Tenim la nostra comanda _"nextflow run"_ aquí. I ara faré _"--input"_, que coincideix amb el _"params.input"_ que vam veure abans. Crec que als documents està en francès. A la Geraldine li agrada parlar francès. Ho faré en suec perquè visc a Suècia. així que diré, "_Hej Världen_" i premo enter.

Puc utilitzar cometes simples o dobles, només afecta com Bash ho interpreta.

Executa el pipeline de Nextflow exactament de la mateixa manera. Podeu veure el directori de treball i tot és el mateix. Però ara si vaig a _"results/hello_world/output"_. Podem veure el nostre bonic suec aquí en lloc d'això.

Així que hem passat dinàmicament una entrada d'una CLI a un paràmetre. Hem passat això com una entrada a un procés i el procés ha interpretat això i l'ha posat en un bloc script, que després ha canviat dinàmicament la sortida d'aquest resultat de script. Bastant interessant.

Lògica bastant complexa amb molt poca sintaxi aquí. I esperem que pugueu veure com això ara comença a escalar. I així és com realment construïm la lògica i la personalització dels nostres pipelines a l'script de Nextflow.

## 3.4. Utilitzar valors per defecte per a paràmetres de línia de comandes

D'acord, això és genial. El problema però ara és, cada vegada que executo aquest pipeline, necessito fer guió, input perquè s'executi.

Si intento executar sense aquest paràmetre, ara Nextflow llançarà un error dient que necessitava aquest paràmetre i no estava establert. i així que no sabia què fer.

Això és una cosa nova interessant, per cert. En el passat, Nextflow simplement s'hauria executat amb una cadena buida, i hauríeu tingut tot tipus d'errors estranys, que haurien estat difícils d'entendre. Però en el nou analitzador de sintaxi de Nextflow, és una mica més acurat i us ho diu immediatament.

Així que no sempre volem especificar cada opció individual. És una bona pràctica especificar valors per defecte sensats. Així que com ho fem al nostre script?

Notareu que quan vam escriure això, només vam posar _params.input_ directament on l'estem utilitzant. Així que la solució òbvia és que definim un valor per defecte, i ho fem a la part superior de l'script aquí en un bloc params especial al workflow. Això està al script del workflow aquí.

De nou, alguna sintaxi nova aquí, així que pareu atenció. Això és realment interessant. Tenim el nom del paràmetre, que s'esperarà aquí.

I després després d'aquest caràcter de dos punts, estem definint un tipus de la variable. No heu de fer això, podeu deixar-ho en blanc, però és realment agradable. Diu a Nextflow que estem esperant una cadena i tractar-la com a tal.

Si volem un número en lloc d'això, per exemple, podríem escriure float, i això diria que volem un número de punt flotant. I si intentem executar amb això, llavors llançarà un error. Si li donem una cadena, que no és un float. I també el passarà com a tal. Com si fem string, llavors sap que és una cadena. I fins i tot si té zeros inicials i és tot numèric, encara el passarà com una cadena real.

Així que aquesta seguretat de tipus és una característica molt nova de Nextflow, però realment potent per fer el vostre codi més segur d'escriure i executar.

Després d'això tenim un símbol igual i després el valor per defecte aquí. Nextflow es va escriure a Barcelona originalment, així que sembla apropiat que tinguem una mica d'espanyol aquí, _"Holà mundo!"_ com a valor per defecte.

D'acord, desaré aquest script, tornaré, executaré l'script de nou sense _--input_. I aquesta vegada hauria d'executar-se i crearà el nostre nou fitxer a _results_. I en aquest fitxer ara diu _"Holà mundo!"_.

Això és només un valor per defecte però, així que no significa que no puguem fer encara el mateix que abans. Si torno i trobo el meu antic script aquí, _"Hej Världen"_, perquè faig _--input_ a la línia de comandes, això sobreescriurà aquest valor per defecte i utilitzarà això de nou al fitxer output.txt.

Així que això a l'script és només el valor per defecte que estic establint.

A mesura que construïm el nostre workflow per ser més complex i incloure més paràmetres, aquest bloc params a la part superior de l'script començarà a recollir-los tots en un sol lloc.

I acabeu amb aquesta simetria bastant agradable al vostre script, on efectivament teniu totes les vostres entrades de workflow aquí i les vostres sortides de workflow a la part inferior. I és molt clar quina és la interfície del vostre workflow amb el món exterior. Així que podeu agafar un nou pipeline molt ràpidament amb la nova sintaxi i entendre com utilitzar-lo.

Una última cosa interessant. No hem d'establir un valor per defecte amb això. Si fem params input però no establim un valor per defecte, llavors diu a Nextflow que aquest paràmetre és obligatori, i de nou, el pipeline fallarà a executar-se sense ell, però us donarà un missatge d'error més útil en lloc d'alguna cosa sobre que és nul.

Així que diu que estem esperant que la seva entrada sigui obligatòria, però no es va especificar a la línia de comandes. Molt bé.

D'acord, així que esperem que ara estigui clar sobre com configurar el vostre pipeline de Nextflow amb entrades variables i paràmetres, com establir el valor per defecte, establir els tipus, podria ser un booleà true false flag o un enter o diferents tipus aquí. Com passar-los al vostre workflow, on passa, i després també interpola al vostre procés. I després també sabeu com personalitzar-los a la línia de comandes quan llancem Nextflow. Això està començant a semblar més interessant que la nostra simple comanda bash.

## 4. Gestionar execucions de workflow

D'acord. Què segueix? Per a la part final d'aquest capítol, parlarem una mica sobre com gestionar totes les diferents execucions de workflow. Si mireu a la meva barra lateral aquí i l'Explorador sota work, veureu que he executat un munt de pipelines diferents i aquests directoris work s'estan fent força llargs, hi ha molts d'ells.

I l'altra cosa és, com he dit abans, cada vegada que torno a executar aquest pipeline, està creant un nou conjunt de directoris work, i està tornant a executar tots els processos des de zero, que és una cosa bona. Això és comportament intencionat. És reproduïble i està regenerant tot de nou. Però òbviament, si esteu executant processos de molt llarga durada, és molest haver de començar sempre el vostre pipeline des del principi si es va bloquejar a la meitat, o si canvieu alguna cosa al final del pipeline.

## 4.1. Tornar a llançar un workflow amb -resume

Per sort, Nextflow és realment bo sabent què s'ha executat prèviament i què està disponible, i reutilitzar aquests resultats antics és molt simple. Només afegim una nova bandera al final de la comanda _"-resume"_.

Ara, tingueu en compte que hi ha dos guions a input perquè això és el paràmetre. Només hi ha un guió a resume perquè això és una opció bàsica de Nextflow.

Això fa ensopegar la gent tot el temps, fins i tot si heu estat utilitzant Nextflow durant molt de temps. Així que recordeu sempre un o dos guions. Depèn de si és una opció bàsica de Nextflow.

D'acord, així que ara faig _-resume_ i executo exactament el mateix workflow de nou. I aquesta vegada hauria de semblar bastant exactament el mateix amb una diferència clau.

A la sortida aquí, podeu veure que els resultats es van emmagatzemar a la memòria cau. I de fet, aquest hash de tasca aquí és exactament el mateix que l'execució anterior, i només ha reutilitzat aquest directori work en la seva totalitat. Les entrades i les sortides i l'script no es van modificar. I així que només pren aquest fitxer d'això i si hi ha passos posteriors al procés, els passaria al següent pas del pipeline.

Així que encara està executant tot el pipeline de principi a fi, però està utilitzant resultats emmagatzemats a la memòria cau per a cadascuna d'aquestes tasques, on pot.

Ara, quan feu _-resume_, només reprèn l'última execució de pipeline al vostre directori de treball, el que sigui que fos. Però realment podeu reprendre des de qualsevol execució anterior que hàgiu fet allà. I n'hem fet força ara.

## 4.2. Inspeccionar el registre d'execucions passades

Per mirar-les totes, podem fer _"nextflow log"_ en lloc de _"nextflow run"_, i això ens donarà una sortida agradable mostrant totes aquestes diferents.. Necessito fer la meva pantalla una mica més petita perquè puguem veure-ho, totes aquestes diferents execucions quan les vam fer, l'id de sessió, la comanda i tot.

I podem mirar aquí i podem agafar el nom d'execució de qualsevol d'aquestes i després reprendre una d'aquestes específiques. Així que puc tornar i puc reprendre aquella anomenada _hungry_ekeblad_. I només poso això després del _resume_.

Si teniu curiositat, per cert, tots aquests adjectius i noms de científics estan al codi font de Nextflow. És una manera realment bona d'aconseguir la vostra primera sol·licitud de pull a Nextflow anant i trobant-ho i afegint el vostre científic favorit.

I de totes maneres, així que vaig fer això i va tornar i va mirar els resultats emmagatzemats a la memòria cau d'aquesta execució de workflow, es va adonar que encara podia reutilitzar-los, i ho va fer. Així que vaig obtenir els resultats emmagatzemats a la memòria cau de nou.

## 4.3. Eliminar directoris work més antics

Això és genial. Què passa si vull netejar aquests directoris work? Hi ha càrregues d'ells aquí. Hi ha càrregues de fitxers. Potser sé amb certesa que vull reprendre des de les últimes parelles d'execucions de pipeline, però no m'importen totes les anteriors.

Llavors puc triar-ne una aquí i puc utilitzar una altra comanda de Nextflow, que és _"nextflow clean"_, i puc fer _"nextflow clean"_, faré _"-before"_, i el nom d'execució particular, que en aquest cas era _reverent_pike_ i faré _"-n"_, que diu a Nextflow que només faci una execució de prova. Així que només em diu el que eliminarà. Sense fer res realment, així que eliminaria aquests directoris work.

Això sembla sensible. Així que faré la mateixa comanda de nou, però en lloc de _"-n"_ faré _"-f"_ per fer realment la neteja. I aquesta vegada ha eliminat realment tots aquests directoris. I si entro i miro els directoris work, ara està sembla molt més lleuger. Fantàstic.

Així que així és com netejar tots els vostres directoris work locals d'una manera bastant segura sense destruir completament la memòria cau. Així que encara podeu reprendre si voleu.

Si mai oblideu quines són aquestes banderes per a cada comanda de Nextflow podeu fer _"nextflow help"_, i després el nom de la comanda. Així que si faig _"nextflow help clean"_, podeu veure totes les diferents opcions: _-after, -before, -but_, totes diferents maneres de configurar aquest comportament de neteja. Bastant interessant.

## Conclusió

D'acord, això és el final de la part 1 de Hello Nextflow. És un començament força intens del curs, però esperem que ara tingueu una comprensió bastant bona de com es veu un script de Nextflow; amb diferents parts clau, els processos, els workflows, les sortides, i els paràmetres. Sabeu com configurar-los amb sobreescriptures bàsiques des de la línia de comandes, com fer un bloc d'entrada dinàmic amb un script dinàmic i sabeu com gestionar totes les vostres execucions de càrrega de treball: veient el que ja heu executat, reprenent, netejant. Hi ha moltes coses. Heu arribat molt lluny. Així que si voleu fer una pausa i fer una petita passejada i una tassa de te, ara és probablement un bon moment. Us ho heu guanyat.

A partir d'aquí, bàsicament estem construint sobre aquesta base. Com podem fer això més complex, més potent? Com podem fer-ho més flexible? Fer les coses que volem fer la nostra anàlisi a escala.

## Qüestionari

Ara si feu scroll cap avall a la part 1, hello world, a la pàgina web veureu un petit qüestionari i això és alguna cosa nova que hem fet per a aquesta versió de la formació de Nextflow. I podeu passar i posar-vos a prova per comprovar que heu entès tot el material que hem fet en aquest capítol.

Això no s'envia a nosaltres ni res, només s'emmagatzema al vostre navegador. Així que no sabem quines són les vostres respostes, però és només una petita autocomprovació per assegurar-vos que no us heu perdut res o no heu entès malament res. I podeu provar-ho tantes vegades com vulgueu.

Si sou com jo, potser voleu quedar-vos al terminal a la vostra instància de VS Code, en aquest cas podeu escriure la comanda _quiz_ i després només dir-li quin capítol esteu. Així que fem _"Hello World"_, i després podeu fer exactament les mateixes preguntes del qüestionari, que estan al navegador web, però només al vostre terminal.

Genial. D'acord. Espero que gaudiu d'això. Divertiu-vos una mica i, us veurem al següent capítol en només un minut per parlar tot sobre els canals de Nextflow.
