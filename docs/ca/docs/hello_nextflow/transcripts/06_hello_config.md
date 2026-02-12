# Part 6: Hello Config - Transcripció del vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importants"

    Aquesta pàgina mostra només la transcripció. Per a instruccions completes pas a pas, torneu al [material del curs](../06_hello_config.md).

    Els números de secció mostrats a la transcripció es proporcionen només amb finalitats indicatives i poden no incloure tots els números de secció dels materials.

## Benvinguda

Hola, i benvinguts de nou a la Part Sis de Hello Nextflow. Aquesta secció tracta sobre configuracions, i és l'última part d'aquest curs.

Nextflow és particularment bo en dues coses: reproducibilitat i portabilitat. Les configuracions són on veiem la segona d'aquestes brillar realment. La capacitat de configurar un pipeline de Nextflow per executar-se de diferents maneres i funcionar en diferents sistemes, sense haver d'editar el codi subjacent del pipeline.

Aquest superpoder permet que els pipelines de Nextflow siguin reutilitzats per altres persones en diferents llocs, o a través de diferents infraestructures a les quals podeu tenir accés vosaltres mateixos.

Això significa que podeu desenvolupar codi de pipeline al vostre portàtil, pujar-lo al núvol, executar-lo al vostre HPC, i és el mateix codi de pipeline i s'executa a tot arreu.

En aquesta secció, passarem per alguns temes. Començarem amb com Nextflow gestiona els fitxers de configuració, d'on els carrega, i com els escriviu i com els estructureu, i aquesta separació entre el pipeline mateix i el que hauria d'anar en un fitxer de configuració.

Després passarem a alguns casos d'ús comuns com canviar on s'emmagatzemen els fitxers de sortida, i també com aconseguir que el pipeline funcioni en diferents infraestructures, tant utilitzant diferents tipus d'empaquetament de programari com enviant tasques a diferents infraestructures.

## Jerarquies de fitxers de configuració

D'acord, comencem. Quan es tracta de carregar fitxers de configuració, Nextflow pot obtenir-los de molts llocs diferents, que és una cosa bona i també pot ser una cosa lleugerament arriscada perquè de vegades pot ser una mica difícil saber d'on està obtenint un fitxer de configuració i en quin ordre carrega les coses.

Així que recomano molt que feu clic en aquest enllaç aquí, que ens porta a la documentació de Nextflow. I en aquesta pàgina de configuració, llista els llocs clau des d'on es carrega la configuració, i el més important, l'ordre en què es carreguen aquestes coses.

Així que podeu veure, podeu posar un fitxer de configuració al vostre directori home de Nextflow, que típicament és ".nextflow" al vostre directori home. I aquest fitxer sempre serà carregat per cada execució de Nextflow al vostre sistema.

El següent lloc a mirar és un fitxer a l'arrel del vostre repositori o directori de pipeline anomenat "nextflow.config".

Després d'això, un altre fitxer anomenat "nextflow.config", però aquesta vegada al directori des d'on esteu llançant Nextflow: el directori de llançament.

Finalment, podeu proporcionar camins de fitxers de configuració a la línia de comandes amb un argument "-c", i podeu fer-ho múltiples vegades. I s'apliquen en l'ordre que els especifiqueu.

Podeu proporcionar fitxers de configuració en tots aquests llocs si voleu, i es carregaran iterativament, cada un sobreescrivint l'anterior només en els àmbits de configuració on entrin en conflicte.

Aquest és un sistema realment potent perquè significa que podeu establir valors per defecte sensats i després anar sent gradualment més i més específics a mesura que us centreu en aquesta configuració.

## 0. Escalfament: Executar hello-config.nf

D'acord, tanquem això i saltem als nostres Codespaces i comencem. Com abans, he netejat aquí, he eliminat els meus directoris de resultats anteriors, el meu Nextflow i els meus directoris de treball i així successivament. No us preocupeu si encara teniu aquests fitxers per aquí. És només perquè estic molt ampliat i així les coses es tornen desordenades molt ràpidament altrament.

Treballarem amb hello-config.nf, l'últim fitxer del nostre directori, i això hauria de seguir des d'on vam deixar-ho a la secció anterior.

Així que tenim els nostres quatre processos diferents, que s'estan incloent des de fitxers de mòdul. Tenim els nostres paràmetres de pipeline, el nostre bloc de workflow on estem cridant els diferents processos i unint els canals, publicant els canals de sortida, i després el bloc de sortida a la part inferior on definim on s'haurien d'emmagatzemar aquests fitxers i com s'haurien de copiar.

També ja tenim un fitxer "nextflow.config" del capítol anterior, on habilitem Docker, i avui anirem construint sobre aquest fitxer.

Com abans, hem canviat el camí de sortida en aquest script principal a hello config, només perquè no entri en conflicte amb resultats anteriors que heu generat.

D'acord, comprovem ràpidament que tot encara funciona com esperem. Obro un terminal i faig nextflow run hello-config.nf. Nextflow es carrega. Hauria d'executar els nostres quatre processos diferents. Generar alguna bonica art ascii utilitzant cowpy i després desar els nostres resultats als nostres fitxers de resultats en aquest directori.

Puc fer una ullada ràpida aquí només per assegurar-me que aquests fitxers es veuen com esperem, i efectivament, aquí està el nostre gall dindi gegant. Genial.

## 1.1. Moure valors per defecte a nextflow.config

Ara el primer que farem és començar a moure algunes coses del nostre script al nostre fitxer de configuració.

I el que ens importa són principalment els paràmetres en aquesta etapa. Volem portar els valors per defecte al fitxer de configuració, perquè sigui més clar quins són els valors per defecte i sigui més fàcil per a la gent sobreescriure'ls.

Agafaré aquest bloc de params aquí del script i el posaré al fitxer de configuració. I hem de tenir una mica de cura aquí, perquè ara mateix la sintaxi és lleugerament diferent entre configuració i scripts. El fitxer de configuració no pot acceptar declaracions de tipus perquè no estem realment definint aquests params, només els estem referenciant. Així que eliminaré aquests.

Però altrament és molt similar. Tenim un bloc de params i després tenim els nostres diferents paràmetres d'entrada, paràmetre batch, paràmetre character.

Ara puc tornar al meu script i ja no necessito definir aquests valors per defecte perquè aquests valors ara estan al meu fitxer de configuració de Nextflow.

No obstant això, deixo els noms dels paràmetres i els seus tipus, perquè Nextflow conegui aquesta informació i encara pugui fer tota la seguretat de tipus i tot.

D'acord. Desem aquests fitxers i comprovem ràpidament que tot encara funciona igual que abans. No hi hauria d'haver cap canvi aquí. Hem mantingut els valors iguals. Només hem mogut on s'han definit.

Genial.

## 1.2. Utilitzar un fitxer de configuració específic per a l'execució

Ara, fins ara hem estat llançant Nextflow des del mateix directori on tenim el nostre script de pipeline. Així que el nostre directori de llançament i el nostre directori de pipeline són una mica la mateixa cosa.

Per mostrar com podem tenir diferents fitxers de configuració amb diferents directoris de llançament, crearem un nou subdirectori ara.

Així que diré mkdir, i l'anomenarem tux-run.

I després faré cd, canviar de directori a tux-run. I noteu que ara estem en el nostre directori de treball ja no està en el mateix directori que els scripts de pipeline.

D'acord, creem un nou fitxer "nextflow.config". Així que touch nextflow.config, i obrim-lo a VS Code. També podeu veure a la barra lateral aquí que ara estem en aquest subdirectori.

Ara podem agafar el mateix bloc de params que teníem al nextflow.config de nivell superior, copiar-lo aquí i ara podem canviar aquests valors.

Primer, les dades ara són un camí relatiu diferent perquè estem en un subdirectori, així que hem d'actualitzar això. I després canviarem batch a experiment, i canviarem el character de Turkey a tux.

Ara cliquem desar aquí, i provem-ho. Igual que amb data, ara he de dir ../ per arribar al script. Així que és hello-config. I premo enter.

El codi del pipeline no ha canviat gens, però ara tindrem dos conjunts de configuració carregant-se, i el fitxer de configuració del directori de llançament hauria de sobreescriure els valors per defecte, que es van establir al nextflow.config del pipeline, i hauríem d'obtenir diferents conjunts de resultats.

Efectivament, dins del nostre directori aquí, dins de tux-run, podeu veure que tenim un directori punt Nextflow i un directori de treball i això és perquè aquests es creen sempre al vostre directori de llançament. Així que aquests són diferents dels de treball i resultats que teníem d'execucions anteriors.

Ara, si miro a results, podem veure el nostre collected i aquí està el nostre petit personatge tux. Així que podeu veure que aquests paràmetres es van interpretar correctament.

## 1.3. Utilitzar un fitxer de paràmetres

D'acord. Abans quan estava parlant dels diferents fitxers de configuració que es podien carregar, vaig ometre un altre lloc d'on podem obtenir configuració.

Podeu obtenir-la d'una línia de comandes com hem vist amb noms de paràmetres amb doble guió, però també podem proporcionar un fitxer YAML o JSON, només de params.

El fitxer de configuració pot tenir tots els diferents tipus d'àmbits, però aquests fitxers són només paràmetres, i és una manera agradable i amigable per a l'usuari de proporcionar molts paràmetres alhora, i potser una manera una mica més reproduïble perquè els escriviu a un fitxer, així que és fàcil obtenir-los en una etapa posterior.

Així que tornem al nostre terminal i just abans que ho oblidem, assegurem-nos que pugem un directori, així que ja no estic al subdirectori, i miraré el fitxer YAML que tenim aquí anomenat test-params.yaml.

Així que si només faig code test-params.yaml, podeu veure que això és només un fitxer YAML regular. Res especial. Amb les claus sent els nostres noms de paràmetres, amb el format YAML així que dos punts aquí, i després un valor.

Noteu que això no és codi Nextflow, així que no podem posar coses com variables aquí. Aquests són només valors estàtics.

També perquè JSON realment s'analitza com YAML, també podem tenir un fitxer test-params.json, que es veu molt similar. És només un format de dades diferent.

Així que tenim dos fitxers de prova diferents aquí i tenim variables lleugerament diferents.

D'acord, així que com els donem a Nextflow? És molt simple. Fem nextflow run hello-config, com abans. I en lloc de "-c" per a fitxer de configuració, o carregar aquests noms de fitxer per defecte, fem -params-file. Guió únic perquè és una opció central de Nextflow.

I després passem el camí per a aquest fitxer. Així que faré "-params-file test-params.yaml", i veurem si aquests es carreguen correctament.

D'acord. S'ha executat. Recordem-nos què hi havia en aquest fitxer YAML. Així que el batch estava establert a YAML, així que així és com s'hauria d'anomenar, i hauria de tenir un stegosaurus. Així que anem amunt i mirem a results. I tenim COLLECTED-yaml. Així que vegem si tenim un Stegosaurus. Fantàstic, un Stegosaurus portant un barret. Això és el que ens agrada.

Així que això ha funcionat molt bé, i és exactament igual amb el fitxer JSON. Només canviem l'extensió del fitxer aquí i Nextflow sap com llegir-lo.

I en aquest cas, hauríem de tenir un batch anomenat JSON i hauríem de tenir una tortuga. Vegem-ho. Meravellós. Una de les meves eines CLI favorites.

## 2.1. Personalitzar el directori de sortida amb -output-dir

D'acord, així que això ha estat principalment pensant en les entrades al pipeline i canviant paràmetres. Què passa amb les sortides?

Ara, encara que hem estat canviant els subdirectoris utilitzant params, potser heu notat que tots els nostres fitxers encara van a results.

Podem canviar aquest directori base on es publiquen tots els fitxers amb una opció de línia de comandes anomenada -output-dir. Així que si faig nextflow run hello-config, i després faig -output-dir, i l'anomenarem "custom-outdir-cli". No puc escriure. Només perquè recordem d'on venien aquests fitxers.

Aquesta és una opció central de Nextflow i és molt nova. Això només es va afegir recentment, i aquesta és una de les coses que podem fer amb el nou analitzador de llenguatge i tot.

És una mica llarg d'escriure. També podeu anomenar-lo "-o" si voleu. Així que si només torno enrere. Puc escurçar això a "-o", que és una mica més simple.

D'acord. Executem això. No hem canviat res al nostre pipeline o fins i tot a la nostra configuració en aquest punt, i hauria d'esperar desar tots els nostres resultats en un directori de nivell superior diferent. I podeu imaginar que podeu establir això bàsicament a qualsevol camí que vulgueu.

Acaba d'arribar a la part superior. Tenim un custom-outdir-cli, i tots els fitxers estan organitzats allà exactament de la mateixa manera, amb els seus mateixos subdirectoris i noms de fitxer. Així que aquesta és una manera realment fàcil de canviar on el pipeline publica els seus resultats, sense pensar massa en com s'organitzen aquests resultats.

## 2.1.2. Eliminar camins codificats del bloc de sortida

Si miro dins d'aquest directori, podem veure que encara tenim un subdirectori anomenat Hello Config, que sembla una mica redundant ara.

Així que carreguem el nostre script de nou i ara podem eliminar aquest subdirectori del bloc de sortida a la part inferior. Perquè ja no el necessitem realment. Així que podem fer-ho ara, eliminar-lo d'aquí. I després si només és això, podeu eliminar-ho completament o deixar-ho com una cadena buida. Deixaré-ho com una cadena buida per ara, perquè tornarem i posarem algunes coses diferents al seu lloc en el futur. Però si no us importen els subdirectoris, és més net eliminar completament la declaració de camí allà.

D'acord, desem. Provem-ho ràpidament de nou. De fet eliminaré el meu directori "custom-outdir-cli" perquè no ens confonguin els fitxers existents allà. Perquè recordeu, quan publiqueu coses, no elimina els fitxers que hi havia abans. Només afegeix de nous. Executem aquesta comanda de nou, custom-outdir-cli.

I ara si feu "ls custom-outdir-cli", ja no hi ha més directori allà anomenat Hello Config.

## 2.2.1. Establir outputDir al fitxer de configuració

D'acord, l'opció de línia de comandes aquí, "-o" o "-output-dir" està bé. Però què passa amb establir valors per defecte per a això a la configuració? Com ho fem?

Obro el fitxer "nextflow.config", tanco tot el demés i elimino això. Podem afegir una nova opció de configuració aquí, que he copiat del lloc web del material de formació, i s'anomena outputDir.

No està sota cap àmbit. No està sota params ni res. És de nivell superior, i podem establir això a una cadena. Ara una cosa simple a fer és només canviar-ho a qualsevol cosa diferent de results com una cadena codificada. Però perquè això està en un fitxer de configuració de Nextflow, podem ser una mica intel·ligents aquí i també incloure variables.

I podeu veure aquí que hem inclòs una variable params, params.batch, que és part d'aquesta cadena. Això significa que podem reutilitzar variables que vénen d'altres llocs. I en aquest cas, si fem --batch, quan executem el Pipeline de Nextflow, obtindrem un subdirectori al nostre camí personalitzat basat en quin era el nom del batch.

D'acord, així que provem això i fem una ullada ràpida per veure com es veuen els resultats. Així que si faig nextflow run hello-config i --batch my_run. Recordem-nos com es veia la configuració. Així que és custom-outdir-config.

Tree custom-outdir-config. I podeu veure que el batch s'anomenava my_run. I després tenim aquest subdirectori anomenat my_run. Així que aquest camí de fitxer dinàmic va funcionar.

I no només això, ja no va al directori de results per defecte, i no vaig haver d'especificar res a la línia de comandes per canviar el directori base. Així que hem restablert amb èxit el valor per defecte per al outputDir per defecte.

## 2.2.2. Subdirectoris amb noms de batch i procés

D'acord, portem això una mica més lluny. Això és una variable dinàmica dins del fitxer de configuració. Què passa amb el script mateix? Ara, fins ara hem tingut aquests camins aquí i aquests també poden ser dinàmics. Així que en lloc de només codificar alguna cosa, podem posar alguns claudàtors i posar alguna cosa dinàmica.

Així que per exemple, tenim els nostres processos anomenats sayHello. Podríem fer sayHello.name, que és un atribut del procés, que és una mica avorrit perquè és només "sayHello" en aquest cas. Però és variable.

Així que això us dóna una idea. Així que podem posar això aquí i dir convertToUpper.name, collectGreetings.name, collectGreetings.name de nou, i cowpy.

Ara quan executem, el directori base encara serà custom-outdir-config. I serà en un subdirectori anomenat params.batch, però els subdirectoris sota això haurien d'estar organitzats per nom de procés.

Provem això i vegem si funciona. Així que eliminaré el directori anterior perquè no ens confongui, i utilitzaré exactament la mateixa comanda Nextflow Run.

Hauria d'executar-se de la mateixa manera. Podria estar utilitzant dash resume en tots aquests per fer-ho una mica més ràpid i utilitzar els resultats calculats prèviament. Ara, si faig tree custom-outdir-config, podeu veure que no està a results, està al nostre directori base amb el nom del batch. I podeu veure que tots els resultats ara estan organitzats dins de subdirectoris anomenats segons el procés. Així que tenim dos llocs diferents on estem definint camins de sortida dinàmics aquí.

D'acord. Última cosa, afegim de nou aquestes carpetes intermèdies, que teníem abans perquè eren una mica agradables. Intermediates.

I també podem pensar una mica sobre aquest params.batch, potser com a desenvolupador de pipeline realment m'agradava tenir això al subdirectori, però si els usuaris finals del pipeline estan establint "-o" o -output-dir al CLI, està sobreescrivint completament aquesta declaració sencera, i perdem aquest subdirectori.

Així que el que podem fer és podem treure aquest camí dinàmic del outputDir config, que seria esborrat, i posar-lo al camí de sortida, que no és esborrat.

Així que podem fer params.batch barra intermediates barra sayHello.name, i fer tot això en una cadena entre cometes dobles, perquè sigui interpolat per Nextflow.

Ara puc copiar, ups. Copiar aquests cap avall als altres processos. Recordeu posar-los tots entre cometes. I eliminar intermediates d'aquestes sortides particulars.

D'acord? Està semblant lleugerament més complex ara, però podeu veure que realment estem començant a construir una estructura de directori de sortida ben organitzada al nostre codi.

I el que és realment agradable és que aquesta complexitat extra al codi no passa al CLI. Així que podem executar la nostra comanda amb -output-dir i qualsevol variable batch, només pensant en com executar el pipeline i no pensant realment massa en què hi ha al codi. I els nostres fitxers de sortida es construiran realment bé d'una manera molt ben organitzada, que és agradable per a la gent que utilitza el pipeline bàsicament.

Genial. Mentre escric això, m'adono que he comès un error. Vegem si algú m'ha enxampat aquí. Tenim collectGreetings.name, així que alguna cosa ha anat una mica malament. I sí, efectivament, vaig oblidar accidentalment posar aquests entre claudàtors.

Així que recordeu, aneu amb compte quan estigueu escrivint el vostre codi i assegureu-vos que dieu a Nextflow què és una variable i què és només una cadena. Perquè farà exactament el que li digueu que faci. I res més. Com tots els bons ordinadors. D'acord, això hauria d'arreglar-ho.

## 2.3. Establir el mode de publicació a nivell de workflow

Hi ha una part d'aquest script, que encara no m'agrada, que és el fet que estem escrivint mode copy una vegada i una altra, i si hi ha una cosa que no ens agrada, és repetir-nos.

Així que podem netejar això una mica agafant això i movent-ho a la configuració. I de fet, podem establir-ho per a tot el pipeline d'un cop. Així que no hem de dir-ho múltiples vegades.

Anem al nostre fitxer de configuració i tenim un nou àmbit aquí anomenat workflow. I podem fer claudàtors o podem fer notació de punt. No fa cap diferència. El lloc web del material de formació utilitza notació de punt. Puc dir output i podem barrejar i combinar, així que mode equals copy. Genial.

I ara podem tornar aquí i eliminar aquests. Ara podríem deixar-los al seu lloc. La configuració bàsicament està sobreescrivint el que està escrit aquí, però com ho tenim a la configuració de nivell de pipeline, i aquests dos fitxers s'envien junts, no hi ha cap raó per fer-ho realment dues vegades.

D'acord. Només una comprovació de sanitat, perquè aparentment cometem errors. Executem això de nou i només comprovem que estem utilitzant correctament el mode copy per publicar fitxers. Així que executarem el script de nou i aquesta vegada hem posat els resultats en un directori anomenat config-output-mode, vegem com es veuen els fitxers allà.

I després si faig "ls -l" per mirar batch, i podem mirar cowpy, per exemple. I hauríem de veure, sí, que aquest és un fitxer adequat aquí, que no és un enllaç simbòlic, així que aquest atribut de configuració s'ha aplicat correctament.

## 3. Seleccionar una tecnologia d'empaquetament de programari

D'acord. Fins ara ens hem estat centrant en les entrades i les sortides, els fitxers amb els quals s'executa el workflow. Però què passa amb la infraestructura? Vaig dir al principi que Nextflow us permet executar el mateix pipeline en diferents configuracions informàtiques. Així que com es veu això?

Per mostrar això, canviarem d'utilitzar Docker per executar cowpy, i en lloc d'això utilitzarem Conda per fer el mateix.

Puc fer això molt simplement. Si vaig a code, "nextflow.config". Si recordeu a la part superior, vam definir docker.enabled anteriorment, i el capítol anterior perquè poguéssim utilitzar el contenidor amb cowpy dins.

Diré a Nextflow que no utilitzi Docker. Establir això a false. I diré conda enabled equals true. Així que dir a Nextflow, si us plau utilitzeu Conda.

Ara només habilitar Conda no és suficient per si mateix. Igual que vam fer amb Docker, hem de dir a Nextflow on pot obtenir el programari que necessita.

Així que si saltem als mòduls aquí. I obrim l'script de cowpy. Podem veure que tenim una declaració de container a la part superior. I el container és utilitzat per Docker, però també Singularity, Apptainer, i moltes de les altres eines de programari.

Però no pot ser utilitzat per Conda, així que tenim una declaració separada anomenada "conda", i podríem només escriure "cowpy". I això deixarà a la resolució de paquets de conda esbrinar la millor manera de resoldre això, segons el vostre entorn conda local.

O és una bona pràctica fer el que el lloc web del material de formació diu que feu, que és definir un canal conda específic amb la seva notació de doble dos punts, i definitivament definir una versió específica del programari perquè cada persona que executi el pipeline obtingui la mateixa versió.

Noteu que els contenidors són una mica superiors en aquest aspecte, perquè quan instal·leu alguna cosa amb Conda, encara resoldrà totes les dependències per a aquest paquet, i poden canviar amb el temps. S'anomena deriva de dependències.

Així que els contenidors, no obstant això, bloquegen tota la pila de dependències de programari fins al final, així que podeu estar una mica més segurs que A, funcionarà, i B, serà reproduïble.

Així que si podeu utilitzar Docker o Singularity o Apptainer, definitivament ho recomanaria.

Ara el que és agradable d'això és que el fitxer de mòdul, que està escrit pel desenvolupador del pipeline, ara té tant Container com Conda, i així estem dient a la persona que està executant aquest pipeline, no ens importa quina solució d'empaquetament de programari utilitzeu. Funcionarà tant amb Docker com amb Conda, i aquí és on obtenir el programari en ambdós casos.

Podem obrir el terminal i provem això. Així que nextflow run hello-config --batch conda. I la primera vegada que això s'executa amb conda, serà una mica lent quan arribi a aquest procés particular, perquè ha d'executar "conda install".

I està creant un entorn conda especial només per a aquest procés. Així que no està utilitzant el meu entorn conda global, que tinc al meu terminal. Està creant un només per a aquest procés. Això és bo perquè evita coses com xocs de dependències entre diferents processos al vostre workflow. Si els vostres processos tenen eines que necessiten diferents versions de Python o coses així, això està bé perquè estan utilitzant diferents entorns conda.

Nextflow emmagatzema aquests entorns conda localment, podeu veure que us diu on és aquest camí, està al directori de treball aquí. I així la propera vegada que executi aquest script amb Conda, serà molt més ràpid perquè trobarà aquest entorn conda existent i només el reutilitzarà. Però la primera vegada que ho fem, ha d'anar i buscar-lo, resoldre'l, descarregar totes les dependències, i configurar-ho tot.

D'acord, genial, s'ha executat. Podem recordar-nos què està configurat actualment el pipeline per utilitzar. Si mirem al fitxer de configuració, era "custom-outdir-config" ara mateix per a mi. Vegem si vaig fins a aquest directori base. I vaig fer --batch conda. Aquí està el nostre subdirectori conda. Així que va funcionar i aquí està la nostra sortida de cowpy.

Així que va buscar cowpy, el va instal·lar al meu sistema local utilitzant conda, i va executar el procés. I el que és genial és que, com a usuari final, no vaig haver de pensar gens en cap de la gestió de programari allà. Nextflow només ho va resoldre per a mi. Vaig dir, necessito utilitzar conda en aquest sistema. El desenvolupador del pipeline va dir quins paquets necessitava. I Nextflow va fer la resta. Molt potent.

Noteu que podeu utilitzar realment una barreja de diferents tecnologies. Així que puc habilitar Docker per a processos específics, i conda per a altres processos, o dir que alguns processos haurien de només utilitzar qualsevol programari local que tingués instal·lat. Això és força inusual, però és possible, i en alguns casos, per exemple, si esteu utilitzant cert programari que podria ser difícil d'empaquetar a Docker, teniu una escapatòria.

## 4. Seleccionar una plataforma d'execució

Així que això és empaquetament de programari. L'altra part de la portabilitat a altres sistemes és on s'executen realment les tasques. En aquest moment, estic executant bàsicament al meu portàtil o en aquests Codespaces, que és un sol ordinador. No hi ha res sofisticat. Nextflow està sent una mica intel·ligent sobre paral·lelitzar les tasques tan bé com pot, però tot està en un sistema.

Ara, si esteu executant en un HPC, probablement teniu algun tipus de planificador de tasques com SLURM o PBS o alguna cosa, i enviareu tasques a aquest planificador i distribuirà totes les tasques a diferents nodes de càlcul.

Una altra manera d'executar és al núvol. Així que potser esteu utilitzant AWS Batch, o Azure Cloud, o Google. I tots aquests funcionen en un sistema similar on teniu un planificador i envieu tasques i s'envien a diferents llocs per ser calculades.

Ara en el llunyà passat quan vaig començar a fer bioinformàtica, el programari de tothom per executar anàlisis estava molt lligat a la seva infraestructura computacional, que feia gairebé impossible replicar.

Però amb aquesta separació de configuració a Nextflow, i amb la capacitat de Nextflow d'interactuar amb molts backends d'infraestructura de càlcul diferents, és molt simple agafar el nostre pipeline sense modificar el codi del pipeline gens i només canviar això.

## 4.1. Apuntar a un backend diferent

Així que si anem al nostre fitxer "nextflow.config", i ara podem posar alguna configuració de nivell de procés. Així que si poso a la part superior l'àmbit process i puc establir l'executor, i aquí està establert a local, que és el valor per defecte.

Noteu que perquè això és a nivell de procés, podem apuntar coses a diferents processos. I així podeu configurar realment executors per ser específics del procés i tenir una execució híbrida, on algunes tasques podrien executar-se localment, allà on s'està executant la tasca de Nextflow. Algunes s'envien a diferents HPC i algunes podrien enviar-se al núvol. Podeu ser tan intel·ligents com vulgueu.

Ara, és molt difícil demostrar això en un entorn de formació com aquest perquè no tinc un HPC al qual enviar. Però el que puc fer és si escric slurm, podem fer trampes una mica i podeu tenir una sensació d'això.

I això és realment més interessant per a persones que estan acostumades a executar a SLURM i saben com es veuen les capçaleres de SLURM. Però si faig nextflow run, hello-config. Fallarà perquè intentarà enviar tasques a un clúster que no existeix. Així que obtindrem algun tipus d'error sobre sbatch no estar disponible.

Sí, escrit. Aquesta és l'eina. Aquesta és l'eina CLI que utilitzeu per enviar tasques a un clúster slurm. Però el que podem fer és podem anar i mirar al nostre directori de treball aquí fent command clic, obrir aquest directori i mirar el .command.run. I podeu veure a la part superior del fitxer .command.run, tenim les nostres capçaleres sbatch, dient a un clúster SLURM teòric com gestionar aquest enviament de tasca.

I així podeu veure que Nextflow està sent intel·ligent, està fent totes les coses correctes. És només que no teníem un clúster al qual enviar.

## 5. Controlar les assignacions de recursos de càlcul

Què més és diferent entre diferents infraestructures informàtiques? Una altra cosa és quants recursos disponibles teniu, i de fet, en molts entorns de càlcul, és un requisit que hàgiu d'especificar quantes CPUs i quanta memòria necessita una tasca.

De nou, Nextflow abstreu això per a nosaltres, perquè ja no és específic d'un sol tipus d'entorn de càlcul, i podem escriure a l'àmbit de nivell de procés aquí. CPUs equals one, memory equals two gigabytes. El nostre pipeline no és molt exigent, així que això hauria d'estar bé.

Ara, només he endevinat aquests números aquí, però com sabeu quina és una quantitat sensata de recursos a utilitzar? És una feina força difícil anar i excavar a través de tots aquests diferents processos d'un gran pipeline de moltes mostres i entendre quina va ser la utilització de recursos.

Així que un bon enfocament per a això és establir aquests valors a números alts per començar, només perquè el vostre pipeline s'executi sense cap error, i després demanar a Nextflow que generi un informe d'ús per a vosaltres.

Això és súper fàcil de fer, així que tornaré a un terminal. Oh, necessito recordar establir això de nou a local perquè el meu pipeline realment s'executi. I diré nextflow run, i utilitzaré una opció de línia de comandes -with-report.

I puc deixar això en blanc i donarà un nom de fitxer per defecte, però li donaré un nom de fitxer específic perquè això es desi a un lloc específic.

Premo Enter, i el pipeline s'executa exactament com de normal, però quan acabi, generarà un bon informe HTML per a mi.

Així que a la barra lateral aquí, tinc aquest fitxer HTML. Si estigués executant això localment, només l'obriria. Estic, perquè estic a Codespaces, faré clic dret sobre això i faré clic a descarregar, que el descarregarà al meu ordinador local. I puc obrir-lo fàcilment al navegador web.

Nextflow pot generar un informe com aquest per a qualsevol pipeline i té informació realment agradable. Així que és una bona pràctica desar sempre aquestes coses. Ens diu quan vam executar, on vam executar, si va tenir èxit o no, quins paràmetres es van utilitzar, quina va ser la comanda CLI, coses així.

I també hi ha aquests gràfics sobre l'ús de recursos. Així que ens diu quin percentatge de crides de CPU es van utilitzar per a cada procés com un gràfic de caixa aquí, perquè hi ha moltes tasques per a cada procés, així que podem veure la distribució.

Podeu veure els nostres processos aquí, cowpy i collectGreetings només tenien una sola tasca, així que és només una sola línia. I tenim tant CPU com memòria i durada de la tasca, i van ser molt ràpids.

Si esteu utilitzant Seqera Platform, per cert, obteniu els mateixos gràfics integrats a la interfície de Platform sense haver de fer res. Així que sempre teniu aquesta informació a l'abast.

D'acord, així que podem utilitzar aquest informe i en una execució real, i tenir una sensació de quantes CPUs i quanta memòria està sent utilitzada pel nostre pipeline i tornar i posar aquests valors de nou al nostre fitxer de configuració, perquè la propera vegada potser no demanem tant. I podem ser una mica més eficients.

Ara podeu ser realment intel·ligents sobre configurar fitxers de configuració de pipeline. I de nou, si esteu utilitzant Seqera Platform, busqueu un petit botó que sembla una bombeta. Perquè si cliqueu això, generarà un fitxer de configuració altament optimitzat, que està adaptat específicament per a les vostres dades, la vostra execució i el vostre pipeline. Per executar-lo de la manera més eficient possible.

Però per ara, diré que en realitat el nombre per defecte de CPUs que Nextflow estava donant estava bé i però només necessitem un gigabyte de memòria.

## 5.3. Establir assignacions de recursos per a un procés específic

Ara, a la vida real, és força inusual que tots els processos del vostre pipeline necessitin els mateixos requisits. Podríeu tenir alguna cosa com MultiQC com a eina d'informes, que necessita molt poc en termes de recursos i s'executa força ràpidament.

I després potser teniu alguna cosa que està indexant un genoma de referència o fent algun alineament o fent alguna altra tasca. No importa què sigui, que pren molts recursos. I així per a aquests diferents enviaments de tasques a un planificador, voleu donar diferents quantitats de recursos.

Sota aquest àmbit de procés, podem definir una configuració, que apunta a processos específics de diferents maneres.

Aquí estem utilitzant withName, també podem utilitzar etiquetes, i aquestes poden utilitzar un patró per apuntar a un o múltiples processos. Aquí només estem dient qualsevol procés que tingui un nom cowpy establert a dos gigabytes de memòria i dues CPUs, i perquè aquest és un selector més específic que el procés de nivell superior, això és sobreescrit en aquests casos, així que podeu construir un bon fitxer de configuració aquí, que realment adapta tots els vostres diferents processos al vostre pipeline per fer-los realment eficients.

## 5.5. Afegir límits de recursos

Ara com a desenvolupador de pipeline, probablement conec les eines força bé, i vull que tot s'executi tan ràpid i tan bé com sigui possible. Així que podria ser que posi números força alts per a alguns d'aquests perquè sé que s'executarà molt més ràpid si dono a cowpy 20 CPUs en lloc de dues.

Això està bé fins que aneu a executar al vostre portàtil o a GitHub Actions Continuous Integration test, o algun altre sistema, que potser no té 20 CPUs disponibles.

Ara quan intenteu executar el pipeline, fallarà perquè Nextflow dirà, no puc enviar aquesta tasca a cap lloc. No tinc els recursos disponibles.

Ara per evitar aquest error dur, podem afegir una mica més de configuració, que és específica del nostre sistema ara, anomenada límits de recursos. I això es veu així. Està sota l'àmbit de procés de nou.

I límits de recursos, podeu especificar bàsicament el sostre del que teniu disponible. És un mapa aquí, i podeu, dins d'aquest mapa, podeu establir la memòria, les CPUs, i el temps.

Ara el que passa és quan Nextflow envia una tasca d'un procés, mira el que es demana i bàsicament només fa un mínim entre això i això. Així que si vam demanar 20 CPUs, però només quatre estan disponibles, demanarà quatre. El pipeline no falla i utilitza tan a prop del que va ser dissenyat pel desenvolupador del pipeline com sigui possible.

## 6. Utilitzar perfils per canviar entre configuracions predefinides

D'acord. Vaig dir que els límits de recursos aquí podrien ser específics del sistema, i potser tinc un fitxer de configuració de Nextflow al meu pipeline, i sé que la gent utilitzarà això en una gamma de llocs diferents. Ara, en lloc de forçar tothom a crear el seu propi fitxer de configuració de Nextflow cada vegada, el que puc fer és que puc agrupar diferents preestabliments de configuració junts en perfils de configuració.

Baixaré una mica aquí i després just després de params, perquè l'ordre del fitxer de configuració aquí és important, el fitxer de configuració es carrega seqüencialment, així que posaré aquests perfils després de tot el demés perquè sobreescrigui els params definits prèviament. I enganxaré aquests perfils del material de formació.

Així que hi ha un nou àmbit de nivell superior anomenat profiles. Podem tenir noms arbitraris aquí. Així que tenim my_laptop i univ_hpc. I aquí podem veure que estem establint els mateixos paràmetres de configuració que estàvem abans. Ara només dins d'un perfil. Així que tenim un executor local per executar al meu portàtil i estic enviant a un clúster SLURM a l'HPC.

Estic utilitzant Docker localment, conda a l'HPC, i el sistema HPC té límits de recursos molt més alts.

Ara puc executar el pipeline amb l'opció CLI -profile, dir quin perfil vull utilitzar. Així que utilitzaré my_laptop, i Nextflow aplicarà tota la configuració dins d'aquest àmbit de perfil. Així que puc provar això ara. És la mateixa comanda que abans. Nextflow run hello-config, i faig dash profile, guió únic perquè és l'opció central de Nextflow, dash profile my_laptop.

Ara aplicarà per lots tota aquesta opció de configuració. Oh, i podeu veure, vaig dir abans que això podria passar que el requisit del procés, va demanar quatre CPUs i només en tinc dues en aquesta instància de Codespaces.

Així que aquesta és una bona oportunitat només per provar els límits de recursos del procés, i dir que només tinc dues CPUs al meu portàtil, o en aquests Codespaces. Ara si ho executem de nou, hauria de limitar aquest requisit a dues i esperem que el pipeline s'executi. Genial.

## 6.2. Crear un perfil de paràmetres de prova

Noteu que aquests perfils no han de tenir només configuració sobre la seva infraestructura. Podeu tenir agrupacions de qualsevol configuració aquí, incloent paràmetres.

Així que una altra cosa que veureu molt sovint als pipelines de la gent és un perfil de prova, que inclou paràmetres, que normalment enviaríeu per usuari. Però aquí tenim, bàsicament diferents valors per defecte sensats per quan vull executar casos de prova.

I això és genial perquè no he d'anar necessàriament i especificar totes aquestes coses, que podrien ser paràmetres requerits. Altrament puc només dir dash profile test i s'executarà directament.

Ara una cosa a notar és que els perfils també es poden combinar més d'un. Així que puc fer profile my_laptop aquí, i després també afegir test. No faig profile dues vegades. Només faig una llista separada per comes aquí sense espais. I aplicarà aquests perfils en ordre. Així que prendrà la configuració del perfil my_laptop, i després aplicarà la configuració de test a sobre.

Realment convenient i podeu veure com podeu configurar molts grups per defecte sensats aquí per facilitar l'execució del vostre pipeline.

## 6.3. Utilitzar nextflow config per veure la configuració resolta

Esperem que us hagi convençut que la resolució de configuració de Nextflow és potent, però no us culparia si esteu anant una mica amb els ulls creuats en aquest punt després que hagi dit unes 20 maneres diferents de proporcionar configuració i donar totes aquestes diferents capes com una pell de cebolla.

Així que si mai us sentiu insegurs sobre quina és la configuració resolta final per a Nextflow, sabeu que hi ha una comanda anomenada "nextflow config", i podem executar això i ens dirà quina és la configuració resolta a la nostra ubicació actual.

Així que quan l'executo aquí, troba el fitxer "nextflow.config" al directori de treball actual, i processa tota la configuració diferent, i em dóna la sortida resolta.

Noteu que el fitxer de configuració de Nextflow també pot prendre l'opció CLI de perfil. Així que si li dic que resolgui en els perfils my_laptop i test, i podeu veure que també ha aplicat els límits de recursos aquí de l'opció de configuració my_laptop i també ha establert els params, que estaven a la prova.

Així que aquesta és una manera agradable només d'explorar com està funcionant la resolució de configuració, si teniu algun dubte.

## Conclusió

D'acord, això és tot. Això és la configuració de Nextflow en poques paraules. Podeu fer moltes coses amb la configuració. És realment potent. Però aquests són la majoria dels casos d'ús comuns que us trobareu fent, i aquests conceptes s'apliquen a totes les diferents opcions.

Doneu-vos una palmadeta a l'esquena perquè aquest és el final del curs de formació Hello Nextflow. Esperem que ara tingueu confiança tant en escriure el vostre propi pipeline de Nextflow des de zero, configurar-lo i executar-lo, i coneixeu tots els detalls i les coses a tenir en compte.

Hi ha un qüestionari més que podeu provar a la pàgina de formació de configuració. Així que baixeu i proveu-ho i assegureu-vos que heu entès totes aquestes parts sobre la configuració.

I uniu-vos a nosaltres a l'últim vídeo només per a una ràpida conclusió sobre alguns dels propers passos que podrien ser bons de fer després d'aquest curs de formació.

Gràcies per quedar-vos amb nosaltres. Ben fet i us veuré al proper vídeo.
