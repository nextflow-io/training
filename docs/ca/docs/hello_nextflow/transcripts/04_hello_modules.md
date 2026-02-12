# Part 4: Hello Modules - Transcripció del vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importants"

    Aquesta pàgina mostra només la transcripció. Per a instruccions completes pas a pas, torneu al [material del curs](../04_hello_modules.md).

    Els números de secció mostrats a la transcripció es proporcionen només amb finalitats indicatives i poden no incloure tots els números de secció dels materials.

## Benvinguda

Hola, i benvinguts de nou a la part quatre de Hello Nextflow. Aquesta secció tracta sobre mòduls, i és una secció força curta del curs. En realitat no escriurem gaire codi, es tracta més de com organitzem el codi al nostre pipeline.

Fins ara, hem estat posant tot en un únic fitxer, cosa que està bé, i de fet així és com solíem construir pipelines de Nextflow en els vells temps.

Però a mesura que el pipeline creix, l'script es fa més i més llarg i cada cop més difícil de navegar, de mantenir, i també significa que no podem compartir realment cap del codi.

Els mòduls de Nextflow ens permeten extreure processos de l'script principal i després importar-los. Això significa que el codi és més fàcil de navegar i també significa que podem compartir aquest codi de mòdul entre diferents pipelines.

Aquest petit diagrama a la pàgina principal de la documentació mostra el concepte de manera clara. En lloc d'un script enorme, inclourem aquests fitxers de mòdul separats, de diferents scripts de mòdul, i tot s'incorporarà al workflow, però encara s'executarà exactament de la mateixa manera.

Així que anem a GitHub Codespaces i fem una ullada. Com abans, he netejat una mica el meu espai de treball aquí. He eliminat els directoris antics de Nextflow i el directori work i així successivament. Però no importa si encara teniu aquests fitxers.

Començaré a treballar al fitxer hello modules, que bàsicament és on el vam deixar al final del capítol anterior. Tenim els nostres tres processos aquí. Tenim un parell de params, el bloc workflow, on estem executant aquests tres processos i unint-los amb canals. Després publiquem els canals de sortida i tenim el bloc output que diu com publicar aquests fitxers.

## 1. Crear un directori per emmagatzemar mòduls

Ara, com he dit, en realitat no escriurem ni editarem gaire codi. Simplement mourem el codi que ja tenim. Els fitxers de mòdul de Nextflow normalment tenen un únic procés, i per convenció normalment els mantenim en un directori anomenat modules. Però podeu anomenar-lo com vulgueu. Però mantindré un directori modules al meu repositori aquí, i després crearé un fitxer per a cada procés. Així que diré new file, sayHello.nf.

## 2. Crear un mòdul per a sayHello()

Ara agafaré el meu procés i simplement seleccionaré aquest codi, el tallaré del fitxer principal hello modules i l'enganxaré aquí.

Òbviament això no fa res per si mateix. El nostre script principal encara necessita aquest procés, així que hem de tornar-lo a incorporar d'alguna manera. I ho fem amb la declaració include.

Així que escric include i unes claus, i després prenc el nom del procés. I dic from, i després dono una ruta de fitxer relativa. Així que diu, comença amb ./ perquè és relatiu des d'on es desa aquest script. Així que és modules sayHello.nf.

Noteu que l'extensió de VS Code és força útil aquí. Ens diu si pot trobar aquest fitxer i si pot trobar un procés, que estic anomenant. Si deliberadament poso un error tipogràfic aquí, em dóna un error immediatament i em dirà que no pot trobar aquest procés que estic intentant importar. Així que estigueu atents a qualsevol error que trobeu.

I això és realment tot. Encara tenim el nostre procés aquí. No calen canvis aquí a baix. El procés té el mateix nom i s'executa exactament de la mateixa manera. Només que el codi real del procés ara està en un fitxer separat.

Podem executar el workflow de Nextflow de nou, funcionarà exactament de la mateixa manera. I això és bàsicament la resta d'aquest capítol del curs, només moure aquests tres processos als seus propis fitxers.

Així que fem-ho ara. Crearé ràpidament un nou fitxer de mòdul per al segon procés: convertToUpper.nf. Tallaré aquest codi, l'enganxaré aquí. I després inclouré aquest. Genial.

I després crearé un nou fitxer per a collectGreetings.nf. Tallaré això.

Molt de tallar, tallar i copiar i enganxar.

I ara el nostre script de workflow principal de sobte sembla molt, molt més curt, molt més accessible i molt més fàcil de llegir.

I podeu veure com el projecte ara comença a construir-se amb els nostres diferents fitxers. Podem aprofundir en el detall als llocs que volem. Navegar pel nostre camí per trobar passos específics al pipeline molt més fàcilment, i obtenir una visió general del que fa el pipeline ràpidament.

## Navegar pels mòduls amb VS Code

Ara, per descomptat, l'inconvenient de fer això és que si teniu un pipeline gran, tindreu molts fitxers de mòdul i podrien estar organitzats en múltiples subdirectoris o tot tipus de coses. Ara, de nou, un petit consell aquí. L'extensió de VS Code és força bona navegant per la vostra base de codi i també informant-vos sobre el codi allà.

Podeu veure que VS Code entén què és aquest procés i em dóna una petita visió general quan passo el cursor per sobre, així que puc veure sense haver d'anar a buscar el codi font, quines són les entrades i les sortides, que normalment és el més important quan l'estic utilitzant en un workflow.

I també si mantinc premut command, estic en un Mac, i faig clic al nom del procés, obre el fitxer directament immediatament. L'incorpora. Així que puc saltar directament allà sense ni tan sols pensar en quines són les rutes de fitxer reals. I això funciona a qualsevol lloc, també ho puc fer on s'estan cridant els processos. Així que això fa que sigui molt ràpid.

## 4.4. Executar el workflow

D'acord, comprovem simplement que el pipeline encara s'executa com esperem. Així que obro el terminal. Fem "nextflow run hello modules", i vegem si s'executa sense cap problema.

Esperem que tot el sentit d'això sigui que el pipeline bàsicament no ha canviat, així que no hauríeu de veure realment cap canvi respecte a quan el vam executar abans. La sortida aquí sembla exactament la mateixa, i podeu veure el nostre directori results amb tots els mateixos fitxers, així que això és genial. Cap canvi és bo.

## Una nota sobre nf-core/modules

Just abans d'acabar, vull tocar breument el poder de la col·laboració quan es tracta de mòduls. Aquests fitxers estan al meu repositori, així que no és obvi immediatament com podríem col·laborar-hi. I hi ha moltes maneres diferents de fer-ho, però probablement l'exemple més gran i més conegut d'això és nf-core.

Si vaig al lloc web d'nf-core, vaig a resources, i modules. Podeu veure que nf-core té una enorme biblioteca de mòduls, gairebé just per sota de 1700 mòduls quan ho veig. I així puc escriure el nom de qualsevol de les meves eines favorites, anar i trobar si algú altre ja ha escrit un mòdul per a això, i veure aquest procés de mòdul pre-escrit aquí amb totes les entrades, les sortides, els contenidors de programari, tota aquesta informació, i podeu veure al costat aquí, quants pipelines diferents d'nf-core estan tots utilitzant aquest únic procés compartit.

Aquest és un exemple una mica extrem, però podeu veure que això és realment reutilitzar aquest codi. I si faig clic a la font de GitHub per a això, és exactament el mateix que estem fent. És només un gran procés en un fitxer.

Ara al costat d'nf-core, fem alguns trucs per poder compartir aquests fitxers i portar-los a diferents repositoris. I si voleu saber més sobre això, aneu i consulteu el curs que tenim sobre l'ús i la construcció amb nf-core específicament. Però volia només donar-vos una idea de fins a quin punt pot ser poderós aquest concepte de reutilització de codi.

## Conclusió

Bé, això és tot per als mòduls. Us vaig dir que era una secció curta del curs. Consulteu el qüestionari, assegureu-vos que ho enteneu i assegureu-vos que tot encara funciona correctament. I us veuré de nou al proper vídeo, que tracta tot sobre contenidors de programari. Moltes gràcies.
