# Başlangıç

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Eğitim ortamını başlatın

GitHub Codespaces'te sağladığımız önceden oluşturulmuş ortamı kullanmak için aşağıdaki "Open in GitHub Codespaces" butonuna tıklayın. Diğer seçenekler için [Ortam seçenekleri](../envsetup/index.md) bölümüne bakın.

Ortam yüklenirken okumaya devam edebilmeniz için eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (cihazınıza bağlı olarak sağ tıklama, ctrl-tıklama veya cmd-tıklama kullanın).
Kurs boyunca çalışmak için bu talimatları paralel olarak açık tutmanız gerekecektir.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam temelleri

Bu eğitim ortamı, eğitim kursunda çalışmak için gerekli tüm yazılım, kod ve verileri içerir; böylece kendiniz hiçbir şey yüklemeniz gerekmez.

Codespace, dosya sistemi gezgini, kod editörü ve terminal shell içeren bir VSCode arayüzüyle kurulmuştur.
Kurs sırasında verilen tüm talimatlar (örn. 'dosyayı açın', 'kodu düzenleyin' veya 'bu komutu çalıştırın') aksi belirtilmedikçe VSCode arayüzünün bu üç bölümüne atıfta bulunur.

Bu kursu kendi başınıza çalışıyorsanız, daha fazla ayrıntı için [ortam temelleri](../envsetup/01_setup.md) ile tanışın.

### Versiyon gereksinimleri

Bu kurs, v2 sözdizimi ayrıştırıcısı etkin (25.10+ sürümünde varsayılan) olarak Nextflow 25.10.2 veya üstünü gerektirir.
Yerel veya özel bir ortam kullanıyorsanız, [burada](../info/nxf_versions.md) belgelendiği gibi doğru ayarları kullandığınızdan emin olun.

## Çalışmaya hazırlanın

Codespace'iniz çalışmaya başladığında, eğitime dalmadan önce yapmanız gereken iki şey var: çalışma dizininizi ayarlayın ve sağlanan materyallere göz atın.

### Çalışma dizinini ayarlayın

Varsayılan olarak, codespace tüm eğitim kurslarının kökünde açılır.
Bu kurs için `nextflow-run/` dizinine geçin:

```bash
cd nextflow-run/
```

Ardından VSCode'u bu dizine odaklanacak şekilde ayarlayın; böylece dosya gezgini kenar çubuğunda yalnızca ilgili dosyalar görünür:

```bash
code .
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örn. codespace'iniz uyursa), Github Codespaces eğitim ortamında çalıştığınızı varsayarak, ona geri dönmek için her zaman tam yolu kullanabilirsiniz:

    ```bash
    cd /workspaces/training/nextflow-run
    ```

### Sağlanan materyalleri keşfedin

Kurs materyallerini soldaki dosya gezginini kullanarak veya `tree` komutuyla keşfedebilirsiniz.
Tam yapıyı görmek için terminalden şunu çalıştırın:

```bash
tree . -L 2
```

??? abstract "Dizin içeriği"

    ```console
    .
    ├── 1-hello.nf
    ├── 2-inputs.nf
    ├── data
    │   ├── greetings-extended.csv
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

**`.nf` dosyaları**, artan karmaşıklıkta iş akışı betikleridir ve kurs boyunca bu sırayla kullanılır.

**`data/`** dizini, 2. bölümden itibaren kullanacağımız CSV girdi dosyalarını içerir.

**`modules/`** dizini, `main.nf` tarafından kullanılan süreç tanımlarını içerir.

**`nextflow.config`** dosyası, minimum ortam özelliklerini ayarlayan bir yapılandırma dosyasıdır. Şimdilik görmezden gelebilirsiniz; 4. bölümde üzerinden geçeceğiz.

## Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu mu düşünüyorsunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Ortamım çalışıyor
- [ ] Çalışma dizinimi uygun şekilde ayarladım

Tüm kutuları işaretleyebiliyorsanız, hazırsınız.

**[Bölüm 1: Nextflow'u Çalıştırma](./01_run_nextflow.md) bölümüne devam etmek için bu sayfanın sağ alt köşesindeki oka tıklayın.**
