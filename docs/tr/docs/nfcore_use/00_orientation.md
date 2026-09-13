---

# Başlarken

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Eğitim ortamını başlatın

GitHub Codespaces üzerinde sunduğumuz hazır ortamı kullanmak için aşağıdaki "Open in GitHub Codespaces" düğmesine tıklayın. Diğer seçenekler için [Ortam seçenekleri](../envsetup/index.md) sayfasına bakın.

Ortam yüklenirken okumaya devam edebilmek için eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (kullandığınız cihaza göre sağ tıklama, ctrl+tıklama veya cmd+tıklama kullanabilirsiniz). Kursu tamamlamak için bu talimatları paralel olarak açık tutmanız gerekecektir.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam temelleri

Bu eğitim ortamı, kursu tamamlamak için gereken tüm yazılım, kod ve verileri içermektedir; bu nedenle herhangi bir şey kurmanıza gerek yoktur.

Codespace, VSCode arayüzüyle yapılandırılmıştır; bu arayüz bir dosya sistemi gezgini, bir kod düzenleyici ve bir terminal kabuğu içerir. Kurs boyunca verilen tüm talimatlar (örneğin "dosyayı açın", "kodu düzenleyin" veya "bu komutu çalıştırın"), aksi belirtilmedikçe VSCode arayüzünün bu üç bölümüne atıfta bulunur.

Bu kursu kendi başınıza tamamlıyorsanız, daha fazla ayrıntı için lütfen [ortam temelleri](../envsetup/01_setup.md) sayfasını inceleyin.

### Sürüm gereksinimleri

Bu eğitim, Nextflow 25.10.2 veya sonraki sürümlerle **v2 sözdizimi ayrıştırıcısıyla** birlikte çalışmaktadır; v2 ayrıştırıcısı, Nextflow 26.04'ten itibaren varsayılan olarak kullanılmaktadır. Eğitim ortamımızda herhangi bir işlem yapmanıza gerek yoktur: ortam, v2 ayrıştırıcısıyla birlikte Nextflow 26.04.4 çalıştırmaktadır. Yerel veya özel bir ortam kullanıyorsanız [sürüm notlarına](../info/nxf_versions.md) bakın.

!!! warning "nf-core/demo, Nextflow 25.10.4 veya sonrasını gerektirir"

    Part 1'de kullanılan `nf-core/demo` pipeline'ı, genel eğitim tabanı olan 25.10.2'den daha katı olan kendi minimum Nextflow sürümünü (`>=25.10.4`) zorunlu kılmaktadır.
    Eğitim ortamımız bu gereksinimi zaten karşılamaktadır; yerel veya özel bir ortam kullanıyorsanız Nextflow 25.10.4 veya sonraki bir sürümde olduğunuzdan emin olun.

Bu eğitim ayrıca **nf-core tools 4.0.2** gerektirmektedir. Farklı bir nf-core araç sürümü kullanıyorsanız kursu takip etmekte güçlük çekebilirsiniz.

Ortamınızda hangi sürümün yüklü olduğunu `nf-core --version` komutunu kullanarak kontrol edebilirsiniz.

!!! warning "v2 ayrıştırıcı uyumluluğu"

    Birçok nf-core pipeline'ı henüz v2 sözdizimi ayrıştırıcısını desteklememektedir.
    Bu kursta kullanılanlar dışında bir nf-core pipeline'ı çalıştırır ve hatalarla karşılaşırsanız, `export NXF_SYNTAX_PARSER=v1` ayarını yaparak v1 ayrıştırıcısına geçmeniz gerekebilir.
    Ayrıntılar için [sürüm notlarına](../info/nxf_versions.md) bakın.

## Çalışmaya hazırlanın

Codespace'iniz çalışmaya başladıktan sonra, eğitime dalmadan önce yapmanız gereken iki şey vardır: bu kursa özgü çalışma dizinini ayarlamak ve sağlanan materyallere göz atmak.

### Çalışma dizinini ayarlayın

Codespace varsayılan olarak tüm eğitim kurslarının kök dizininde açılır; ancak bu kurs için `nfcore-use/` dizininde çalışacağız.

Terminalde şu komutu çalıştırarak dizini şimdi değiştirin:

```bash
cd nfcore-use/
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örneğin codespace'iniz uyku moduna geçerse), Github Codespaces eğitim ortamında çalıştığınızı varsayarak tam yolu kullanarak her zaman geri dönebilirsiniz:

    ```bash
    cd /workspaces/training/nfcore-use
    ```

Ardından bu dizinin içeriğini inceleyin.

### Sağlanan materyalleri keşfedin

Bu dizinin içeriğini, eğitim çalışma alanının sol tarafındaki dosya gezginini kullanarak inceleyebilirsiniz. Alternatif olarak `tree` komutunu kullanabilirsiniz.

```bash
tree .
```

??? abstract "Dizin içeriği"

    ```console
    .
    ├── custom.config
    ├── laptop.config
    ├── malformed_samplesheet.csv
    └── my_params.yml
    ```

- **`laptop.config` dosyası**, üretim ölçeğinde bir pipeline'ı yerel olarak çalıştırırken kaynak kullanımını sınırlamak amacıyla 4. bölümde kullanacağımız bir yapılandırma dosyasıdır. O zamana kadar görmezden gelebilirsiniz.
- **`my_params.yml`, `malformed_samplesheet.csv` ve `custom.config` dosyaları**, Part 2'de parametrelerin bir dosyadan ayarlanmasını, girdi doğrulamasını ve süreç düzeyinde yapılandırma geçersiz kılmalarını göstermek için kullanılmaktadır. Bunları da o zamana kadar görmezden gelebilirsiniz.

## Hazırlık kontrol listesi

Başlamaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Ortamım çalışıyor
- [ ] nf-core tools 4.0.2 kullanıyorum (`nf-core --version` ile kontrol edin)
- [ ] Çalışma dizinimi uygun şekilde ayarladım

Tüm kutuları işaretleyebildiyseniz, başlamaya hazırsınız.

**Part 1'e devam etmek için bu sayfanın sağ alt köşesindeki oka tıklayın.**
