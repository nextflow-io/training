# Başlarken

## Eğitim ortamını başlatın

GitHub Codespaces üzerinde sağladığımız önceden hazırlanmış ortamı kullanmak için aşağıdaki "GitHub Codespaces'te Aç" düğmesine tıklayın. Diğer seçenekler için [Ortam seçenekleri](../envsetup/index.md) sayfasına bakın.

Eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (ekipmanınıza bağlı olarak sağ tıklama, ctrl-tıklama veya cmd-tıklama kullanın), böylece ortam yüklenirken okumaya devam edebilirsiniz.
Kursu tamamlamak için bu talimatları paralel olarak açık tutmanız gerekecektir.

[![GitHub Codespaces'te Aç](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam temelleri

Bu eğitim ortamı, eğitim kursunda çalışmak için gerekli tüm yazılımları, kodu ve veriyi içerir, bu nedenle kendiniz herhangi bir şey yüklemenize gerek yoktur.

Codespace, bir dosya sistemi gezgini, bir kod düzenleyici ve bir terminal kabuğu içeren bir VSCode arayüzü ile kurulmuştur.
Kurs boyunca verilen tüm talimatlar (örneğin 'dosyayı açın', 'kodu düzenleyin' veya 'bu komutu çalıştırın'), aksi belirtilmedikçe VSCode arayüzünün bu üç bölümüne atıfta bulunur.

Bu kursu kendi başınıza çalışıyorsanız, lütfen daha fazla ayrıntı için [ortam temelleri](../envsetup/01_setup.md) ile kendinizi tanıştırın.

### Sürüm gereksinimleri

Bu eğitim, **Nextflow 25.10.2** veya sonrası için **v2 sözdizimi ayrıştırıcısı DEVRE DIŞI** olacak şekilde tasarlanmıştır.

#### Eğitim ortamımızı kullanıyorsanız:

Daha ileri gitmeden önce aşağıdaki komutu çalıştırmanız GEREKİR:

```bash
export NXF_SYNTAX_PARSER=v1
```

#### Yerel veya özel bir ortam kullanıyorsanız:

Lütfen [burada](../info/nxf_versions.md) belgelenen doğru ayarları kullandığınızdan emin olun.

Eğitim ayrıca **nf-core tools 3.4.1** gerektirir.
Farklı bir nf-core araç sürümü kullanırsanız, takip etmekte zorluk yaşayabilirsiniz.

Ortamınızda hangi sürümün yüklü olduğunu `nf-core --version` komutunu kullanarak kontrol edebilirsiniz.

## Çalışmaya hazırlanın

Codespace'iniz çalışmaya başladığında, eğitime dalmadan önce yapmanız gereken iki şey var: bu özel kurs için çalışma dizininizi ayarlayın ve sağlanan materyallere bir göz atın.

### Çalışma dizinini ayarlayın

Varsayılan olarak, codespace tüm eğitim kurslarının kök dizininde ayarlanmış çalışma dizini ile açılır, ancak bu kurs için `hello-nf-core/` dizininde çalışacağız.

Terminalde şu komutu çalıştırarak şimdi dizini değiştirin:

```bash
cd hello-nf-core/
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örneğin codespace'iniz uyku moduna geçerse), Github Codespaces eğitim ortamında çalıştığınızı varsayarak, tam yolu kullanarak her zaman geri dönebilirsiniz:

    ```bash
    cd /workspaces/training/hello-nf-core
    ```

Şimdi bu dizinin içeriğine bir göz atalım.

### Sağlanan materyalleri keşfedin

Bu dizinin içeriğini, eğitim çalışma alanının sol tarafındaki dosya gezginini kullanarak keşfedebilirsiniz.
Alternatif olarak, `tree` komutunu kullanabilirsiniz.

Kurs boyunca, dizin yapısını ve içeriğini okunabilir bir biçimde temsil etmek için `tree` çıktısını kullanırız, bazen netlik için küçük değişikliklerle.

Burada ikinci seviyeye kadar bir içindekiler tablosu oluşturuyoruz:

```bash
tree . -L 2
```

??? abstract "Dizin içeriği"

    ```console
    .
    ├── greetings.csv
    ├── original-hello
    │   ├── hello.nf
    │   ├── modules
    │   └── nextflow.config
    └── solutions
        ├── composable-hello
        ├── core-hello-part2
        ├── core-hello-part3
        ├── core-hello-part4
        ├── core-hello-part5
        └── core-hello-start
    ```

Bölümü genişletmek ve içeriğini görüntülemek için renkli kutuya tıklayın.
Beklenen komut çıktısını özlü bir şekilde dahil etmek için bunun gibi daraltılabilir bölümler kullanırız.

- **`greetings.csv` dosyası**, test amaçları için kullandığımız bazı minimal sütunlu veriler içeren bir CSV'dir.

- **`original-hello` dizini**, Hello Nextflow eğitim serisinin tamamında çalışılarak üretilen kaynak kodun bir kopyasını içerir (Docker etkin).

- **`solutions` dizini**, kursun her adımından kaynaklanan tamamlanmış iş akışı betiklerini içerir.
  Çalışmanızı kontrol etmek ve herhangi bir sorunu gidermek için referans olarak kullanılmaları amaçlanmıştır.

## Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu mu düşünüyorsunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Ortamım çalışır durumda
- [ ] Sözdizimi ayrıştırıcısının **v1** olarak ayarlandığından emin oldum
- [ ] Çalışma dizinini uygun şekilde ayarladım

Tüm kutuları işaretleyebiliyorsanız, hazırsınız demektir.

**Bölüm 1'e devam etmek için bu sayfanın sağ alt köşesindeki oka tıklayın.**
