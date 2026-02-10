---
title: Hello pipeline'ı
description: Hello pipeline'ının ne yaptığı ve nasıl yapılandırıldığının özeti.
hide:
  - toc
  - footer
---

# Hello pipeline'ı

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Eğitim kurslarımızın çoğu, Nextflow kavramlarını ve mekanizmalarını göstermek için basit, alan bağımsız bir pipeline kullanır.
Hello Nextflow kursu, bu pipeline'ı her tasarım ve uygulama kararını açıklayan adım adım bir şekilde nasıl geliştireceğinizi gösterir.
Diğer eğitimler bu pipeline'ı veya parçalarını başlangıç noktası olarak kullanır.

Bu sayfa, Hello Nextflow kursunun tamamlanmasındaki pipeline'ın durumunu özetlemektedir.

### Özet açıklama

Hello workflow'u, selamlamalar içeren bir CSV dosyası alır, bunları ayrı dosyalara yazar, her birini büyük harfe dönüştürür, hepsini tekrar bir araya toplar ve selamlamaları söyleyen eğlenceli bir karakterin ASCII resmini içeren tek bir metin dosyası çıktılar.

### Workflow adımları (process'ler)

Dört adım, ayrı modül dosyalarında saklanan Nextflow process'leri (`sayHello`, `convertToUpper`, `collectGreetings` ve `cowpy`) olarak uygulanmıştır.

1. **`sayHello`:** Her selamlamayı kendi çıktı dosyasına yazar (örn. "Hello-output.txt")
2. **`convertToUpper`:** Her selamlamayı büyük harfe dönüştürür (örn. "HELLO")
3. **`collectGreetings`:** Tüm büyük harfli selamlamaları tek bir toplu dosyada toplar
4. **`cowpy`:** `cowpy` aracını kullanarak ASCII sanatı oluşturur

### Diyagram

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### Sonuçlar

Sonuçlar `results/` adlı bir dizine yayınlanır ve pipeline'ın son çıktısı (varsayılan parametrelerle çalıştırıldığında), büyük harfli selamlamaları söyleyen bir hindinin ASCII sanatını içeren düz metin dosyasıdır.

```txt title="results/cowpy-COLLECTED-test-batch-output.txt"
  _________
/ BONJOUR \
| HELLO   |
\ HOLà    /
---------
  \                                  ,+*^^*+___+++_
  \                           ,*^^^^              )
    \                       _+*                     ^**+_
    \                    +^       _ _++*+_+++_,         )
              _+^^*+_    (     ,+*^ ^          \+_        )
            {       )  (    ,(    ,_+--+--,      ^)      ^\
            { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
          {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
          ( /  (    (        ,___    ^*+_+* )   <    <      \
          U _/     )    *--<  ) ^\-----++__)   )    )       )
            (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
          (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
        (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
          *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
          \             \_)^)_)) ))^^^^^^^^^^))^^^^)
          (_             ^\__^^^^^^^^^^^^))^^^^^^^)
            ^\___            ^\__^^^^^^))^^^^^^^^)\\
                  ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                    ___) >____) >___   ^\_\_\_\_\_\_\)
                    ^^^//\\_^^//\\_^       ^(\_\_\_\)
                      ^^^ ^^ ^^^ ^
```

Pipeline'ın yer aldığı kursa bağlı olarak ayrıntılarda bazı farklılıklarla karşılaşabilirsiniz.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
