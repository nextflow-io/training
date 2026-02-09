---
title: Hello boru hattı
description: Hello boru hattının ne yaptığı ve nasıl yapılandırıldığının özeti.
hide:
  - toc
  - footer
---

# Hello boru hattı

Eğitim kurslarımızın çoğu, Nextflow kavramlarını ve mekanizmalarını göstermek için basit, alana özgü olmayan bir boru hattı kullanır.
Hello Nextflow kursu, bu boru hattının her tasarım ve uygulama kararını açıklayan adım adım bir şekilde nasıl geliştirileceğini gösterir.
Diğer eğitimler bu boru hattını veya parçalarını başlangıç noktası olarak kullanır.

Bu sayfa, Hello Nextflow kursunun tamamlanmasıyla boru hattının ulaştığı durumu özetlemektedir.

### Özet açıklama

Hello iş akışı, selamlamaları içeren bir CSV dosyası alır, bunları ayrı dosyalara yazar, her birini büyük harfe dönüştürür, tekrar bir araya toplar ve selamlamaları söyleyen eğlenceli bir karakterin ASCII resmini içeren tek bir metin dosyası çıktısı verir.

### İş akışı adımları (süreçler)

Dört adım, ayrı modül dosyalarında saklanan Nextflow süreçleri (`sayHello`, `convertToUpper`, `collectGreetings` ve `cowpy`) olarak uygulanmıştır.

1. **`sayHello`:** Her selamlamayı kendi çıktı dosyasına yazar (örn., "Hello-output.txt")
2. **`convertToUpper`:** Her selamlamayı büyük harfe dönüştürür (örn., "HELLO")
3. **`collectGreetings`:** Tüm büyük harfli selamlamaları tek bir toplu dosyada toplar
4. **`cowpy`:** `cowpy` aracını kullanarak ASCII sanatı oluşturur

### Diyagram

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### Sonuçlar

Sonuçlar `results/` adlı bir dizinde yayınlanır ve boru hattının nihai çıktısı (varsayılan parametrelerle çalıştırıldığında) büyük harfli selamlamaları söyleyen bir hindinin ASCII sanatını içeren düz bir metin dosyasıdır.

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

Boru hattının yer aldığı kursa bağlı olarak ayrıntılarda bazı farklılıklarla karşılaşabilirsiniz.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
