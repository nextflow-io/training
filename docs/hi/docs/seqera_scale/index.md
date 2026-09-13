---
title: Scale with Seqera
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Sign up for Seqera Platform and explore the Community Showcase
    - Add a pipeline to a workspace and launch it from the web interface
    - Authenticate and launch pipelines from the command line with the `tw` CLI
    - Register a GitHub-hosted pipeline and launch it both ways
  audience_prerequisites:
    - "**Audience:** This course is designed for learners who want to run Nextflow pipelines at scale using Seqera Platform."
    - "**Skills:** Familiarity with running nf-core pipelines from the command line is assumed."
    - "**Courses:** Must have completed [Nextflow Run](../nextflow_run/index.md) and [Use nf-core](../nfcore_use/index.md), or otherwise be comfortable running local and `nf-core/rnaseq` pipelines."
---

---
title: Seqera के साथ स्केल करें
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Seqera Platform के लिए साइन अप करें और Community Showcase को एक्सप्लोर करें
    - एक workspace में पाइपलाइन जोड़ें और वेब इंटरफ़ेस से उसे लॉन्च करें
    - `tw` CLI के साथ कमांड लाइन से authenticate करें और पाइपलाइन लॉन्च करें
    - GitHub पर होस्ट की गई पाइपलाइन को रजिस्टर करें और दोनों तरीकों से लॉन्च करें
  audience_prerequisites:
    - "**दर्शक:** यह कोर्स उन लोगों के लिए बनाया गया है जो Seqera Platform का उपयोग करके Nextflow पाइपलाइन को बड़े पैमाने पर चलाना चाहते हैं।"
    - "**कौशल:** कमांड लाइन से nf-core पाइपलाइन चलाने की जानकारी होना ज़रूरी है।"
    - "**कोर्स:** [Nextflow Run](../nextflow_run/index.md) और [Use nf-core](../nfcore_use/index.md) पूरे किए हों, या फिर लोकल और `nf-core/rnaseq` पाइपलाइन चलाने में सहज हों।"
---

# Seqera के साथ स्केल करें

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Scale with Seqera, Seqera Platform के साथ Nextflow पाइपलाइन लॉन्च करने और मॉनिटर करने का एक व्यावहारिक परिचय है।**

व्यावहारिक उदाहरणों के ज़रिए, तुम Seqera Platform तक पहुँच सेट अप करोगे, वेब इंटरफ़ेस और कमांड लाइन दोनों से एक प्रोडक्शन-स्केल पाइपलाइन लॉन्च करोगे, और अपने workspace में एक नई पाइपलाइन जोड़ोगे।

तुम Seqera Platform पर अपनी खुद की पाइपलाइन चलाने और मॉनिटर करने के लिए ज़रूरी कौशल और आत्मविश्वास लेकर जाओगे।

<!-- additional_information -->

## कोर्स का अवलोकन

यह कोर्स व्यावहारिक है, और उन पाइपलाइन पर आधारित है जो तुमने [Use nf-core](../nfcore_use/index.md) में पहले से चलाई हैं।

तुम Seqera Platform के लिए साइन अप करके और वेब इंटरफ़ेस से `nf-core/rnaseq`, एक प्रोडक्शन-स्केल पाइपलाइन, लॉन्च करके शुरुआत करोगे।
फिर तुम टर्मिनल से वही काम करने के लिए `tw` कमांड-लाइन टूल पर स्विच करोगे, और अंत में एक नई पाइपलाइन `nf-core/demo` रजिस्टर करके उसे दोनों तरीकों से लॉन्च करोगे।

### पाठ योजना

| कोर्स अध्याय                                                               | सारांश                                                                                                    | अनुमानित समय |
| -------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- | ------------ |
| [भाग 1: वेब इंटरफ़ेस से पाइपलाइन लॉन्च करें](./01_run_with_seqera.md)    | Seqera Platform की पहुँच सेट अप करें और वेब इंटरफ़ेस से एक प्रोडक्शन-स्केल पाइपलाइन लॉन्च करें         | 20 मिनट      |
| [भाग 2: कमांड लाइन से पाइपलाइन लॉन्च करें](./02_launch_from_cli.md)      | `tw` CLI को authenticate करें, एक सेव की गई पाइपलाइन लॉन्च करें, और CLI से एक नई पाइपलाइन रजिस्टर करें | 25 मिनट      |

इस कोर्स के अंत तक, तुम Seqera Platform पर Nextflow पाइपलाइन लॉन्च करने और मॉनिटर करने में सहज हो जाओगे, चाहे तुम वेब इंटरफ़ेस से काम करना पसंद करो या कमांड लाइन से।

क्या तुम कोर्स शुरू करने के लिए तैयार हो?

[सीखना शुरू करें :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
