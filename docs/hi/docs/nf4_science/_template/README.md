# nf4_science कोर्स टेम्पलेट

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

यह `nf4_science/` के अंतर्गत एक नया डोमेन-विशिष्ट कोर्स बनाने के लिए एक सामान्य टेम्पलेट है।
यह Genomics और RNAseq कोर्स में पाए जाने वाले सामान्य पैटर्न पर आधारित है।

## उपयोग कैसे करें

1. अपना कोर्स बनाने के लिए इस डायरेक्टरी और संबंधित स्क्रिप्ट डायरेक्टरी को कॉपी करो:

   ```bash
   # डॉक्स
   cp -r docs/en/docs/nf4_science/_template docs/en/docs/nf4_science/{your_domain}

   # स्क्रिप्ट्स
   cp -r nf4-science/_template nf4-science/{your_domain}
   ```

2. वर्कफ़्लो स्क्रिप्ट का नाम बदलो:

   ```bash
   mv nf4-science/{your_domain}/_template.nf nf4-science/{your_domain}/{your_domain}.nf
   ```

3. मॉड्यूल फ़ाइलों का नाम अपने टूल्स के अनुसार बदलो:

   ```bash
   mv nf4-science/{your_domain}/modules/_tool_a.nf nf4-science/{your_domain}/modules/{tool_a}.nf
   mv nf4-science/{your_domain}/modules/_tool_b.nf nf4-science/{your_domain}/modules/{tool_b}.nf
   ```

4. डॉक्स और स्क्रिप्ट्स में सभी `{PLACEHOLDER}` वैल्यू खोजो और उन्हें अपनी सामग्री से बदलो।

5. `mkdocs.yml` nav में जोड़ो (नीचे स्निपेट देखो)।

6. `data/` डायरेक्टरी को टेस्ट डेटासेट्स से भरो।

7. `solutions/` डायरेक्टरी को प्रत्येक भाग के लिए काम करने वाले कोड के साथ तैयार करो।

## mkdocs.yml nav स्निपेट

```yaml
- { Your Domain }:
    - nf4_science/{your_domain}/index.md
    - nf4_science/{your_domain}/00_orientation.md
    - nf4_science/{your_domain}/01_method.md
    - nf4_science/{your_domain}/02_single_sample.md
    - nf4_science/{your_domain}/03_multi_sample.md
    - nf4_science/{your_domain}/survey.md
    - nf4_science/{your_domain}/next_steps.md
```

## प्लेसहोल्डर संदर्भ

### कोर्स-स्तर के प्लेसहोल्डर

| प्लेसहोल्डर                  | विवरण                               | उदाहरण (Genomics)                          |
| ---------------------------- | ----------------------------------- | ------------------------------------------ |
| `{DOMAIN}`                   | डोमेन नाम (टाइटल केस)               | Genomics                                   |
| `{DOMAIN_DIR}`               | डायरेक्टरी नाम (लोअरकेस)            | genomics                                   |
| `{METHOD}`                   | विश्लेषण विधि का नाम                | variant calling                            |
| `{METHOD_SHORT_DESCRIPTION}` | एक-पंक्ति विधि विवरण                | variant calling with GATK                  |
| `{ACCESSORY_FILES}`          | उपयोग की गई सहायक फ़ाइलों के प्रकार | index files and reference genome resources |

### टूल प्लेसहोल्डर

| प्लेसहोल्डर                           | विवरण                          | उदाहरण (Genomics)                                   |
| ------------------------------------- | ------------------------------ | --------------------------------------------------- |
| `{TOOL_A}` / `{TOOL_B}`               | टूल डिस्प्ले नाम               | Samtools / GATK                                     |
| `{TOOL_A_URL}` / `{TOOL_B_URL}`       | टूल डॉक्यूमेंटेशन URLs         | https://www.htslib.org/                             |
| `{TOOL_A_CONTAINER_URI}`              | पूरा कंटेनर URI                | community.wave.seqera.io/library/samtools:1.20--... |
| `{TOOL_A_MODULE}` / `{TOOL_B_MODULE}` | मॉड्यूल फ़ाइलनाम (.nf के बिना) | samtools_index / gatk_haplotypecaller               |
| `{TOOL_A_PROCESS_NAME}`               | UPPERCASE प्रोसेस नाम          | SAMTOOLS_INDEX                                      |
| `{TOOL_A_COMMAND}`                    | चलाने के लिए शेल कमांड         | samtools index '<input_bam>'                        |

### इनपुट/आउटपुट प्लेसहोल्डर

| प्लेसहोल्डर            | विवरण                        | उदाहरण (Genomics)    |
| ---------------------- | ---------------------------- | -------------------- |
| `{PRIMARY_INPUT_TYPE}` | प्राथमिक इनपुट फ़ाइल प्रकार  | BAM file             |
| `{PRIMARY_PARAM_NAME}` | Nextflow पैरामीटर नाम        | reads_bam            |
| `{TEST_INPUT_PATH}`    | data/ में टेस्ट इनपुट का पाथ | bam/reads_mother.bam |
| `{FINAL_OUTPUT_TYPE}`  | अंतिम पाइपलाइन आउटपुट प्रकार | joint-called VCFs    |

### सामग्री प्लेसहोल्डर

| प्लेसहोल्डर             | विवरण                               |
| ----------------------- | ----------------------------------- |
| `{PART2_SUMMARY}`       | भाग 2 के दायरे का एक-पंक्ति सारांश  |
| `{AGGREGATION_SUMMARY}` | एकत्रीकरण चरण का एक-पंक्ति सारांश   |
| `{PROCESS_LIST}`        | प्रोसेस नामों की कॉमा-सेपरेटेड सूची |
| `{TYPEFORM_ID}`         | सर्वे के लिए Typeform embed ID      |

## शैक्षणिक संरचना

टेम्पलेट तीन-भाग की संरचना का पालन करता है:

1. **भाग 1 (01_method.md)**: पद्धति को समझने के लिए Docker कंटेनरों में मैनुअल टेस्टिंग
2. **भाग 2 (02_single_sample.md)**: कमांड्स को Nextflow में रैप करना; सिंगल-सैंपल, फिर बैच
3. **भाग 3 (03_multi_sample.md)**: चैनल ऑपरेटर्स का उपयोग करके मल्टी-सैंपल एकत्रीकरण जोड़ना

### मुख्य परंपराएं

- कोड परिवर्तन दिखाने के लिए `hl_lines` के साथ **Before/After टैब्ड कोड ब्लॉक** का उपयोग करो
- कोलैप्सिबल अपेक्षित आउटपुट के लिए `??? success "Command output"` का उपयोग करो
- कोलैप्सिबल डायरेक्टरी ट्री के लिए `??? abstract "Directory contents"` का उपयोग करो
- प्रत्येक प्रमुख सेक्शन को **सारांश** और **आगे क्या है?** उपसेक्शन के साथ समाप्त करो
- प्रमुख क्रमांकित सेक्शन को अलग करने के लिए `---` हॉरिजॉन्टल रूल्स का उपयोग करो
- सीखने वालों के लिए भरने के लिए स्केलेटन फ़ाइलें (वर्कफ़्लो + मॉड्यूल) प्रदान करो
- समाधानों को भाग के अनुसार व्यवस्थित करो (`solutions/part2/`, `solutions/part3/`)

## डायरेक्टरी संरचना

```
docs/en/docs/nf4_science/{domain}/
├── index.md                    # Course overview with frontmatter
├── 00_orientation.md           # Environment setup
├── 01_method.md                # Manual testing in containers
├── 02_single_sample.md         # Single-sample Nextflow implementation
├── 03_multi_sample.md          # Multi-sample aggregation
├── survey.md                   # Typeform feedback survey
├── next_steps.md               # Course summary and suggestions
└── img/                        # Diagrams (.excalidraw.svg, .png)

nf4-science/{domain}/
├── {domain}.nf                 # Skeleton workflow file
├── nextflow.config             # Minimal config (docker.enabled = true)
├── data/                       # Test datasets and resources
│   └── samplesheet.csv         # Sample metadata
├── modules/                    # Skeleton module files
│   ├── {tool_a}.nf
│   └── {tool_b}.nf
└── solutions/
    ├── part2/                  # Complete Part 2 solution
    │   ├── {domain}-2.nf
    │   ├── nextflow.config
    │   └── modules/
    └── part3/                  # Complete Part 3 solution
        ├── {domain}-3.nf
        ├── nextflow.config
        └── modules/
```
