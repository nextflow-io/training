#!/usr/bin/env python3
"""
Add AI translation notice to all translated markdown files.

This script adds:
1. An admonition to index.md files (homepage)
2. An inline translation notice after the first heading in all other files
"""

import re
from pathlib import Path

# Inline notices for each language (for non-index pages)
INLINE_NOTICES = {
    "de": '<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>',
    "es": '<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>',
    "fr": '<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>',
    "hi": '<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>',
    "it": '<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>',
    "ko": '<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>',
    "pl": '<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>',
    "pt": '<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>',
    "tr": '<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>',
}

# Admonition notices for index.md (homepage) - these are inserted as-is
INDEX_ADMONITIONS = {
    "de": """!!! note "KI-gestützte Übersetzung"

    Diese Übersetzung wurde mit künstlicher Intelligenz erstellt und von menschlichen Übersetzern überprüft.
    Wir freuen uns über Feedback und Verbesserungsvorschläge.
    Weitere Informationen findest du in unserer [Übersetzungsanleitung](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md).""",
    "es": """!!! note "Traducción asistida por IA"

    Esta traducción fue creada utilizando inteligencia artificial y revisada por traductores humanos.
    Agradecemos tus comentarios y sugerencias de mejora.
    Consulta nuestra [guía de traducción](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) para más información.""",
    "fr": """!!! note "Traduction assistée par IA"

    Cette traduction a été réalisée à l'aide de l'intelligence artificielle et révisée par des traducteurs humains.
    Vos commentaires et suggestions d'amélioration sont les bienvenus.
    Consultez notre [guide de traduction](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) pour plus d'informations.""",
    "hi": """!!! note "AI-सहायता प्राप्त अनुवाद"

    यह अनुवाद कृत्रिम बुद्धिमत्ता का उपयोग करके बनाया गया था और मानव अनुवादकों द्वारा समीक्षा की गई है।
    हम आपकी प्रतिक्रिया और सुधार के सुझावों का स्वागत करते हैं।
    अधिक जानकारी के लिए हमारी [अनुवाद मार्गदर्शिका](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) देखें।""",
    "it": """!!! note "Traduzione assistita da IA"

    Questa traduzione è stata creata utilizzando l'intelligenza artificiale e revisionata da traduttori umani.
    Apprezziamo il tuo feedback e i tuoi suggerimenti per miglioramenti.
    Consulta la nostra [guida alla traduzione](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) per maggiori informazioni.""",
    "ko": """!!! note "AI 지원 번역"

    이 번역은 인공지능을 사용하여 생성되었으며 인간 번역가가 검토했습니다.
    피드백과 개선 제안을 환영합니다.
    자세한 내용은 [번역 가이드](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)를 참조하세요.""",
    "pl": """!!! note "Tłumaczenie wspomagane przez AI"

    To tłumaczenie zostało utworzone przy użyciu sztucznej inteligencji i zweryfikowane przez ludzkich tłumaczy.
    Zachęcamy do przekazywania opinii i sugestii ulepszeń.
    Więcej informacji znajdziesz w naszym [przewodniku tłumaczenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md).""",
    "pt": """!!! note "Tradução assistida por IA"

    Esta tradução foi criada utilizando inteligência artificial e revisada por tradutores humanos.
    Agradecemos seu feedback e sugestões de melhorias.
    Consulte nosso [guia de tradução](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) para mais informações.""",
    "tr": """!!! note "Yapay Zeka Destekli Çeviri"

    Bu çeviri yapay zeka kullanılarak oluşturulmuş ve insan çevirmenler tarafından gözden geçirilmiştir.
    Geri bildirimlerinizi ve iyileştirme önerilerinizi memnuniyetle karşılıyoruz.
    Daha fazla bilgi için [çeviri kılavuzumuza](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) bakın.""",
}

# Pattern to match first H1 heading (# Title)
HEADING_PATTERN = re.compile(r"^(# .+)$", re.MULTILINE)


def add_inline_notice_to_file(filepath: Path, lang: str) -> bool:
    """Add inline translation notice after first heading in a file."""
    content = filepath.read_text(encoding="utf-8")

    # Skip if notice already exists
    if "ai-translation-notice" in content:
        return False

    notice = INLINE_NOTICES[lang]

    # Find the first H1 heading and add notice after it
    match = HEADING_PATTERN.search(content)
    if match:
        heading = match.group(1)
        # Insert notice on the line after the heading
        new_content = content.replace(
            heading + "\n",
            heading + "\n\n" + notice + "\n",
            1,  # Only replace first occurrence
        )
        filepath.write_text(new_content, encoding="utf-8")
        return True

    return False


def add_admonition_to_index(filepath: Path, lang: str) -> bool:
    """Add admonition notice to index.md file."""
    content = filepath.read_text(encoding="utf-8")

    # Skip if notice already exists
    if (
        "KI-gestützte Übersetzung" in content
        or "Traducción asistida" in content
        or "Traduction assistée" in content
    ):
        return False
    if (
        "AI-सहायता प्राप्त" in content
        or "Traduzione assistita" in content
        or "AI 지원 번역" in content
    ):
        return False
    if (
        "Tłumaczenie wspomagane" in content
        or "Tradução assistida" in content
        or "Yapay Zeka Destekli" in content
    ):
        return False

    admonition = INDEX_ADMONITIONS[lang]

    # Find the first H1 heading and add admonition after it
    match = HEADING_PATTERN.search(content)
    if match:
        heading = match.group(1)
        # Insert admonition on the line after the heading
        new_content = content.replace(
            heading + "\n",
            heading + "\n\n" + admonition + "\n",
            1,  # Only replace first occurrence
        )
        filepath.write_text(new_content, encoding="utf-8")
        return True

    return False


def main():
    base_dir = Path(__file__).parent.parent / "docs"

    languages = ["de", "es", "fr", "hi", "it", "ko", "pl", "pt", "tr"]

    total_modified = 0

    for lang in languages:
        lang_dir = base_dir / lang / "docs"
        if not lang_dir.exists():
            print(f"Warning: {lang_dir} does not exist")
            continue

        modified = 0
        for md_file in lang_dir.rglob("*.md"):
            # Skip README files
            if md_file.name == "README.md":
                continue

            # Handle main index.md separately (it gets the full admonition)
            if md_file.name == "index.md" and md_file.parent == lang_dir:
                if add_admonition_to_index(md_file, lang):
                    modified += 1
            else:
                # All other files get the inline notice
                if add_inline_notice_to_file(md_file, lang):
                    modified += 1

        print(f"{lang}: Modified {modified} files")
        total_modified += modified

    print(f"\nTotal: Modified {total_modified} files")


if __name__ == "__main__":
    main()
