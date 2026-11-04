FROM squidfunk/mkdocs-material

RUN pip install \
    mkdocs-enumerate-headings-plugin>=0.6.0 \
    mkdocs-quiz>=1.5.0
