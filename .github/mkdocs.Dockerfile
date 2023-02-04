FROM squidfunk/mkdocs-material

# TODO: Remove this once the PR is merged + released
RUN pip install git+https://github.com/ewels/mkdocs-enumerate-headings-plugin@master#egg=mkdocs-enumerate-headings-plugin
# RUN pip install mkdocs-enumerate-headings-plugin
